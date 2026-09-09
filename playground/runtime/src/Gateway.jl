const MAX_GATEWAY_REQUEST = 256 * 1024
const MAX_UI_RESPONSE = 32 * 1024 * 1024
const MAX_UI_REQUEST = 32 * 1024 * 1024

# Bound buffered HTTP assets independently of untrusted Content-Length.
mutable struct BoundedResponseBody <: IO
    buffer::IOBuffer
    limit::Int
end
Base.isopen(body::BoundedResponseBody) = isopen(body.buffer)
function Base.unsafe_write(body::BoundedResponseBody, pointer::Ptr{UInt8}, count::UInt)
    count <= body.limit - position(body.buffer) ||
        throw(ArgumentError("UI response exceeds the configured byte bound"))
    return Base.unsafe_write(body.buffer, pointer, count)
end

function connection_peer(stream)
    address = HTTP.peeraddr(stream)
    isnothing(address) && throw(AccessDenied(401, "Connection identity unavailable"))
    octets = address.ip
    if length(octets) == 4
        return join(octets, '.')
    end
    length(octets) == 16 || throw(AccessDenied(401, "Connection identity unavailable"))
    groups = ntuple(i -> (UInt16(octets[2i - 1]) << 8) | UInt16(octets[2i]), 8)
    return string(IPv6(groups...))
end

function gateway_response(stream, status::Int, body; content_type="application/json", headers=Pair{String,String}[])
    payload = body isa AbstractString ? Vector{UInt8}(codeunits(body)) : body
    HTTP.setstatus(stream, status)
    for (key,value) in headers
        HTTP.setheader(stream, key, value)
    end
    HTTP.setheader(stream, "Content-Type", content_type)
    HTTP.setheader(stream, "Cache-Control", "no-store")
    HTTP.setheader(stream, "X-Content-Type-Options", "nosniff")
    HTTP.setheader(stream, "Content-Length", string(length(payload)))
    HTTP.startwrite(stream)
    stream.message.method == "HEAD" || write(stream, payload)
    closewrite(stream)
    return nothing
end

run_payload(run::RunRecord) = (
    id=string(run.id), application=run.application, version=string(run.version),
    state=String(run.state), created_at=string(run.created_at),
    updated_at=string(run.updated_at), reason=run.reason,
)
application_payload(app::ApplicationDefinition) = (
    id=app.id, title=app.title, kind=String(app.kind), entrypoint=app.entrypoint,
    version=string(app.version), visibility=String(app.visibility),
    entry_surface=String(app.entry_surface),
    requirements=[(role=r.role, profiles=r.profiles, required=r.required) for r in app.requirements],
)

function validate_gateway_json(value, depth=0)
    depth <= 32 || throw(ArgumentError("request nesting exceeds its bound"))
    if value isa JSON3.Object
        length(value) == length(Set(keys(value))) || throw(ArgumentError("duplicate request fields"))
        foreach(item -> validate_gateway_json(item, depth + 1), values(value))
    elseif value isa JSON3.Array
        foreach(item -> validate_gateway_json(item, depth + 1), value)
    end
    return nothing
end

function read_gateway_json(stream)
    bytes = read(stream, MAX_GATEWAY_REQUEST + 1)
    length(bytes) <= MAX_GATEWAY_REQUEST || throw(AccessDenied(413, "Request is too large"))
    try
        object = JSON3.read(bytes)
        object isa JSON3.Object || throw(ArgumentError("request must be an object"))
        validate_gateway_json(object)
        return JSON3.read(bytes, Dict{String,Any})
    catch
        throw(AccessDenied(400, "Expected a JSON object"))
    end
end

function requested_uuid(value)
    value isa AbstractString || throw(AccessDenied(400, "Invalid request identity"))
    id = tryparse(UUID, value)
    isnothing(id) && throw(AccessDenied(400, "Invalid request identity"))
    return id
end

function proxy_handle(supervisor::UIHostSupervisor, principal::Principal, id::UUID)
    lock(supervisor.lock) do
        record = get_run(supervisor.store, principal, id)
        record.state == :running || throw(AccessDenied(503, "UI host is unavailable"))
        handle = get(supervisor.handles, id, nothing)
        !isnothing(handle) && !isnothing(handle.port) && Base.process_running(handle.process) ||
            throw(AccessDenied(503, "UI host is unavailable"))
        handle.last_seen = time_ns()
        return handle
    end
end

function forward_ui_http(stream, handle::UIHostHandle, target::String)
    method = stream.message.method
    method in ("GET", "HEAD", "POST", "PUT", "PATCH", "DELETE") ||
        throw(AccessDenied(405, "Method is not supported"))
    body = read(stream, MAX_UI_REQUEST + 1)
    length(body) <= MAX_UI_REQUEST || throw(AccessDenied(413, "Request is too large"))
    allowed_request = ("content-type", "accept", "accept-encoding", "if-none-match",
        "if-modified-since", "range", "x-lcm-original-name", "x-lcm-upload-generation", "x-lcm-upload")
    headers = [String(k)=>String(v) for (k,v) in stream.message.headers if lowercase(k) in allowed_request]
    push!(headers, "X-LCM-Host-Key" => handle.key)
    output = BoundedResponseBody(IOBuffer(), MAX_UI_RESPONSE)
    response = HTTP.request(method, "http://127.0.0.1:$(handle.port)" * target, headers, body;
        response_stream=output, request_timeout=15, connect_timeout=2, proxy=nothing,
        redirect=false, retry=false, status_exception=false, decompress=false)
    allowed_response = ("content-type", "content-encoding", "etag", "last-modified",
        "content-range", "accept-ranges")
    forwarded = [String(k)=>String(v) for (k,v) in response.headers if lowercase(k) in allowed_response]
    gateway_response(stream, response.status, take!(output.buffer);
        content_type=HTTP.header(response, "Content-Type", "application/octet-stream"), headers=forwarded)
end

function relay_ui_socket(stream, supervisor::UIHostSupervisor, principal::Principal,
        handle::UIHostHandle, target::String)
    lock(supervisor.lock) do
        proxy_handle(supervisor, principal, handle.id)
        handle.connections += 1
    end
    try
        relay_authorized_ui_socket(stream, supervisor, principal, handle, target)
    finally
        lock(supervisor.lock) do
            handle.connections -= 1
            handle.last_seen = time_ns()
        end
    end
end

function relay_authorized_ui_socket(stream, supervisor::UIHostSupervisor, principal::Principal,
        handle::UIHostHandle, target::String)
    # Authenticate before either upgrade. Upstream gets only its private key,
    # never the browser's proxy credential, cookies or asserted identity.
    HTTP.WebSockets.open("ws://127.0.0.1:$(handle.port)" * target;
            headers=["X-LCM-Host-Key"=>handle.key], cookies=false, proxy=nothing,
            redirect=false, request_timeout=3, maxframesize=MAX_UI_RESPONSE) do upstream
        HTTP.WebSockets.upgrade(stream;
                check_origin=(_...) -> true, maxframesize=MAX_UI_RESPONSE) do downstream
            @sync begin
                @async try
                    for message in downstream
                        proxy_handle(supervisor, principal, handle.id)
                        HTTP.WebSockets.send(upstream, message)
                    end
                catch
                    # Close both halves; protocol/backend details stay private.
                finally
                    close(upstream)
                end
                try
                    for message in upstream
                        proxy_handle(supervisor, principal, handle.id)
                        HTTP.WebSockets.send(downstream, message)
                    end
                catch
                finally
                    close(downstream)
                end
            end
        end
    end
end

function gateway_request(supervisor::UIHostSupervisor, policy::AbstractIdentityPolicy, stream;
        site::Union{Nothing,PublishedSite}=nothing, control::Union{Nothing,ControlService}=nothing)
    request = stream.message
    target = String(request.target)
    startswith(target, "/") && !startswith(target, "//") &&
        !occursin(r"[\\\x00\r\n]", target) ||
        throw(AccessDenied(400, "Invalid request target"))
    path = String(URIs.URI(target).path)
    if site !== nothing && request.method in ("GET", "HEAD") && serve_published(stream, site, path)
        return nothing
    end
    asset = match(r"^/runtime/assets/([a-z.-]+)$", path)
    if asset !== nothing && request.method in ("GET", "HEAD")
        haskey(RUNTIME_ASSETS, asset[1]) || throw(AccessDenied(404, "Asset not found"))
        file, mime = RUNTIME_ASSETS[asset[1]]
        return gateway_response(stream, 200, read(file); content_type=mime)
    end
    if path == "/health" && request.method == "GET"
        return gateway_response(stream, 200, "ok"; content_type="text/plain")
    elseif path == "/runtime/api/capabilities" && request.method == "GET"
        return gateway_response(stream, 200, JSON3.write((; schema_version=1, ui_hosts=true,
            control_capabilities(control)...,
            applications=sort!(collect(keys(supervisor.registry.applications))))))
    elseif path == "/runtime/api/applications" && request.method == "GET"
        descriptions = sort!([app for app in values(supervisor.registry.definitions)
            if app.visibility == :public]; by=app -> app.id)
        return gateway_response(stream, 200, JSON3.write(application_payload.(descriptions)))
    end
    principal = authorize_request(policy, request.headers, connection_peer(stream);
        method=request.method, websocket=HTTP.WebSockets.isupgrade(request))
    if policy isa LocalIdentity
        # A loopback peer alone is insufficient against browser DNS rebinding.
        # Local development accepts only its explicitly configured authority.
        origin = URIs.URI(policy.origin)
        authority = origin.host * (isempty(origin.port) ? "" : ":" * origin.port)
        lowercase(only_header(request.headers, "Host")) == lowercase(authority) ||
            throw(AccessDenied(403, "Local development authority does not match"))
    end
    if path == "/runtime/control" && request.method == "GET"
        return control_surface(stream, supervisor, principal)
    end
    terminal_gateway_request(control,principal,stream,path) && return nothing
    control_request(control, supervisor, principal, stream, path) && return nothing
    if path == "/runtime/api/runs"
        if request.method == "GET"
            return gateway_response(stream, 200,
                JSON3.write(run_payload.(list_runs(supervisor.store, principal))))
        elseif request.method == "POST"
            data = read_gateway_json(stream)
            Set(keys(data)) == Set(("application", "request_id")) ||
                throw(AccessDenied(400, "Expected application and request_id only"))
            data["application"] isa AbstractString ||
                throw(AccessDenied(400, "Invalid application identity"))
            run = start_ui!(supervisor, principal, data["application"];
                request_id=requested_uuid(data["request_id"]))
            return gateway_response(stream, 202, JSON3.write(run_payload(run)))
        end
        throw(AccessDenied(405, "Method is not supported"))
    end
    surface = match(r"^/runtime/runs/([a-f0-9-]{36})$", path)
    if surface !== nothing && request.method == "GET"
        run = get_run(supervisor.store, principal, requested_uuid(surface[1]))
        return run_surface(stream, supervisor, run)
    end
    run_route = match(r"^/runtime/api/runs/([a-f0-9-]{36})$", path)
    if !isnothing(run_route)
        id = requested_uuid(run_route[1])
        run = if request.method == "GET"
            get_run(supervisor.store, principal, id)
        elseif request.method == "DELETE"
            stop_ui!(supervisor, principal, id)
        else
            throw(AccessDenied(405, "Method is not supported"))
        end
        return gateway_response(stream, 200, JSON3.write(run_payload(run)))
    end
    ui_route = match(r"^/applications/runs/([a-f0-9-]{36})/", path)
    if !isnothing(ui_route)
        id = requested_uuid(ui_route[1])
        handle = try
            proxy_handle(supervisor, principal, id)
        catch error
            if error isa AccessDenied && error.status == 503 && request.method == "GET" &&
                    occursin("text/html", HTTP.header(request, "Accept", ""))
                run = get_run(supervisor.store, principal, id)
                return run_surface(stream, supervisor, run; status=503, automatic=false)
            end
            rethrow()
        end
        upstream_target = "/" * target[(ncodeunits(ui_route.match) + 1):end]
        # No decoding/normalization may turn a namespaced target into an
        # authority or a path outside this particular upstream process.
        occursin(r"(?i)%2f|%5c|%2e|/\.\.?(/|$)", upstream_target) &&
            throw(AccessDenied(400, "Ambiguous UI path"))
        return HTTP.WebSockets.isupgrade(request) ?
            relay_ui_socket(stream, supervisor, principal, handle, upstream_target) :
            forward_ui_http(stream, handle, upstream_target)
    end
    throw(AccessDenied(404, "Route not found"))
end

function compile_gateway_paths(supervisor, policy, site, control)
    # Compile this configured HTTP/1 entry before opening the listener. The
    # latest-world boundary below intentionally does not compile it transitively;
    # charging that cold specialization to the first browser request can exceed
    # its deadline. This requests code generation only, without invoking routes,
    # starting control services, allocating runs or reading private state.
    keywords = NamedTuple{(:site, :control),Tuple{typeof(site),typeof(control)}}
    stream_type = HTTP.Stream{false,HTTP.Request{HTTP.EmptyBody}}
    precompile(Core.kwcall, (keywords, typeof(gateway_request),
        typeof(supervisor), typeof(policy), stream_type))
    if control !== nothing
        precompile(relay_terminal_socket!, (typeof(control.terminals), TerminalAttachment))
        precompile(terminal_browser_request, (typeof(control.terminals), TerminalAttachment, Dict{String,Any}))
    end
    return nothing
end

"""
    start_gateway(supervisor, policy; host="127.0.0.1", port=0) -> HTTP.Server

Expose the owned runtime API and run-prefixed HTTP/WebSocket proxy on a private
loopback listener. HTTP identity, mutation-origin and run ownership checks happen
before forwarding; public catalogue access never launches an application.
"""
function start_gateway(supervisor::UIHostSupervisor, policy::AbstractIdentityPolicy;
        host::AbstractString="127.0.0.1", port::Integer=0,
        site::Union{Nothing,PublishedSite}=nothing, control::Union{Nothing,ControlService}=nothing)
    host in ("127.0.0.1", "::1") || throw(ArgumentError("gateway listener must be private loopback"))
    compile_gateway_paths(supervisor, policy, site, control)
    # Stream handlers enforce body bounds through the capped reads above;
    # HTTP.listen! does not accept the buffered serve! max_body_bytes option.
    return HTTP.listen!(host, port; max_header_bytes=32 * 1024, read_header_timeout=10,
            read_timeout=20, write_timeout=20) do stream
        try
            Base.invokelatest(gateway_request, supervisor, policy, stream; site, control)
        catch error
            status, reason = if error isa AccessDenied
                error.status, error.reason
            elseif error isa CapacityUnavailable
                429, "Application capacity is unavailable"
            elseif error isa ArgumentError
                400, "Invalid runtime request"
            else
                503, "Runtime request could not be completed"
            end
            gateway_response(stream, status, JSON3.write((error=reason,)))
        end
    end
end
