const UPLOAD_ROUTE_PATTERN = r"^/uploads/[0-9a-f-]{36}$"
const UPLOAD_KEY_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9._-]*$"
const UPLOAD_RETENTIONS = (:persistent, :session)
const UPLOAD_STATES = (
    :empty, :selected, :uploading, :validating, :ready, :removing, :failed,
)

const UPLOAD_STYLES_PATH = normpath(joinpath(@__DIR__, "..", "assets", "upload.css"))
include_dependency(UPLOAD_STYLES_PATH)
const UPLOAD_STYLES = read(UPLOAD_STYLES_PATH, String)

abstract type AbstractUploadStore end

"Store one logical upload per target beneath a server-owned directory."
struct DirectoryUploadStore <: AbstractUploadStore
    root::String
    staging::String
end

function DirectoryUploadStore(root::AbstractString)
    resolved_root = abspath(expanduser(string(root)))
    store = DirectoryUploadStore(
        resolved_root,
        joinpath(resolved_root, ".lcm-upload-staging")
    )
    mkpath(store.staging)
    sweep_upload_staging!(store)
    return store
end

"Validation and retention rules applied independently of browser hints."
struct UploadPolicy
    max_bytes::Int
    extensions::Set{String}
    media_types::Set{String}
    retention::Symbol
    allow_empty::Bool
end

function normalized_upload_extension(extension::AbstractString)
    value = lowercase(strip(string(extension)))
    isempty(value) && throw(ArgumentError("upload extensions cannot be empty"))
    startswith(value, '.') || (value = ".$value")
    occursin(r"^\.[a-z0-9][a-z0-9._+-]*$", value) || throw(ArgumentError(
        "invalid upload extension: $extension"
    ))
    return value
end

function normalized_media_type(media_type::AbstractString)
    value = lowercase(strip(first(split(string(media_type), ';'; limit=2))))
    occursin(r"^[a-z0-9][a-z0-9!#$&^_.+-]*/[a-z0-9][a-z0-9!#$&^_.+-]*$", value) ||
        throw(ArgumentError("invalid upload media type: $media_type"))
    return value
end

function UploadPolicy(;
        max_bytes::Integer=16 * 1024 * 1024,
        extensions=String[],
        media_types=String[],
        retention::Symbol=:persistent,
        allow_empty::Bool=false
    )
    max_bytes > 0 || throw(ArgumentError("upload max_bytes must be positive"))
    retention in UPLOAD_RETENTIONS || throw(ArgumentError(
        "upload retention must be :persistent or :session"
    ))
    return UploadPolicy(
        Int(max_bytes),
        Set(normalized_upload_extension.(string.(collect(extensions)))),
        Set(normalized_media_type.(string.(collect(media_types)))),
        retention,
        allow_empty
    )
end

_accept_upload(_path) = nothing

"A server-owned destination and its validation policy."
struct UploadTarget{Store<:AbstractUploadStore,Validator}
    store::Store
    key::String
    policy::UploadPolicy
    validate::Validator
end

function UploadTarget(
        store::AbstractUploadStore,
        key::Union{Symbol,AbstractString};
        max_bytes::Integer=16 * 1024 * 1024,
        extensions=String[],
        media_types=String[],
        retention::Symbol=:persistent,
        allow_empty::Bool=false,
        validate=_accept_upload
    )
    resolved_key = string(key)
    occursin(UPLOAD_KEY_PATTERN, resolved_key) || throw(ArgumentError(
        "upload key must be one filename-like identifier without path segments"
    ))
    applicable(validate, "candidate-path") || throw(ArgumentError(
        "upload validator must accept the staged file path"
    ))
    policy = UploadPolicy(;
        max_bytes,
        extensions,
        media_types,
        retention,
        allow_empty
    )
    return UploadTarget(store, resolved_key, policy, validate)
end

"Metadata for the file currently committed to an upload slot."
struct UploadedFile
    key::String
    original_name::String
    media_type::String
    bytes::Int
    sha256::String
    committed_at::DateTime
end

struct UploadRejected <: Exception
    status::Int
    message::String
end

Base.showerror(io::IO, error::UploadRejected) = print(io, error.message)

upload_path(target::UploadTarget{<:DirectoryUploadStore}) =
    joinpath(target.store.root, target.key)

function clean_original_name(name::AbstractString)
    normalized = replace(string(name), '\\' => '/')
    cleaned = strip(replace(basename(normalized), r"[\x00-\x1f\x7f]" => ""))
    isempty(cleaned) && throw(UploadRejected(400, "The selected file has no usable name"))
    ncodeunits(cleaned) <= 255 || throw(UploadRejected(400, "The selected filename is too long"))
    return cleaned
end

function check_upload!(target::UploadTarget, name, media_type, byte_count)
    policy = target.policy
    byte_count <= policy.max_bytes || throw(UploadRejected(
        413,
        "File exceeds the $(policy.max_bytes)-byte limit"
    ))
    (policy.allow_empty || byte_count > 0) || throw(UploadRejected(
        400,
        "Empty files are not accepted"
    ))
    if !isempty(policy.extensions)
        extension = lowercase(splitext(name)[2])
        extension_label = isempty(extension) ? "<none>" : extension
        extension in policy.extensions || throw(UploadRejected(
            415,
            "File extension $extension_label is not accepted"
        ))
    end
    if !isempty(policy.media_types)
        media_type in policy.media_types || throw(UploadRejected(
            415,
            "Media type $media_type is not accepted"
        ))
    end
    return nothing
end

function sweep_upload_staging!(
        store::DirectoryUploadStore;
        older_than_seconds::Real=24 * 60 * 60
    )
    older_than_seconds >= 0 || throw(ArgumentError(
        "staging-file lifetime must be nonnegative"
    ))
    isdir(store.staging) || return 0
    removed = 0
    cutoff = time() - Float64(older_than_seconds)
    for name in readdir(store.staging)
        path = joinpath(store.staging, name)
        isfile(path) || continue
        stat(path).mtime <= cutoff || continue
        rm(path; force=true)
        removed += 1
    end
    return removed
end

function validate_staged_upload!(target::UploadTarget, path::AbstractString)
    try
        result = target.validate(path)
        isnothing(result) || result === true || throw(ArgumentError(
            "upload validator must return nothing or true"
        ))
    catch error
        error isa UploadRejected && rethrow()
        throw(UploadRejected(
            422,
            "File validation failed: $(sprint(showerror, error))"
        ))
    end
    return nothing
end

function store_upload!(
        target::UploadTarget{<:DirectoryUploadStore},
        contents::AbstractVector{UInt8};
        original_name::AbstractString,
        media_type::AbstractString="application/octet-stream"
    )
    name = clean_original_name(original_name)
    normalized_type = try
        normalized_media_type(media_type)
    catch error
        error isa ArgumentError || rethrow()
        throw(UploadRejected(415, "The selected file has an invalid media type"))
    end
    check_upload!(target, name, normalized_type, length(contents))

    mkpath(target.store.root)
    mkpath(target.store.staging)
    temporary_path, io = mktemp(target.store.staging; cleanup=false)
    committed = false
    try
        write(io, contents)
        flush(io)
        close(io)
        validate_staged_upload!(target, temporary_path)
        digest = bytes2hex(SHA.sha256(contents))
        Base.Filesystem.rename(temporary_path, upload_path(target))
        committed = true
        return UploadedFile(
            target.key,
            name,
            normalized_type,
            length(contents),
            digest,
            now(UTC)
        )
    finally
        isopen(io) && close(io)
        !committed && ispath(temporary_path) && rm(temporary_path; force=true)
    end
end

function delete_upload!(target::UploadTarget{<:DirectoryUploadStore})
    path = upload_path(target)
    ispath(path) && rm(path; force=true)
    return nothing
end

mutable struct UploadRegistry
    entries::Dict{String,Any}
    lock::ReentrantLock
end

UploadRegistry() = UploadRegistry(Dict{String,Any}(), ReentrantLock())

const DEFAULT_UPLOAD_REGISTRY = UploadRegistry()
default_upload_registry() = DEFAULT_UPLOAD_REGISTRY

"A reusable, single-file upload slot rendered by Bonito."
mutable struct UploadField{Target<:UploadTarget}
    target::Target
    registry::UploadRegistry
    value::Observable{Union{Nothing,UploadedFile}}
    status::Observable{Symbol}
    progress::Observable{Float64}
    error::Observable{String}
    label::String
    hint::String
    state_lock::ReentrantLock
    generation::Int
    phase::Int
end

function UploadField(
        target::UploadTarget;
        label::AbstractString="File",
        hint::AbstractString="Select a file to upload.",
        initial::Union{Nothing,UploadedFile}=nothing,
        registry::UploadRegistry=default_upload_registry()
    )
    initial_status = isnothing(initial) ? :empty : :ready
    return UploadField(
        target,
        registry,
        Observable{Union{Nothing,UploadedFile}}(initial),
        Observable(initial_status),
        Observable(isnothing(initial) ? 0.0 : 1.0),
        Observable(""),
        string(label),
        string(hint),
        ReentrantLock(),
        0,
        0
    )
end

function apply_upload_state!(
        field::UploadField,
        generation::Integer,
        phase::Integer,
        status::Symbol;
        progress::Real=field.progress[],
        error::AbstractString=""
    )
    status in UPLOAD_STATES || return false
    return lock(field.state_lock) do
        generation < field.generation && return false
        if generation > field.generation
            field.generation = Int(generation)
            field.phase = -1
        end
        phase < field.phase && return false
        field.phase = Int(phase)
        field.status[] = status
        field.progress[] = clamp(Float64(progress), 0.0, 1.0)
        field.error[] = string(error)
        return true
    end
end

mutable struct UploadRegistration
    target::UploadTarget
    field::UploadField
    latest_generation::Int
    lock::ReentrantLock
end

struct UploadGateway
    registry::UploadRegistry
end

function register_upload!(registry::UploadRegistry, target::UploadTarget, field::UploadField)
    token = string(uuid4())
    registration = UploadRegistration(
        target,
        field,
        field.generation,
        ReentrantLock()
    )
    lock(registry.lock) do
        registry.entries[token] = registration
    end
    return token
end

function unregister_upload!(registry::UploadRegistry, token::AbstractString)
    return lock(registry.lock) do
        pop!(registry.entries, string(token), nothing)
    end
end

function registered_upload(registry::UploadRegistry, token::AbstractString)
    return lock(registry.lock) do
        get(registry.entries, string(token), nothing)
    end
end

function upload_token_from_target(target::AbstractString)
    path = first(split(string(target), '?'; limit=2))
    matched = match(r"^/uploads/([0-9a-f-]{36})$", path)
    return isnothing(matched) ? nothing : matched.captures[1]
end

function upload_request_bytes(request)
    body = request.body
    body isa Bonito.HTTP.EmptyBody && return UInt8[]
    if body isa Bonito.HTTP.BytesBody
        body.next_index > length(body.data) && return UInt8[]
        return @view body.data[body.next_index:end]
    end
    throw(UploadRejected(400, "Unsupported upload request body"))
end

function upload_generation(request)
    value = Bonito.HTTP.header(request, "X-LCM-Upload-Generation", "")
    generation = tryparse(Int, value)
    isnothing(generation) && throw(UploadRejected(400, "Missing upload generation"))
    generation > 0 || throw(UploadRejected(400, "Invalid upload generation"))
    return generation
end

function json_upload_response(status::Integer, value)
    return Bonito.HTTP.Response(
        status,
        [
            "Content-Type" => "application/json; charset=utf-8",
            "Cache-Control" => "no-store",
            "X-Content-Type-Options" => "nosniff",
        ];
        body=JSON3.write(value)
    )
end

uploaded_file_payload(file::UploadedFile) = (
    key=file.key,
    original_name=file.original_name,
    media_type=file.media_type,
    bytes=file.bytes,
    sha256=file.sha256,
    committed_at=string(file.committed_at),
)

function perform_upload!(registration::UploadRegistration, request, generation::Int)
    field = registration.field
    apply_upload_state!(field, generation, 3, :validating; progress=1.0)
    encoded_name = Bonito.HTTP.header(request, "X-LCM-Original-Name", "")
    isempty(encoded_name) && throw(UploadRejected(400, "Missing original filename"))
    name = try
        URIs.unescapeuri(encoded_name)
    catch
        throw(UploadRejected(400, "Invalid original filename encoding"))
    end
    media_type = Bonito.HTTP.header(
        request,
        "Content-Type",
        "application/octet-stream"
    )
    file = store_upload!(
        registration.target,
        upload_request_bytes(request);
        original_name=name,
        media_type
    )
    field.value[] = file
    apply_upload_state!(field, generation, 4, :ready; progress=1.0)
    return json_upload_response(201, (status="ready", file=uploaded_file_payload(file)))
end

function perform_upload_removal!(registration::UploadRegistration, generation::Int)
    field = registration.field
    apply_upload_state!(field, generation, 2, :removing; progress=0.0)
    delete_upload!(registration.target)
    field.value[] = nothing
    apply_upload_state!(field, generation, 4, :empty; progress=0.0)
    return json_upload_response(200, (status="empty",))
end

function perform_upload_download(registration::UploadRegistration)
    file = registration.field.value[]
    path = upload_path(registration.target)
    (isnothing(file) || !isfile(path)) && return json_upload_response(
        404,
        (error="No file is committed to this upload slot",)
    )
    return Bonito.HTTP.Response(
        200,
        [
            "Content-Type" => file.media_type,
            "Content-Length" => string(file.bytes),
            "Cache-Control" => "no-store",
            "X-Content-Type-Options" => "nosniff",
        ];
        body=read(path)
    )
end

function Bonito.HTTPServer.apply_handler(gateway::UploadGateway, context)
    request = context.request
    request.method in ("GET", "POST", "DELETE") || return Bonito.HTTP.Response(
        405,
        ["Allow" => "GET, POST, DELETE"];
        body="Method not allowed"
    )
    token = upload_token_from_target(request.target)
    isnothing(token) && return json_upload_response(404, (error="Upload slot not found",))
    registration = registered_upload(gateway.registry, token)
    isnothing(registration) && return json_upload_response(404, (error="Upload slot not found",))

    request.method == "GET" && return lock(registration.lock) do
        perform_upload_download(registration)
    end

    generation = try
        upload_generation(request)
    catch error
        error isa UploadRejected || rethrow()
        return json_upload_response(error.status, (error=error.message,))
    end

    return lock(registration.lock) do
        generation > registration.latest_generation || return json_upload_response(
            409,
            (error="Stale upload operation",)
        )
        registration.latest_generation = generation
        try
            if request.method == "POST"
                Bonito.HTTP.header(request, "X-LCM-Upload", "") == "1" ||
                    throw(UploadRejected(400, "Missing upload request marker"))
                return perform_upload!(registration, request, generation)
            end
            return perform_upload_removal!(registration, generation)
        catch error
            if error isa UploadRejected
                apply_upload_state!(
                    registration.field,
                    generation,
                    4,
                    :failed;
                    progress=0.0,
                    error=error.message
                )
                return json_upload_response(error.status, (error=error.message,))
            end
            diagnostic_id = string(uuid4())
            @error "Upload request failed" diagnostic_id exception=(error, catch_backtrace())
            message = "Upload failed (diagnostic $diagnostic_id)"
            apply_upload_state!(
                registration.field,
                generation,
                4,
                :failed;
                progress=0.0,
                error=message
            )
            return json_upload_response(500, (error=message,))
        end
    end
end

function register_upload_route!(
        server,
        registry::UploadRegistry=default_upload_registry()
    )
    Bonito.route!(server, UPLOAD_ROUTE_PATTERN => UploadGateway(registry))
    return registry
end

function upload_accept_attribute(policy::UploadPolicy)
    values = sort!(collect(union(policy.extensions, policy.media_types)))
    return join(values, ',')
end

function upload_size(bytes::Integer)
    bytes < 1024 && return "$bytes B"
    bytes < 1024^2 && return "$(round(bytes / 1024; digits=1)) KiB"
    return "$(round(bytes / 1024^2; digits=1)) MiB"
end

function upload_status_label(status::Symbol)
    status == :empty && return "Awaiting file"
    status == :selected && return "Selected"
    status == :uploading && return "Uploading"
    status == :validating && return "Validating"
    status == :ready && return "Ready"
    status == :removing && return "Removing"
    return "Upload failed"
end

function upload_icon()
    return SVG.svg(
        SVG.path(; d="M12 16V4"),
        SVG.path(; d="m7 9 5-5 5 5"),
        SVG.path(; d="M5 20h14");
        viewBox="0 0 24 24",
        fill="none",
        stroke="currentColor",
        var"stroke-width"="1.8",
        var"stroke-linecap"="round",
        var"stroke-linejoin"="round",
        var"aria-hidden"="true",
        focusable="false",
        class="lc-upload-icon"
    )
end

function lifecycle_message!(field::UploadField, message::AbstractString)
    isempty(message) && return nothing
    parsed = try
        JSON3.read(message)
    catch
        return nothing
    end
    generation = try Int(parsed.generation) catch; return nothing end
    phase = try Int(parsed.phase) catch; return nothing end
    status = try Symbol(String(parsed.status)) catch; return nothing end
    progress = try Float64(parsed.progress) catch; field.progress[] end
    error = try String(parsed.error) catch; "" end
    apply_upload_state!(field, generation, phase, status; progress, error)
    return nothing
end

function Bonito.jsrender(session::Session, field::UploadField)
    token = register_upload!(field.registry, field.target, field)
    endpoint = "/uploads/$token"
    lifecycle = Observable("")
    on(session, lifecycle) do message
        lifecycle_message!(field, message)
    end
    on(session.on_close) do _
        registration = unregister_upload!(field.registry, token)
        if !isnothing(registration)
            lock(registration.lock) do
                registration.latest_generation += 1
                if field.target.policy.retention == :session
                    delete_upload!(field.target)
                    # The browser-side DOM is already detached while the session
                    # close hooks run. Keep the Julia value accurate without
                    # notifying render listeners that no longer have a target.
                    Observables.setexcludinghandlers!(field.value, nothing)
                end
            end
        end
        return nothing
    end

    current_file = field.value[]
    filename = isnothing(current_file) ? "No file selected" : current_file.original_name
    file_details = isnothing(current_file) ? "" :
        "$(upload_size(current_file.bytes)) · SHA-256 $(first(current_file.sha256, 12))…"
    choose_label = isnothing(current_file) ? "Choose file" : "Replace file"
    initial_busy = field.status[] in (:uploading, :validating, :removing)
    initial_remove_disabled = isnothing(field.value[]) || initial_busy

    input = DOM.input(
        type="file",
        accept=upload_accept_attribute(field.target.policy),
        class="lc-upload-native-input",
        tabindex="-1",
        var"aria-hidden"="true"
    )
    choose = DOM.button(
        upload_icon(),
        DOM.span(choose_label);
        type="button",
        class="lc-upload-choose",
        disabled=initial_busy
    )
    remove = DOM.button(
        "Remove";
        type="button",
        class="lc-upload-remove",
        disabled=initial_remove_disabled
    )
    filename_node = DOM.span(filename; class="lc-upload-filename")
    details_node = DOM.span(file_details; class="lc-upload-details")
    status_node = DOM.span(upload_status_label(field.status[]); class="lc-upload-status")
    progress_node = DOM.progress(
        value=field.progress[],
        max=1,
        class="lc-upload-progress"
    )
    error_node = DOM.div(field.error[]; class="lc-upload-error", role="alert")
    root = DOM.div(
        DOM.style(UPLOAD_STYLES),
        input,
        DOM.div(
            DOM.div(
                DOM.span(field.label; class="lc-upload-label"),
                DOM.span(field.hint; class="lc-upload-hint");
                class="lc-upload-copy"
            ),
            status_node;
            class="lc-upload-header"
        ),
        DOM.div(
            DOM.div(filename_node, details_node; class="lc-upload-file"),
            DOM.div(choose, remove; class="lc-upload-actions");
            class="lc-upload-main"
        ),
        progress_node,
        error_node;
        class="lc-upload-field",
        var"data-upload-state"=string(field.status[]),
        var"data-upload-url"=endpoint,
        var"data-upload-name"=isnothing(field.value[]) ? "" : field.value[].original_name,
        role="group",
        var"aria-label"=field.label
    )

    Bonito.onload(session, root, js"""
        (element) => {
            const input = $(input);
            const choose = $(choose);
            const remove = $(remove);
            const filename = $(filename_node);
            const details = $(details_node);
            const status = $(status_node);
            const progress = $(progress_node);
            const error = $(error_node);
            const lifecycle = $(lifecycle);
            const endpoint = $endpoint;
            const maxBytes = $(field.target.policy.max_bytes);
            const allowEmpty = $(field.target.policy.allow_empty);
            let generation = $(field.generation);
            let request = null;
            let committedName = filename.textContent;
            let committedDetails = details.textContent;
            let lastProgressNotice = 0;

            const publish = (phase, state, ratio, message = "") => {
                lifecycle.notify(JSON.stringify({
                    generation,
                    phase,
                    status: state,
                    progress: ratio,
                    error: message,
                }));
            };
            const show = (state, label, ratio, message = "") => {
                element.dataset.uploadState = state;
                status.textContent = label;
                progress.value = ratio;
                error.textContent = message;
                const busy = state === "uploading" || state === "validating" || state === "removing";
                choose.disabled = busy;
                remove.disabled = busy || committedName === "No file selected";
            };
            const responsePayload = (xhr) => {
                try { return JSON.parse(xhr.responseText || "{}"); }
                catch (_) { return {}; }
            };

            choose.addEventListener("click", () => input.click());
            input.addEventListener("change", () => {
                const file = input.files && input.files[0];
                if (!file) return;
                generation += 1;
                const operation = generation;
                if (request) request.abort();
                filename.textContent = file.name;
                details.textContent = `${file.size.toLocaleString()} B selected`;
                if (file.size > maxBytes || (!allowEmpty && file.size === 0)) {
                    const message = file.size > maxBytes
                        ? `File exceeds the ${maxBytes.toLocaleString()}-byte limit`
                        : "Empty files are not accepted";
                    filename.textContent = committedName;
                    details.textContent = committedDetails;
                    input.value = "";
                    show("failed", "Upload failed", 0, message);
                    publish(4, "failed", 0, message);
                    return;
                }
                show("uploading", "Uploading", 0);
                publish(2, "uploading", 0);
                lastProgressNotice = 0;

                const xhr = new XMLHttpRequest();
                request = xhr;
                xhr.open("POST", endpoint, true);
                xhr.setRequestHeader("Content-Type", file.type || "application/octet-stream");
                xhr.setRequestHeader("X-LCM-Original-Name", encodeURIComponent(file.name));
                xhr.setRequestHeader("X-LCM-Upload-Generation", String(operation));
                xhr.setRequestHeader("X-LCM-Upload", "1");
                xhr.upload.addEventListener("progress", (event) => {
                    if (operation !== generation || !event.lengthComputable) return;
                    const ratio = event.total === 0 ? 0 : event.loaded / event.total;
                    progress.value = ratio;
                    const now = performance.now();
                    if (ratio === 1 || now - lastProgressNotice >= 120) {
                        lastProgressNotice = now;
                        publish(2, "uploading", ratio);
                    }
                });
                xhr.addEventListener("load", () => {
                    if (operation !== generation) return;
                    request = null;
                    input.value = "";
                    const payload = responsePayload(xhr);
                    if (xhr.status >= 200 && xhr.status < 300 && payload.file) {
                        committedName = payload.file.original_name;
                        const shortHash = String(payload.file.sha256 || "").slice(0, 12);
                        committedDetails = `${Number(payload.file.bytes).toLocaleString()} B · SHA-256 ${shortHash}…`;
                        filename.textContent = committedName;
                        details.textContent = committedDetails;
                        choose.querySelector("span").textContent = "Replace file";
                        element.dataset.uploadName = committedName;
                        show("ready", "Ready", 1);
                        publish(4, "ready", 1);
                        element.dispatchEvent(new CustomEvent("lcm:upload-committed", {
                            bubbles: true,
                            detail: {url: endpoint, file: payload.file},
                        }));
                    } else {
                        filename.textContent = committedName;
                        details.textContent = committedDetails;
                        const message = payload.error || `Upload failed (${xhr.status})`;
                        show("failed", "Upload failed", 0, message);
                        publish(4, "failed", 0, message);
                    }
                });
                xhr.addEventListener("error", () => {
                    if (operation !== generation) return;
                    request = null;
                    input.value = "";
                    filename.textContent = committedName;
                    details.textContent = committedDetails;
                    const message = "The upload connection failed";
                    show("failed", "Upload failed", 0, message);
                    publish(4, "failed", 0, message);
                });
                xhr.send(file);
            });

            remove.addEventListener("click", () => {
                generation += 1;
                const operation = generation;
                if (request) request.abort();
                show("removing", "Removing", 0);
                publish(2, "removing", 0);
                const xhr = new XMLHttpRequest();
                request = xhr;
                xhr.open("DELETE", endpoint, true);
                xhr.setRequestHeader("X-LCM-Upload-Generation", String(operation));
                xhr.addEventListener("load", () => {
                    if (operation !== generation) return;
                    request = null;
                    const payload = responsePayload(xhr);
                    if (xhr.status >= 200 && xhr.status < 300) {
                        committedName = "No file selected";
                        committedDetails = "";
                        filename.textContent = committedName;
                        details.textContent = committedDetails;
                        choose.querySelector("span").textContent = "Choose file";
                        element.dataset.uploadName = "";
                        show("empty", "Awaiting file", 0);
                        publish(4, "empty", 0);
                        element.dispatchEvent(new CustomEvent("lcm:upload-removed", {
                            bubbles: true,
                        }));
                    } else {
                        const message = payload.error || `Removal failed (${xhr.status})`;
                        show("failed", "Removal failed", 0, message);
                        publish(4, "failed", 0, message);
                    }
                });
                xhr.addEventListener("error", () => {
                    if (operation !== generation) return;
                    request = null;
                    const message = "The removal connection failed";
                    show("failed", "Removal failed", 0, message);
                    publish(4, "failed", 0, message);
                });
                xhr.send();
            });

            window.addEventListener("pagehide", () => {
                if (request) request.abort();
            }, {once: true});
        }
    """)
    return Bonito.jsrender(session, root)
end
