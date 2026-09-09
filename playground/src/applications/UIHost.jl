"""
    serve_owned_ui(routes; startup_seconds=90)

Serve approved route/App pairs in the current, isolated UI process. Validate
the supervisor-provided run namespace, admit only its private proxy credential,
and render each route once before publishing an atomic readiness receipt.
Close probe sessions immediately. No worker preparation or numerical work is
permitted in these route constructors.

This entry point is for a supervised child script, not the publisher process.
It blocks until shutdown and closes its Bonito server on exit.
"""
function serve_owned_ui(routes; startup_seconds::Real=90)
    startup_seconds > 0 && isfinite(startup_seconds) ||
        throw(ArgumentError("UI preparation deadline must be finite and positive"))
    run_id = UUID(ENV["LCM_RUN_ID"])
    parent_pid = parse(Int, ENV["LCM_UI_PARENT_PID"])
    Sys.islinux() && ccall(:getppid, Cint, ()) == parent_pid ||
        throw(ArgumentError("UI host must remain attached to its supervisor"))
    prefix = ENV["LCM_RUN_PREFIX"]
    prefix == "/applications/runs/$run_id/" || throw(ArgumentError("invalid owned UI namespace"))
    key = ENV["LCM_UI_HOST_KEY"]
    occursin(r"^[a-f0-9]{64}$", key) || throw(ArgumentError("invalid private UI credential"))
    file = abspath(ENV["LCM_UI_READY_FILE"])
    basename(dirname(file)) == string(run_id) || throw(ArgumentError("invalid owned UI receipt directory"))
    basename(file) == "ready.json" && !ispath(file) || throw(ArgumentError("UI receipt already exists"))
    entries = collect(routes)
    !isempty(entries) && allunique(first.(entries)) || throw(ArgumentError("UI routes must be nonempty and unique"))
    for (path, app) in entries
        path isa AbstractString && startswith(path, "/") && !startswith(path, "//") &&
            !occursin(r"[\\?#\s]|/\.\.?(/|$)", path) && path != "/health" ||
            throw(ArgumentError("invalid owned UI route"))
        app isa Bonito.App || throw(ArgumentError("owned UI routes must contain Bonito apps"))
    end
    digest = sha256(key)
    function admit_private_connection(request, peer)
        keys = [last(h) for h in request.headers if lowercase(first(h)) == "x-lcm-host-key"]
        length(keys) == 1 || error("Private UI connection denied")
        offered = sha256(only(keys))
        difference = zero(UInt8)
        for i in eachindex(digest)
            difference |= digest[i] ⊻ offered[i]
        end
        iszero(difference) || error("Private UI connection denied")
    end
    server = Bonito.Server("127.0.0.1", 0; proxy_url=prefix, access_log=admit_private_connection)
    try
        Bonito.route!(server, "/health" => (_ -> Bonito.HTTP.Response(200, string(run_id))))
        register_upload_route!(server)
        for (route, app) in entries
            Bonito.route!(server, route => app)
        end
        for (route, app) in entries
            try
                response = Bonito.HTTP.get("http://127.0.0.1:$(server.port)$route";
                    headers=["X-LCM-Host-Key"=>key], proxy=nothing,
                    request_timeout=startup_seconds, retry=false, redirect=false)
                response.status == 200 || error("UI preparation request failed")
                app.session[].init_error[] === nothing || error("UI preparation render failed")
            finally
                app.session[] === nothing || close(app.session[])
            end
        end
        ccall(:getppid, Cint, ()) == parent_pid || error("UI supervisor was lost during startup")
        write(file * ".pending", JSON3.write((schema_version=1, run_id=string(run_id),
            pid=getpid(), port=server.port)))
        mv(file * ".pending", file)
        wait(server)
    finally
        close(server)
    end
end
