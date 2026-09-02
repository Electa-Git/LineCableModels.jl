@testset "architectural invariants" begin
    root = normpath(joinpath(@__DIR__, ".."))

    function files_below(directory, extensions)
        isdir(directory) || return String[]
        return sort!(filter(
            path -> any(extension -> endswith(path, extension), extensions),
            collect(Iterators.flatten(
                joinpath.(walk_root, files)
                for (walk_root, _, files) in walkdir(directory)
            ))
        ))
    end

    function joined_sources(directory, extensions=(".jl",))
        return join((read(path, String) for path in files_below(directory, extensions)), '\n')
    end

    publisher_project = TOML.parsefile(joinpath(root, "Project.toml"))
    @test !haskey(publisher_project, "apps")
    forbidden = Set((
        "LineCableModels",
        "PowerImpedance",
        "PowerModels",
        "PowerModelsACDC",
        "InfrastructureModels",
        "JuMP",
        "Ipopt",
        "NonlinearSolve",
    ))
    @test isempty(intersect(forbidden, Set(keys(publisher_project["deps"]))))

    launcher = read(joinpath(root, "lcm"), String)
    @test occursin("-m LineCableModelsPlayground", launcher)
    @test occursin("-m LineCableModelsWorker", launcher)
    @test occursin("--project=\"\${SCRIPT_DIR}/worker\"", launcher)

    protocol_project = TOML.parsefile(joinpath(root, "protocol", "Project.toml"))
    @test Set(keys(protocol_project["deps"])) == Set((
        "Dates", "JSON3", "SHA", "StructTypes", "UUIDs"
    ))
    protocol_source = joined_sources(joinpath(root, "protocol", "src"))
    @test !occursin("Serialization", protocol_source)
    @test !occursin(r"\beval\s*\(", protocol_source)
    @test !occursin(
        r"(?m)^\s*(?:using|import)\s+(?:LineCableModels|PowerImpedance|NATS|Bonito)\b",
        protocol_source
    )
    @test !occursin("Function", read(joinpath(root, "protocol", "src", "Jobs.jl"), String))

    worker_source = joined_sources(joinpath(root, "worker", "src"))
    @test !occursin("JetChannel", worker_source)
    @test occursin("consumer_ack", worker_source)
    @test occursin("stream_publish", worker_source)
    @test occursin("execute_supervised!", worker_source)
    @test !occursin(r"register!\([^\n]+\"(?:eval|repl)\"", worker_source)

    browser_files = unique!(vcat(
        files_below(root, (".qmd",)),
        files_below(joinpath(root, "_site"), (".html", ".js", ".css", ".json")),
        files_below(joinpath(root, "assets"), (".html", ".js", ".css", ".svg")),
        files_below(joinpath(root, "_extensions"), (".lua", ".js", ".css", ".html")),
    ))
    browser_source = join((read(path, String) for path in browser_files), '\n')
    @test !occursin(r"(?i)nats(?:\+tls)?://", browser_source)
    @test !occursin(r"NATS_(?:CONNECT|ADMIN)_URL", browser_source)
    @test !occursin(r"NATS_[A-Z_]*PASSWORD", browser_source)
    @test !occursin("AWS_SECRET_ACCESS_KEY", browser_source)
    @test !occursin(r"\bsubmit!\s*\(", browser_source)

    client = BrokerClient(autostart=false)
    panel = JobPanel(
        client;
        operation="system.echo",
        parameters=Dict("message" => "not submitted"),
        session_id="architecture-test"
    )
    @test isempty(client.handles)
    @test panel.handle.state[] == :editable
    @test client.connection_status[] == :offline
    for (_, factory) in LineCableModelsPlayground.WIDGET_ROUTES
        @test !isnothing(factory())
    end
    @test !isnothing(LineCableModelsPlayground.runtime_job_widget(client))
    @test isempty(client.handles)

    workbench_source = joined_sources(joinpath(root, "src", "workbench"))
    template_workbench_source = joined_sources(
        joinpath(root, "src", "workbenches")
    )
    @test !occursin(
        r"(?m)^\s*(?:using|import)\s+(?:NATS|LineCableModels|PowerImpedance)\b",
        workbench_source
    )
    @test !occursin(
        r"(?m)^\s*(?:using|import)\s+(?:NATS|LineCableModels|PowerImpedance)\b",
        template_workbench_source
    )
    @test occursin("abstract type AbstractWorkbenchAction", workbench_source)
    @test occursin("handle!(runtime.application", workbench_source)
    @test occursin("ViewStack", workbench_source)
    @test occursin("/workbenches/template", read(
        joinpath(root, "src", "LineCableModelsPlayground.jl"),
        String
    ))

    job_controls = read(
        joinpath(root, "src", "widgets", "JobControls.jl"),
        String
    )
    @test occursin("submit_async!", job_controls)
    @test !occursin(r"\bwait\s*\(", job_controls)

    patch_notes = read(joinpath(root, "vendor", "NATS", "LCM_PATCHES.md"), String)
    @test occursin("No certificate verification bypass", patch_notes)

    local_compose = read(joinpath(root, "deploy", "compose.yml"), String)
    remote_compose = read(joinpath(root, "deploy", "remote", "compose.yml"), String)
    local_cpu_limits = read(
        joinpath(root, "deploy", "compose.cpu-limits.yml"),
        String
    )
    remote_cpu_limits = read(
        joinpath(root, "deploy", "remote", "compose.cpu-limits.yml"),
        String
    )
    @test count("read_only: true", local_compose) == 4
    @test count("no-new-privileges:true", local_compose) == 4
    @test count("pids_limit:", local_compose) == 4
    @test count("mem_limit:", local_compose) == 4
    @test !occursin("cpus:", local_compose)
    @test count("cpus:", local_cpu_limits) == 4
    @test count("read_only: true", remote_compose) == 6
    @test count("no-new-privileges:true", remote_compose) == 6
    @test count("pids_limit:", remote_compose) == 6
    @test count("mem_limit:", remote_compose) == 6
    @test !occursin("cpus:", remote_compose)
    @test count("cpus:", remote_cpu_limits) == 6
    @test !occursin(r"(?m)^\s+-\s+\.\./.*:/workspace", local_compose)
    @test !occursin(r"(?m)^\s+-\s+\.\./.*:/workspace", remote_compose)
    @test occursin("NATS_TLS_CERT_PATH", remote_compose)
    @test occursin("MINIO_WORKER_ACCESS_KEY", remote_compose)
    @test occursin("MINIO_PUBLISHER_ACCESS_KEY", remote_compose)
    @test occursin(r"(?s)worker:.*?runtime-init:.*?service_completed_successfully", local_compose)
    @test occursin(r"(?s)worker:.*?runtime-init:.*?service_completed_successfully", remote_compose)

    local_ignore = read(joinpath(root, ".gitignore"), String)
    docker_ignore = read(joinpath(root, "..", ".dockerignore"), String)
    @test occursin("/deploy/.env", local_ignore)
    @test occursin("/deploy/remote/.env", local_ignore)
    @test occursin("/deploy/remote/certs/", local_ignore)
    @test occursin("playground/deploy/.env", docker_ignore)
    @test occursin("playground/deploy/remote/.env", docker_ignore)
    @test occursin("playground/deploy/remote/certs", docker_ignore)

    image_sources = join((
        local_compose,
        remote_compose,
        read(joinpath(root, "deploy", "publisher.Dockerfile"), String),
        read(joinpath(root, "deploy", "worker.Dockerfile"), String),
    ), '\n')
    for line in split(image_sources, '\n')
        stripped = strip(line)
        if startswith(stripped, "image:") || startswith(stripped, "FROM ")
            @test occursin("@sha256:", stripped)
        end
    end
end
