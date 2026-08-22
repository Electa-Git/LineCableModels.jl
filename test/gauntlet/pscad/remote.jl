struct RemoteConfig
    host::String
    shared_root::String
    remote_root::String
    julia_executable::String
    python_executable::String
    pscad_version::String
    transport::Symbol
    verbosity::Int
    timeout_seconds::Int
end

function RemoteConfig(
        host::AbstractString,
        shared_root::AbstractString,
        remote_root::AbstractString,
        julia_executable::AbstractString,
        python_executable::AbstractString;
        pscad_version::AbstractString = "5.1.0",
        transport::Symbol = :ssh,
        verbosity::Integer = 1,
        timeout_seconds::Integer = 1800
)
    isempty(strip(host)) && throw(ArgumentError("PSCAD host cannot be empty"))
    isempty(strip(shared_root)) && throw(ArgumentError(
        "PSCAD shared root cannot be empty",
    ))
    isempty(strip(remote_root)) && throw(ArgumentError("PSCAD remote root cannot be empty"))
    isempty(strip(julia_executable)) && throw(ArgumentError(
        "PSCAD-host Julia executable cannot be empty",
    ))
    isempty(strip(python_executable)) && throw(ArgumentError(
        "PSCAD-host Python executable cannot be empty",
    ))
    pscad_version == "5.1.0" || throw(ArgumentError(
        "the gauntlet supports PSCAD 5.1.0 only",
    ))
    verbosity in 0:2 || throw(ArgumentError(
        "PSCAD verbosity must be 0, 1, or 2",
    ))
    timeout_seconds > 0 || throw(ArgumentError(
        "PSCAD timeout_seconds must be positive",
    ))
    return RemoteConfig(
        String(host),
        String(shared_root),
        String(remote_root),
        String(julia_executable),
        String(python_executable),
        String(pscad_version),
        transport,
        Int(verbosity),
        Int(timeout_seconds)
    )
end

function _powershell_argv(powershell::AbstractString)
    utf16 = htol.(transcode(UInt16, String(powershell)))
    encoded = base64encode(reinterpret(UInt8, utf16))
    return [
        "powershell.exe", "-NoProfile", "-NonInteractive",
        "-EncodedCommand", encoded
    ]
end

function remote_command(::Val{:ssh}, config::RemoteConfig, powershell::AbstractString)
    return Cmd(vcat(["ssh", config.host], _powershell_argv(powershell)))
end

function remote_command(config::RemoteConfig, powershell::AbstractString)
    # `local.jl` is loaded only for an explicit live/record run and may add the
    # one user-owned transport method after this file has been compiled.
    return Base.invokelatest(
        remote_command,
        Val(config.transport),
        config,
        powershell
    )
end

function _ps_quote(value::AbstractString)
    return "'" * replace(String(value), "'" => "''") * "'"
end

function _remote_path(root::AbstractString, parts::AbstractString...)
    clean = rstrip(String(root), ('\\', '/'))
    return join((clean, parts...), '\\')
end

function _remote_project_name(path::AbstractString)
    return splitext(basename(replace(String(path), '\\' => '/')))[1]
end

function _run_remote(
        config::RemoteConfig,
        powershell::AbstractString;
        stdout_path::Union{Nothing, AbstractString} = nothing,
        stderr_path::Union{Nothing, AbstractString} = nothing,
        stream::Bool = false,
        on_interrupt::Function = () -> nothing
)
    stdout_path === nothing || mkpath(dirname(stdout_path))
    stderr_path === nothing || mkpath(dirname(stderr_path))
    output = Pipe()
    errors = Pipe()
    process = run(
        pipeline(remote_command(config, powershell); stdout = output, stderr = errors);
        wait = false
    )
    close(output.in)
    close(errors.in)
    output_buffer = IOBuffer()
    error_buffer = IOBuffer()
    output_file = stdout_path === nothing ? nothing : open(stdout_path, "w")
    error_file = stderr_path === nothing ? nothing : open(stderr_path, "w")
    output_task = @async for line in eachline(output)
        println(output_buffer, line)
        output_file === nothing || println(output_file, line)
        output_file === nothing || flush(output_file)
        if stream
            println(stdout, line)
            flush(stdout)
        end
    end
    error_task = @async for line in eachline(errors)
        println(error_buffer, line)
        error_file === nothing || println(error_file, line)
        error_file === nothing || flush(error_file)
        if stream
            println(stderr, line)
            flush(stderr)
        end
    end
    interrupted = nothing
    try
        wait(process)
    catch error
        interrupted = error
        process_running(process) && kill(process)
        if error isa InterruptException
            try
                on_interrupt()
            catch cancellation_error
                @warn "Remote PSCAD cancellation could not be confirmed" exception = (
                    cancellation_error, catch_backtrace())
            end
        end
        try
            wait(process)
        catch
        end
    finally
        if interrupted !== nothing
            close(output)
            close(errors)
        end
        try
            wait(output_task)
        catch
        end
        try
            wait(error_task)
        catch
        end
        output_file === nothing || close(output_file)
        error_file === nothing || close(error_file)
    end
    stdout_value = String(take!(output_buffer))
    stderr_value = String(take!(error_buffer))
    interrupted === nothing || throw(interrupted)
    success(process) || throw(ErrorException(
        "remote PSCAD command failed with exit code $(process.exitcode)\n" *
        "stdout:\n$stdout_value\nstderr:\n$stderr_value",
    ))
    return stdout_value
end

function _load_config()
    path = joinpath(GAUNTLET_ROOT, "local.jl")
    isfile(path) || throw(ArgumentError(
        "live PSCAD validation requires $path; copy and edit local.example",
    ))
    config = Base.include(@__MODULE__, path)
    config isa RemoteConfig || throw(ArgumentError(
        "$path must return PSCADHarness.RemoteConfig",
    ))
    return config
end

function _validate_frequencies(frequencies_value::AbstractVector)
    length(frequencies_value) >= 101 || throw(ArgumentError(
        "PSCAD 5.1 requires at least 100 frequency increments (101 samples)",
    ))
    all(>(0), frequencies_value) || throw(DomainError(
        frequencies_value,
        "PSCAD frequencies must be positive"
    ))
    expected = 10.0 .^ range(
        log10(first(frequencies_value));
        stop = log10(last(frequencies_value)),
        length = length(frequencies_value)
    )
    collect(frequencies_value) == expected || throw(ArgumentError(
        "PSCAD live validation supports an exact logarithmic frequency vector only",
    ))
    return frequencies_value
end

function _supervisor_command(
        config::RemoteConfig,
        shared_case::AbstractString,
        remote_case::AbstractString,
        project_name::AbstractString,
        formulation::Symbol,
        frequencies_value::AbstractVector
)
    _validate_frequencies(frequencies_value)
    increments = length(frequencies_value) - 1
    shared_supervisor = _remote_path(shared_case, "toolkit", "supervisor.ps1")
    return join(
        (
            "& $(_ps_quote(shared_supervisor))",
            "-SharedCase $(_ps_quote(shared_case))",
            "-LocalCase $(_ps_quote(remote_case))",
            "-Julia $(_ps_quote(config.julia_executable))",
            "-Python $(_ps_quote(config.python_executable))",
            "-ProjectName $(_ps_quote(project_name))",
            "-Formulation $(_ps_quote(string(formulation)))",
            "-FrequencyStart $(_ps_quote(string(first(frequencies_value))))",
            "-FrequencyEnd $(_ps_quote(string(last(frequencies_value))))",
            "-FrequencyIncrements $(_ps_quote(string(increments)))",
            "-PSCADVersion $(_ps_quote(config.pscad_version))",
            "-Verbosity $(_ps_quote(string(config.verbosity)))",
            "-TimeoutSeconds $(_ps_quote(string(config.timeout_seconds)))"
        ),
        ' ')
end

function _cancel_command(remote_case::AbstractString)
    owner_path = _remote_path(remote_case, "owner.txt")
    return "\$ownerPath=$(_ps_quote(owner_path)); " *
           "if (-not (Test-Path -LiteralPath \$ownerPath -PathType Leaf)) { exit 0 }; " *
           "\$owner=@(Get-Content -LiteralPath \$ownerPath); " *
           "if (\$owner.Count -ne 2) { throw 'invalid PSCAD runner owner file' }; " *
           "\$runner=[IO.Path]::GetFullPath(\$owner[1]); " *
           "\$root=[IO.Path]::GetFullPath($(_ps_quote(remote_case))).TrimEnd('\\')+'\\'; " *
           "if (-not \$runner.StartsWith(\$root,[StringComparison]::OrdinalIgnoreCase)) " *
           "{ throw 'refusing to stop a process outside the PSCAD case directory' }; " *
           "\$runnerPid=0; if (-not [int]::TryParse(\$owner[0],[ref]\$runnerPid)) " *
           "{ throw 'invalid PSCAD runner PID' }; " *
           "\$process=Get-CimInstance Win32_Process -Filter \"ProcessId = \$runnerPid\"; " *
           "if (\$null -eq \$process) { Remove-Item -LiteralPath \$ownerPath -Force; exit 0 }; " *
           "if (\$null -eq \$process.CommandLine -or " *
           "\$process.CommandLine.IndexOf(\$runner,[StringComparison]::OrdinalIgnoreCase) -lt 0) " *
           "{ throw 'recorded PID no longer belongs to the PSCAD runner' }; " *
           "& taskkill.exe /PID \$runnerPid /T /F | Out-Null; " *
           "if (\$LASTEXITCODE -ne 0) { throw 'could not stop PSCAD runner process tree' }; " *
           "Remove-Item -LiteralPath \$ownerPath -Force -ErrorAction SilentlyContinue"
end

function _cancel_remote(config::RemoteConfig, remote_case::AbstractString)
    _run_remote(config, _cancel_command(remote_case); stream = config.verbosity >= 2)
    return nothing
end

function _stage_toolkit(local_project::AbstractString, local_output::AbstractString)
    isfile(local_project) || throw(ArgumentError(
        "local PSCAD input is missing: $local_project",
    ))
    variant_root = dirname(local_output)
    expected_project = joinpath(variant_root, "generated.pscx")
    abspath(local_project) == abspath(expected_project) || throw(ArgumentError(
        "PSCAD project must be staged as $expected_project",
    ))
    toolkit_source = joinpath(@__DIR__, "remote")
    toolkit_stage = joinpath(variant_root, "toolkit")
    isdir(toolkit_stage) && rm(toolkit_stage; recursive = true)
    mkpath(toolkit_stage)
    for name in ("Project.toml", "Manifest.toml", "runner.jl", "supervisor.ps1")
        source = joinpath(toolkit_source, name)
        isfile(source) || throw(ArgumentError("PSCAD toolkit file is missing: $source"))
        cp(source, joinpath(toolkit_stage, name); force = true)
    end
    return toolkit_stage
end

function run_remote_pscad(
        config::RemoteConfig,
        local_project::AbstractString,
        local_output::AbstractString,
        formulation::Symbol,
        frequencies_value::AbstractVector
)
    formulation_spec(formulation)
    _validate_frequencies(frequencies_value)
    isdir(local_output) && rm(local_output; recursive = true)
    mkpath(local_output)
    variant = basename(dirname(local_output))
    case_name = basename(dirname(dirname(local_output)))
    _stage_toolkit(local_project, local_output)
    shared_case = _remote_path(config.shared_root, case_name, variant)
    remote_case = _remote_path(
        config.remote_root,
        case_name,
        variant
    )
    stdout_path = joinpath(local_output, "stdout.txt")
    stderr_path = joinpath(local_output, "stderr.txt")
    transport_stdout = joinpath(local_output, "transport-stdout.txt")
    transport_stderr = joinpath(local_output, "transport-stderr.txt")
    command = _supervisor_command(
        config,
        shared_case,
        remote_case,
        _remote_project_name(local_project),
        formulation,
        frequencies_value
    )
    @info "Executing PSCAD frequency scan" host=config.host variant formulation frequencies=length(frequencies_value) timeout_seconds=config.timeout_seconds
    execution_error = try
        _run_remote(
            config,
            command;
            stdout_path = transport_stdout,
            stderr_path = transport_stderr,
            stream = config.verbosity >= 2,
            on_interrupt = () -> _cancel_remote(config, remote_case)
        )
        nothing
    catch error
        if !(error isa InterruptException)
            try
                _cancel_remote(config, remote_case)
            catch cancellation_error
                @warn "Remote PSCAD cancellation could not be confirmed" host=config.host variant exception=(
                    cancellation_error, catch_backtrace())
            end
        end
        error isa InterruptException && rethrow()
        error
    end
    if execution_error !== nothing
        console_path = joinpath(local_output, "pscad-console.txt")
        if !isfile(console_path)
            write(
                console_path,
                "PSCAD did not produce a console log before the remote failure.\n"
            )
        end
        throw(ErrorException(
            sprint(showerror, execution_error) *
            "\nPSCAD diagnostics: $console_path" *
            "\nRemote scratch: $remote_case",
        ))
    end
    @info "Checking PSCAD outputs" host=config.host variant destination=local_output
    required = (
        "pscad-console.txt", "timing.txt", "result_zm.out", "result_zp.out",
        "result_ym.out", "result_yp.out"
    )
    for name in required
        path = joinpath(local_output, name)
        isfile(path) || throw(ArgumentError("required PSCAD output is missing: $path"))
        filesize(path) > 0 || throw(ArgumentError("required PSCAD output is empty: $path"))
    end
    elapsed = parse(Float64, strip(read(joinpath(local_output, "timing.txt"), String)))
    @info "PSCAD frequency scan completed" host=config.host variant elapsed_seconds=elapsed
    return (
        elapsed_seconds = elapsed,
        exit_code = 0,
        stdout_path,
        stderr_path,
        console_path = joinpath(local_output, "pscad-console.txt"),
        output_dir = String(local_output)
    )
end

function _assert_error(error, tolerance)
    error.absolute <= tolerance.absolute || error.relative <= tolerance.relative
end

function _assert_recordable(case::GauntletCase, current, legacy)
    tolerances = case.tolerances
    _assert_error(current.Z, tolerances.local_vs_pscad.Z) || throw(ArgumentError(
        "Example PSCAD/local Z comparison exceeds the case tolerance",
    ))
    _assert_error(current.Y, tolerances.local_vs_pscad.Y) || throw(ArgumentError(
        "Example PSCAD/local Y comparison exceeds the case tolerance",
    ))
    if case.compare_legacy_exporter
        _assert_error(legacy.Z, tolerances.legacy_vs_current.Z) || throw(ArgumentError(
            "legacy/current PSCAD Z comparison exceeds the case tolerance",
        ))
        _assert_error(legacy.Y, tolerances.legacy_vs_current.Y) || throw(ArgumentError(
            "legacy/current PSCAD Y comparison exceeds the case tolerance",
        ))
    end
    return nothing
end

function run_live(case::GauntletCase; record::Bool = false)
    validate_formulation(case)
    config = _load_config()
    root = work_path(case)
    @info "Starting PSCAD gauntlet case" case=case.name mode=record ? :record : :live work_directory=root
    isdir(root) && rm(root; recursive = true)
    current_root = joinpath(root, "current")
    legacy_root = joinpath(root, "legacy")
    mkpath(current_root)
    mkpath(legacy_root)
    system = case.problem.system
    earth = case.problem.earth_props
    @info "Exporting current PSCAD project" case = case.name
    current_project = export_data(
        :pscad,
        system,
        earth;
        file_name = joinpath(current_root, "generated.pscx")
    )
    current_project isa AbstractString && isfile(current_project) || throw(ArgumentError(
        "the current PSCAD exporter did not create the case project",
    ))
    current_copy = joinpath(current_root, "generated.pscx")
    current_project == current_copy || cp(current_project, current_copy; force = true)
    current_project = current_copy
    @info "Exporting legacy PSCAD project" case = case.name
    legacy_project = export_legacy_pscad(
        system,
        earth;
        file_name = joinpath(legacy_root, "generated.pscx")
    )
    legacy_copy = joinpath(legacy_root, "generated.pscx")
    legacy_project == legacy_copy || cp(legacy_project, legacy_copy; force = true)
    legacy_project = legacy_copy
    @info "Running current PSCAD export" case = case.name
    current_run = run_remote_pscad(
        config,
        current_project,
        joinpath(current_root, "outputs"),
        case.pscad_formulation,
        case.problem.frequencies
    )
    @info "Running legacy PSCAD export" case = case.name
    legacy_run = run_remote_pscad(
        config,
        legacy_project,
        joinpath(legacy_root, "outputs"),
        case.pscad_formulation,
        case.problem.frequencies
    )
    @info "Parsing PSCAD line-parameter outputs" case = case.name
    reference = read_pscad_result(
        current_run.output_dir,
        case.problem.frequencies,
        case.expected_size
    )
    legacy_reference = read_pscad_result(
        legacy_run.output_dir,
        case.problem.frequencies,
        case.expected_size
    )
    @info "Computing LineCableModels reference" case = case.name
    candidate = compute!(case.problem, case.formulation)
    validate_structure(case, reference)
    validate_structure(case, legacy_reference)
    validate_structure(case, candidate)
    current_comparison = compare(reference, candidate)
    legacy_comparison = compare(reference, legacy_reference)
    @info "Benchmarking LineCableModels calculation" case = case.name
    local_timing = benchmark_local(case)
    record && _assert_recordable(case, current_comparison, legacy_comparison)
    if record
        @info "Recording PSCAD snapshot" case=case.name destination=snapshot_path(case)
        write_snapshot(
            snapshot_path(case),
            case,
            reference;
            pscad_version = config.pscad_version,
            pscad_elapsed_seconds = current_run.elapsed_seconds,
            julia_benchmark = local_timing
        )
    end
    @info "PSCAD gauntlet case completed" case=case.name mode=record ? :record : :live
    return (
        mode = record ? :record : :live,
        reference,
        candidate,
        comparison = current_comparison,
        metadata = nothing,
        pscad = current_run,
        legacy = (
            reference = legacy_reference,
            comparison = legacy_comparison,
            run = legacy_run
        ),
        timings = (pscad = current_run.elapsed_seconds, julia = local_timing)
    )
end
