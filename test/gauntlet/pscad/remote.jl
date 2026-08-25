struct RemoteConfig
    host::String
    shared_root::String
    remote_root::String
    julia_executable::String
    python_executable::String
    pscad_version::String
    transport::Symbol
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
        verbosity = nothing,
        timeout_seconds::Integer = 1800
)
    verbosity === nothing || throw(ArgumentError(
        "RemoteConfig no longer owns verbosity; set it with " *
        "options=(reference=(verbosity=(default=0, PSCAD=2),),) in the case file",
    ))
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
        Int(timeout_seconds)
    )
end

function computation_options(
        ::Val{PSCADFormulation},
        options::NamedTuple
)::ComputationOptions
    allowed = (:output_stem, :remote, :verbosity, :output_basis)
    unknown = filter(key -> key ∉ allowed, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown PSCAD computation options: $(sort!(collect(unknown)))",
    ))
    haskey(options, :remote) || throw(ArgumentError(
        "PSCAD computation options must define remote::RemoteConfig",
    ))
    options.remote isa RemoteConfig || throw(ArgumentError(
        "PSCAD computation option remote must be a RemoteConfig",
    ))
    normalized = merge(
        (
            output_stem = "gauntlet",
            verbosity = (default = 0,),
            output_basis = :pul
        ),
        options
    )
    normalized.output_stem isa AbstractString || throw(ArgumentError(
        "PSCAD output_stem must be a string",
    ))
    output_stem = String(normalized.output_stem)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_]{0,19}$", output_stem) ||
        throw(ArgumentError(
            "PSCAD output_stem must contain 1–20 ASCII letters, digits, or underscores",
        ))
    verbosity_values = normalized.verbosity
    verbosity_values isa NamedTuple || throw(ArgumentError(
        "verbosity must be a named tuple",
    ))
    haskey(verbosity_values, :default) || throw(ArgumentError(
        "verbosity must define a default level",
    ))
    all(value -> value isa Integer && value in 0:2, values(verbosity_values)) ||
        throw(ArgumentError("verbosity levels must be integers from 0 to 2"))
    basis_value = normalized.output_basis
    basis_value in (:pul, :total) || throw(ArgumentError(
        "output_basis must be :pul or :total; got $(repr(basis_value))",
    ))
    levels = NamedTuple{keys(verbosity_values)}(Int.(values(verbosity_values)))
    return (
        output_stem,
        remote = options.remote,
        verbosity = levels,
        output_basis = Val(basis_value)
    )
end

function _powershell_argv(powershell::AbstractString)
    command="\$ProgressPreference='SilentlyContinue'; $(String(powershell))"
    utf16 = htol.(transcode(UInt16, command))
    encoded = base64encode(reinterpret(UInt8, utf16))
    return [
        "powershell.exe", "-NoProfile", "-NonInteractive",
        "-EncodedCommand", encoded
    ]
end

function _formulation_label(formulation::PSCADFormulation)
    return join(
        (
            description(formulation.earth_impedance),
            description(formulation.earth_admittance),
            description(formulation.insulation_admittance)
        ),
        '/')
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
    path = get(
        ENV,
        "LINECABLEMODELS_GAUNTLET_CONFIG",
        joinpath(GAUNTLET_ROOT, "local.jl")
    )
    path = abspath(path)
    isfile(path) || throw(ArgumentError(
        "live PSCAD validation requires $path; copy and edit local.example",
    ))
    config = Base.include(@__MODULE__, path)
    config isa RemoteConfig || throw(ArgumentError(
        "$path must return PSCADBenchmarks.RemoteConfig",
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
        formulation::PSCADFormulation,
        frequencies_value::AbstractVector;
        output_stem::AbstractString,
        verbosity::Integer = 0
)
    _validate_frequencies(frequencies_value)
    earth = formulation.earth_impedance
    label = _formulation_label(formulation)
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
            "-OutputStem $(_ps_quote(output_stem))",
            "-Formulation $(_ps_quote(label))",
            "-EarthField $(_ps_quote(string(pscad_field(earth))))",
            "-EarthValue $(_ps_quote(string(pscad_value(earth))))",
            "-EarthReadback $(_ps_quote(pscad_readback(earth)))",
            "-FrequencyStart $(_ps_quote(string(first(frequencies_value))))",
            "-FrequencyEnd $(_ps_quote(string(last(frequencies_value))))",
            "-FrequencyIncrements $(_ps_quote(string(increments)))",
            "-PSCADVersion $(_ps_quote(config.pscad_version))",
            "-Verbosity $(_ps_quote(string(verbosity)))",
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

function _cancel_remote(
        config::RemoteConfig,
        remote_case::AbstractString;
        verbosity::Integer = 0
)
    _run_remote(config, _cancel_command(remote_case); stream = verbosity >= 2)
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
    for name in (
        "Project.toml", "Manifest.toml", "files.jl", "runner.jl", "supervisor.ps1"
    )
        source = joinpath(toolkit_source, name)
        isfile(source) || throw(ArgumentError("PSCAD toolkit file is missing: $source"))
        cp(source, joinpath(toolkit_stage, name); force = true)
    end
    return toolkit_stage
end

function _diagnostic_tail(path::AbstractString; count::Integer = 12)
    isfile(path) || return "PSCAD produced no diagnostic log."
    lines=filter(!isempty, strip.(readlines(path)))
    isempty(lines) && return "PSCAD diagnostic log is empty."
    return join(last(lines, min(count, length(lines))), '\n')
end

function _work_parts(local_output::AbstractString)
    output=abspath(local_output)
    basename(output) == "outputs" || throw(ArgumentError(
        "PSCAD output directory must end in 'outputs': $local_output",
    ))
    relative=relpath(dirname(output), abspath(WORK_ROOT))
    parts=splitpath(relative)
    length(parts) == 3 && all(part -> part ∉ ("", ".", ".."), parts) ||
        throw(ArgumentError(
            "PSCAD work directory must be WORK_ROOT/<backend>/<case>/<variant>",
        ))
    return parts
end

function run_remote_pscad(
        config::RemoteConfig,
        local_project::AbstractString,
        local_output::AbstractString,
        formulation::PSCADFormulation,
        frequencies_value::AbstractVector;
        output_stem::AbstractString,
        verbosity::Integer = 0
)
    verbosity in 0:2 || throw(ArgumentError("PSCAD verbosity must be 0, 1, or 2"))
    _validate_frequencies(frequencies_value)
    isdir(local_output) && rm(local_output; recursive = true)
    mkpath(local_output)
    work_parts = _work_parts(local_output)
    variant = last(work_parts)
    _stage_toolkit(local_project, local_output)
    shared_case = _remote_path(config.shared_root, work_parts...)
    remote_case = _remote_path(config.remote_root, work_parts...)
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
        frequencies_value;
        output_stem,
        verbosity
    )
    @info "Executing PSCAD frequency scan" host=config.host variant formulation=_formulation_label(formulation) frequencies=length(frequencies_value) timeout_seconds=config.timeout_seconds
    execution_error = try
        _run_remote(
            config,
            command;
            stdout_path = transport_stdout,
            stderr_path = transport_stderr,
            stream = verbosity >= 2,
            on_interrupt = () -> _cancel_remote(config, remote_case; verbosity)
        )
        nothing
    catch error
        if !(error isa InterruptException)
            try
                _cancel_remote(config, remote_case; verbosity)
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
        summary=first(split(sprint(showerror, execution_error), '\n'))
        throw(ErrorException(
            "$summary\nLast PSCAD diagnostics:\n$(_diagnostic_tail(console_path))" *
            "\nFull PSCAD diagnostics: $console_path" *
            "\nTransport stdout: $transport_stdout" *
            "\nTransport stderr: $transport_stderr" *
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

function _pscad_size(problem::LineParametersProblem)
    assignments = collect(Iterators.flatten(
        position.conn for position in problem.system.cables
    ))
    isempty(assignments) && throw(ArgumentError(
        "PSCAD benchmark requires at least one explicit terminal",
    ))
    any(iszero, assignments) && throw(ArgumentError(
        "PSCAD benchmark does not permit conductor elimination",
    ))
    all(>(0), assignments) || throw(ArgumentError(
        "PSCAD benchmark phase assignments must be positive",
    ))
    length(unique(assignments)) == length(assignments) || throw(ArgumentError(
        "PSCAD benchmark does not permit bundled terminals",
    ))
    sort(assignments) == collect(1:length(assignments)) || throw(ArgumentError(
        "PSCAD benchmark phase assignments must be contiguous from 1",
    ))
    return (length(assignments), length(assignments), length(problem.frequencies))
end

function _pscad_root(problem::LineParametersProblem)
    return joinpath(WORK_ROOT, "pscad", problem.system.system_id, "reference")
end

function _pscad_basis(parameters, ::LineParametersProblem, ::Val{:pul})
    parameters
end

function _pscad_basis(
        parameters::LineParameters,
        problem::LineParametersProblem,
        ::Val{:total}
)
    return LineParameters(
        PhaseDomain,
        Z(parameters) .* problem.system.line_length,
        Y(parameters) .* problem.system.line_length,
        frequencies(parameters);
        basis = :total
    )
end

function compute(
        problem::LineParametersProblem,
        formulation::PSCADFormulation;
        options::NamedTuple = (;)
)
    execution_options = computation_options(Val(PSCADFormulation), options)
    config = execution_options.remote
    root = _pscad_root(problem)
    isdir(root) && rm(root; recursive = true)
    mkpath(root)
    @info "Exporting PSCAD benchmark project" system = problem.system.system_id
    project = export_data(
        :pscad,
        problem.system,
        problem.earth_props;
        file_name = joinpath(root, "generated.pscx")
    )
    project isa AbstractString && isfile(project) || throw(ArgumentError(
        "the PSCAD exporter did not create the benchmark project",
    ))
    staged = joinpath(root, "generated.pscx")
    project == staged || cp(project, staged; force = true)
    @info "Computing PSCAD line parameters" system = problem.system.system_id
    execution = run_remote_pscad(
        config,
        staged,
        joinpath(root, "outputs"),
        formulation,
        problem.frequencies;
        output_stem = execution_options.output_stem,
        verbosity = verbosity(execution_options, :PSCAD)
    )
    write(joinpath(root, "pscad-version.txt"), config.pscad_version)
    parameters = try
        read_pscad_result(
            execution.output_dir,
            problem.frequencies,
            _pscad_size(problem)
        )
    catch error
        console_path=execution.console_path
        throw(ErrorException(
            "PSCAD result validation failed: $(sprint(showerror, error))" *
            "\nLast PSCAD diagnostics:\n$(_diagnostic_tail(console_path))" *
            "\nFull PSCAD diagnostics: $console_path",
        ))
    end
    return _pscad_basis(parameters, problem, execution_options.output_basis)
end

function benchmark_metadata(
        problem::LineParametersProblem,
        ::PSCADFormulation
)
    root = _pscad_root(problem)
    output = joinpath(root, "outputs")
    timing_path = joinpath(output, "timing.txt")
    version_path = joinpath(root, "pscad-version.txt")
    isfile(timing_path) || throw(ArgumentError("PSCAD timing is missing: $timing_path"))
    isfile(version_path) || throw(ArgumentError("PSCAD version is missing: $version_path"))
    return (
        backend = :pscad,
        elapsed_seconds = parse(Float64, strip(read(timing_path, String))),
        exit_code = 0,
        stdout_path = joinpath(output, "stdout.txt"),
        stderr_path = joinpath(output, "stderr.txt"),
        console_path = joinpath(output, "pscad-console.txt"),
        output_dir = output,
        pscad_version = strip(read(version_path, String))
    )
end
