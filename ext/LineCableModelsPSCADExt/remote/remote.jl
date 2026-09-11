function _powershell_argv(powershell::AbstractString)
    command="\$ProgressPreference='SilentlyContinue'; $(String(powershell))"
    utf16 = htol.(transcode(UInt16, command))
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
    # A caller may supply a transport method after the package was compiled.
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
    receiver=LineCableModels.progress_receiver()
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
        progress_line=_remote_progress(receiver,line)
        progress_line && continue # Disposable observations are not an event archive.
        println(output_buffer, line)
        output_file === nothing || println(output_file, line)
        output_file === nothing || flush(output_file)
        if stream
            LineCableModels.with_progress_output() do
                println(stdout, line)
                flush(stdout)
            end
        end
    end
    error_task = @async for line in eachline(errors)
        println(error_buffer, line)
        error_file === nothing || println(error_file, line)
        error_file === nothing || flush(error_file)
        if stream
            LineCableModels.with_progress_output() do
                println(stderr, line)
                flush(stderr)
            end
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

function _remote_progress(receiver, line)
    startswith(line,"LCM_PROGRESS_V1\t") || return false
    fields=split(line,'\t')
    length(fields)==3 || return true
    if fields[2]=="stage" && fields[3] in
            ("launching","loading","configuring","compiling","waiting_outputs","transferring","validating")
        LineCableModels.report_progress(receiver,(backend=:pscad,stage=Symbol(fields[3])))
    elseif fields[2]=="heartbeat"
        LineCableModels.report_progress(receiver,(backend=:pscad,heartbeat_unix_seconds=time()))
    end
    return true
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
            "-FrequencyStart $(_ps_quote(string(first(frequencies_value))))",
            "-FrequencyEnd $(_ps_quote(string(last(frequencies_value))))",
            "-FrequencyIncrements $(_ps_quote(string(increments)))",
            "-PSCADVersion $(_ps_quote(config.pscad_version))",
            "-Verbosity $(_ps_quote(string(verbosity)))",
            "-TimeoutSeconds $(_ps_quote(string(config.timeout_seconds)))",
            LineCableModels.progress_receiver() === nothing ? "" : "-TrackProgress"
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

# Keep remotely executed code paired with the loaded Julia adapter. Later
# invocations in a long campaign must not pick up working-tree edits mid-run.
const PSCAD_REMOTE_SOURCES = Dict(name => let
                                      path = joinpath(@__DIR__, name)
                                      Base.include_dependency(path)
                                      read(path, String)
                                  end
for name in ("Project.toml", "Manifest.toml", "files.jl",
    "runner.jl", "supervisor.ps1", "identity.py"))

"""
    identify(config::RemoteConfig)

Read the remote PSCAD installation identity without launching a simulation.

# Arguments

- `config`: Station connection and selected PSCAD installation.

# Returns

- A string dictionary containing application, line-constants executable and
  master-library paths and SHA-256 digests, automation versions, and the
  selected line-constants implementation. A temporary application instance
  reads the station settings and closes without loading a project.

# Notes

Passing this record as the `solver_identity` computation option pins a campaign
to that installation. Each computation verifies the station again; supplying a
record does not bypass the check.
"""
function identify(config::RemoteConfig)
    code = PSCAD_REMOTE_SOURCES["identity.py"] * "\nimport json\n" *
           "for key, value in identify(" * repr(config.pscad_version) * ").items():\n" *
           "    print(json.dumps(key) + ' = ' + json.dumps(value))\n"
    # Use the existing shared work directory; embedding a whole script inside
    # an encoded PowerShell command exceeds Windows' command-line limit.
    directory = mktempdir(mkpath(config.local_root); prefix = "pscad-identity-")
    result = try
        write(joinpath(directory, "identify.py"), code)
        remote = _remote_path(config.shared_root, basename(directory), "identify.py")
        command = "& " * _ps_quote(config.python_executable) * " " * _ps_quote(remote) *
                  "; if (\$LASTEXITCODE -ne 0) { exit \$LASTEXITCODE }"
        TOML.parse(_run_remote(config, command))
    finally
        rm(directory; recursive = true)
    end
    get(result, "version", nothing) == config.pscad_version || throw(ArgumentError(
        "PSCAD station did not return the requested solver identity"))
    get(result, "schema", nothing) == "1" &&
    all(name -> occursin(r"^[0-9a-f]{64}$", get(result, name * "_sha256", "")),
        ("pscad", "line_constants", "master_library")) || throw(ArgumentError(
        "PSCAD station returned an incomplete solver identity"))
    return Dict{String, String}(result)
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
    toolkit_stage = joinpath(variant_root, "toolkit")
    isdir(toolkit_stage) && rm(toolkit_stage; recursive = true)
    mkpath(toolkit_stage)
    for (name, source) in PSCAD_REMOTE_SOURCES
        write(joinpath(toolkit_stage, name), source)
    end
    return toolkit_stage
end

function _diagnostic_tail(path::AbstractString; count::Integer = 12)
    isfile(path) || return "PSCAD produced no diagnostic log."
    lines=filter(!isempty, strip.(readlines(path)))
    isempty(lines) && return "PSCAD diagnostic log is empty."
    return join(last(lines, min(count, length(lines))), '\n')
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
    LineCableModels.performance_sample_active() && (verbosity=0)
    LineCableModels.report_progress(LineCableModels.progress_receiver(),(backend=:pscad,stage=:staging))
    _validate_frequencies(frequencies_value)
    isdir(local_output) && !isempty(readdir(local_output)) &&
        throw(ArgumentError(
            "PSCAD output directory is not empty; select a new run directory: $local_output"))
    mkpath(local_output)
    relative = relpath(dirname(abspath(local_output)), config.local_root)
    work_parts = splitpath(relative)
    first(work_parts) == ".." && throw(ArgumentError(
        "PSCAD native run directory must be inside remote.local_root"))
    variant = last(work_parts)
    _stage_toolkit(local_project, local_output)
    shared_case = _remote_path(config.shared_root, work_parts...)
    remote_case = _remote_path(config.remote_root, work_parts..., output_stem)
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
    verbosity >= 1 && @info "Executing PSCAD frequency scan" host=config.host variant formulation=_formulation_label(formulation) frequencies=length(frequencies_value) timeout_seconds=config.timeout_seconds
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
    verbosity >= 1 && @info "Checking PSCAD outputs" host=config.host variant destination=local_output
    LineCableModels.report_progress(LineCableModels.progress_receiver(),(backend=:pscad,stage=:validating))
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
    verbosity >= 1 && @info "PSCAD frequency scan completed" host=config.host variant compile_call_seconds=elapsed timing_scope=PSCAD_TIMING_SCOPE
    return (
        elapsed_seconds = elapsed,
        elapsed_scope = PSCAD_TIMING_SCOPE,
        exit_code = 0,
        stdout_path,
        stderr_path,
        console_path = joinpath(local_output, "pscad-console.txt"),
        output_dir = String(local_output)
    )
end
