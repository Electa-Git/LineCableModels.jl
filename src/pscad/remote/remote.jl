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

function remote_command(::Val{:local}, ::RemoteConfig, powershell::AbstractString)
    return Cmd(_powershell_argv(powershell))
end

function remote_command(::Val{:command}, config::RemoteConfig, powershell::AbstractString)
    arguments = [value == "{host}" ? config.host : value for value in config.command]
    return Cmd(vcat(arguments, _powershell_argv(powershell)))
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
        on_interrupt::Function = () -> nothing,
        timeout::Real = config.timeout
)
    isfinite(timeout) && timeout > 0 || throw(ArgumentError(
        "PSCAD transport timeout must be positive and finite"))
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
        finished = timedwait(() -> istaskfailed(output_task) || istaskfailed(error_task) ||
            (process_exited(process) && istaskdone(output_task) && istaskdone(error_task)),
            timeout; pollint = min(0.1, timeout / 10))
        finished === :ok || throw(ErrorException(
            "PSCAD transport exceeded its timeout of $timeout seconds"))
        istaskfailed(output_task) && wait(output_task)
        istaskfailed(error_task) && wait(error_task)
        wait(process)
    catch error
        interrupted = error
        # SIGKILL (9): a transport must not keep the caller waiting in its own
        # signal handler after timeout or interruption.
        process_running(process) && kill(process, 9)
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
        catch wait_error
            @warn "PSCAD transport wait failed during cleanup" exception = (
                wait_error, catch_backtrace())
        end
    finally
        if interrupted !== nothing
            close(output)
            close(errors)
        end
        for (stream_name, reader) in ((:stdout, output_task), (:stderr, error_task))
            try
                wait(reader)
            catch reader_error
                if interrupted === nothing
                    interrupted = reader_error
                elseif !(interrupted isa TaskFailedException && interrupted.task === reader)
                    @warn "PSCAD output reader failed during cleanup" stream = stream_name exception = (
                        reader_error, catch_backtrace())
                end
            end
        end
        for (stream_name, file) in ((:stdout, output_file), (:stderr, error_file))
            file === nothing && continue
            try
                close(file)
            catch close_error
                if interrupted === nothing
                    interrupted = close_error
                else
                    @warn "PSCAD output log could not be closed" stream = stream_name exception = (
                        close_error, catch_backtrace())
                end
            end
        end
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

function _remote_protocol(formulation::PSCADFormulation, frequencies_value)
    label = only(description([formulation]; roles = [:none]))
    return (; mode = "constants", label,
        start_message = "Executing PSCAD frequency scan",
        start_fields = (; formulation = label, frequencies = length(frequencies_value)),
        completed_message = "PSCAD frequency scan completed",
        elapsed_field = :compile_call_seconds, elapsed_scope = PSCAD_TIMING_SCOPE)
end

function _remote_protocol(formulation::PSCADTransientFormulation, frequencies_value)
    return (; mode = "transient", label = string(formulation.model),
        start_message = "Executing PSCAD transient run",
        start_fields = (; model = formulation.model),
        completed_message = "PSCAD transient run completed",
        elapsed_field = :run_call_seconds, elapsed_scope = "project.run()")
end

_pscad_port(::PSCADFormulation, ::RemoteConfig) = 0
_pscad_port(::PSCADTransientFormulation, config::RemoteConfig) = something(config.pscad_port, 0)

function _supervisor_command(
        config::RemoteConfig,
        shared_case::AbstractString,
        remote_case::AbstractString,
        project_name::AbstractString,
        formulation::Union{PSCADFormulation, PSCADTransientFormulation},
        frequencies_value::AbstractVector;
        output_stem::AbstractString,
        verbosity::Integer = 0
)
    protocol = _remote_protocol(formulation, frequencies_value)
    increments = length(frequencies_value) - 1
    shared_supervisor = _remote_path(shared_case, "toolkit", "supervisor.ps1")
    return join(
        (
            "& $(_ps_quote(shared_supervisor))",
            "-SharedCase $(_ps_quote(shared_case))",
            "-LocalCase $(_ps_quote(remote_case))",
            "-RemoteRoot $(_ps_quote(config.remote_root))",
            "-Julia $(_ps_quote(config.julia_executable))",
            "-Python $(_ps_quote(config.python_executable))",
            "-ProjectName $(_ps_quote(project_name))",
            "-OutputStem $(_ps_quote(output_stem))",
            "-Mode $(_ps_quote(protocol.mode))",
            "-Formulation $(_ps_quote(protocol.label))",
            "-FrequencyStart $(_ps_quote(string(first(frequencies_value))))",
            "-FrequencyEnd $(_ps_quote(string(last(frequencies_value))))",
            "-FrequencyIncrements $(_ps_quote(string(increments)))",
            "-PSCADVersion $(_ps_quote(config.pscad_version))",
            "-PSCADPort $(_pscad_port(formulation, config))",
            "-Verbosity $(_ps_quote(string(verbosity)))",
            "-TimeoutSeconds $(_ps_quote(string(config.timeout)))"
        ),
        ' ')
end

function _cancel_command(shared_toolkit::AbstractString, remote_case::AbstractString)
    lifecycle = _remote_path(shared_toolkit, "lifecycle.ps1")
    return ". $(_ps_quote(lifecycle)); Stop-PSCADCase -LocalCase $(_ps_quote(remote_case))"
end

function _cancel_remote(
        config::RemoteConfig,
        shared_toolkit::AbstractString,
        remote_case::AbstractString;
        verbosity::Integer = 0
)
    _run_remote(
        config, _cancel_command(shared_toolkit, remote_case); stream = verbosity >= 2,
        timeout = min(config.timeout, 30))
    return nothing
end

# Keep remotely executed code paired with the loaded Julia adapter. Later
# invocations must not pick up working-tree edits mid-run.
const PSCAD_REMOTE_SOURCES = Dict(name => let
                                      path = joinpath(@__DIR__, name)
                                      Base.include_dependency(path)
                                      read(path, String)
                                  end
for name in ("Project.toml", "Manifest.toml", "files.jl",
    "runner.jl", "supervisor.ps1", "lifecycle.ps1", "identity.py", "probe.py"))

"""
$(TYPEDSIGNATURES)

Read the remote PSCAD installation identity without launching a simulation.

# Arguments

- `config`: station connection and selected PSCAD installation.

# Returns

- A string dictionary containing application, line-constants executable and
  master-library paths and SHA-256 digests, automation versions, and the
  selected line-constants implementation. A temporary application instance
  reads the station settings and closes without loading a project.

# Notes

Passing this record as the `solver_identity` computation option requires that
installation. Fresh computations verify identity inside their own application
instance. Completed-run reuse verifies the current station before acceptance.
"""
function identify(config::RemoteConfig)
    code = PSCAD_REMOTE_SOURCES["identity.py"] * "\nimport json\n" *
           "for key, value in identify(" * repr(config.pscad_version) * ").items():\n" *
           "    print(json.dumps(key) + ' = ' + json.dumps(value))\n"
    # Use the existing shared work directory. Embedding a whole script inside
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
    return Dict{String, String}(validate(result, config))
end

# The solver identity that a PSCAD station reports for the requested configuration.
function validate(result::AbstractDict, config::RemoteConfig)
    all(pair -> first(pair) isa AbstractString && last(pair) isa AbstractString, result) ||
        throw(ArgumentError("PSCAD solver identity must contain string fields"))
    get(result, "version", nothing) == config.pscad_version || throw(ArgumentError(
        "PSCAD station did not return the requested solver identity"))
    get(result, "schema", nothing) == "1" &&
    all(name -> occursin(r"^[0-9a-f]{64}$", get(result, name * "_sha256", "")),
        ("pscad", "line_constants", "master_library")) || throw(ArgumentError(
        "PSCAD station returned an incomplete solver identity"))
    return result
end

function _stage_toolkit(local_project::AbstractString, local_output::AbstractString)
    isfile(local_project) || throw(ArgumentError(
        "local PSCAD input is missing: $local_project",
    ))
    run_directory = dirname(local_output)
    expected_project = joinpath(run_directory, "generated.pscx")
    abspath(local_project) == abspath(expected_project) || throw(ArgumentError(
        "PSCAD project must be staged as $expected_project",
    ))
    toolkit_stage = joinpath(run_directory, "toolkit")
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

function validate(directory::AbstractString, ::PSCADFormulation;
        output_stem::AbstractString)
    required = ("pscad-console.txt", "timing.txt", "result_zm.out", "result_zp.out",
        "result_ym.out", "result_yp.out")
    for name in required
        path = joinpath(directory, name)
        isfile(path) || throw(ArgumentError("required PSCAD output is missing: $path"))
        filesize(path) > 0 || throw(ArgumentError("required PSCAD output is empty: $path"))
    end
    return directory
end

function validate(directory::AbstractString, ::PSCADTransientFormulation;
        output_stem::AbstractString)
    required = ("pscad-console.txt", "timing.txt", "execution.toml", output_stem * ".inf")
    for name in required
        path = joinpath(directory, name)
        isfile(path) || throw(ArgumentError("required PSCAD output is missing: $path"))
        filesize(path) > 0 || throw(ArgumentError("required PSCAD output is empty: $path"))
    end
    read_pscad_transient(directory; stem = output_stem)
    return directory
end

function _supervisor_errors(path::AbstractString)
    isfile(path) || return ""
    text = read(path, String)
    # Windows PowerShell serializes stderr strings as CLIXML stream records.
    records = [record[1]
               for record in eachmatch(r"<S\b[^>]*\bS=\"Error\"[^>]*>(.*?)</S>"s, text)]
    isempty(records) && return occursin("#< CLIXML", text) ? "" : strip(text)
    return strip(join(
        (rstrip(replace(record,
             "_x000D__x000A_" => "\n", "_x000D_" => "\n", "_x000A_" => "\n",
             "&lt;" => "<", "&gt;" => ">", "&quot;" => "\"", "&apos;" => "'", "&amp;" => "&"))
        for record in records),
        '\n'))
end

function run_remote_pscad(
        config::RemoteConfig,
        local_project::AbstractString,
        local_output::AbstractString,
        formulation::Union{PSCADFormulation, PSCADTransientFormulation},
        frequencies_value::AbstractVector;
        output_stem::AbstractString,
        verbosity::Integer = 0
)
    verbosity in 0:2 || throw(ArgumentError("PSCAD verbosity must be 0, 1, or 2"))
    isdir(local_output) && !isempty(readdir(local_output)) &&
        throw(ArgumentError(
            "PSCAD output directory is not empty; select a new run directory: $local_output"))
    mkpath(local_output)
    relative = relpath(dirname(abspath(local_output)), config.local_root)
    work_parts = splitpath(relative)
    first(work_parts) == ".." && throw(ArgumentError(
        "PSCAD native run directory must be inside remote.local_root"))
    run_id = last(work_parts)
    _stage_toolkit(local_project, local_output)
    shared_case = _remote_path(config.shared_root, work_parts...)
    shared_toolkit = _remote_path(shared_case, "toolkit")
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
    protocol = _remote_protocol(formulation, frequencies_value)
    verbosity >= 2 && @debug protocol.start_message host=config.host run_id protocol.start_fields... timeout=config.timeout
    execution_error = try
        _run_remote(
            config,
            command;
            stdout_path = transport_stdout,
            stderr_path = transport_stderr,
            stream = verbosity >= 2,
            timeout = config.timeout + (_pscad_port(formulation, config) == 0 ? 60 : 160),
            on_interrupt = () -> _cancel_remote(config, shared_toolkit, remote_case; verbosity)
        )
        nothing
    catch error
        if !(error isa InterruptException)
            try
                _cancel_remote(config, shared_toolkit, remote_case; verbosity)
            catch cancellation_error
                @warn "Remote PSCAD cancellation could not be confirmed" host=config.host run_id exception=(
                    cancellation_error, catch_backtrace())
            end
        end
        error isa InterruptException && rethrow()
        error
    end
    if execution_error !== nothing
        failure_path = joinpath(local_output, "failure.toml")
        if isfile(failure_path)
            failure = TOML.parsefile(failure_path)
            if get(failure, "kind", nothing) == "line_constants_error"
                log_path = failure["log_path"]
                local_log = isempty(log_path) ? "" :
                    abspath(joinpath(local_output, split(replace(log_path, '\\' => '/'), '/')...))
                blocks = failure["blocks"]
                message = isempty(blocks) ? join(failure["errors"], "; ") :
                    join([(isempty(block["stage"]) ? "" : block["stage"] * " ") *
                        join(block["errors"], "; ") for block in blocks], '\n')
                throw(LineConstantsError(message, local_log))
            end
        end
        console_path = joinpath(local_output, "pscad-console.txt")
        if !isfile(console_path)
            write(
                console_path,
                "PSCAD did not produce a console log before the remote failure.\n"
            )
        end
        summary=first(split(sprint(showerror, execution_error), '\n'))
        diagnostics = _diagnostic_tail(console_path)
        supervisor_errors = _supervisor_errors(transport_stderr)
        if !isempty(supervisor_errors)
            if !occursin("PSCAD runner failed during", diagnostics)
                summary *= "\nSupervisor error:\n$supervisor_errors"
            elseif occursin("Cleanup failure:", supervisor_errors)
                diagnostics *= "\nSupervisor error:\n$supervisor_errors"
            end
        end
        throw(ErrorException(
            "$summary\nLast PSCAD diagnostics:\n$diagnostics" *
            "\nFull PSCAD diagnostics: $console_path" *
            "\nTransport stdout: $transport_stdout" *
            "\nTransport stderr: $transport_stderr" *
            "\nRemote scratch: $remote_case",
        ))
    end
    verbosity >= 2 && @debug "Checking PSCAD outputs" host=config.host run_id destination=local_output
    validate(local_output, formulation; output_stem)
    elapsed = parse(Float64, strip(read(joinpath(local_output, "timing.txt"), String)))
    elapsed_scope = protocol.elapsed_scope
    elapsed_fields = NamedTuple{(protocol.elapsed_field,)}((elapsed,))
    verbosity >= 2 && @debug protocol.completed_message host=config.host run_id elapsed_fields... timing_scope=elapsed_scope
    return (
        elapsed_seconds = elapsed,
        elapsed_scope,
        exit_code = 0,
        stdout_path,
        stderr_path,
        console_path = joinpath(local_output, "pscad-console.txt"),
        output_dir = String(local_output)
    )
end
