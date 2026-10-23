"""
    LineCableModelsPSCADExt

Run native PSCAD line-parameter computations through PSCAD's automation server with
`mhi.pscad`, which PythonCall loads. The caller connects to the PSCAD that runs on the
station, loads the staged project from the shared directory and reads the run's files from
its own side of that directory.
"""
module LineCableModelsPSCADExt

using LineCableModels: validate
using LineCableModels.PSCAD: PSCAD, PSCAD_VERSION, RemoteConfig
import LineCableModels.PSCAD: identify, run_remote_pscad
using PythonCall: Py, PyException, pyconvert, pydict, pyimport, pystr
import TOML

# `app._load` takes the station's path verbatim. `app.load` would make it absolute on the
# caller. This version of `mhi.pscad` has the remote method.
const MHI_PSCAD_VERSION = "3.1.2"
const PSCAD_TIMING_SCOPE = "PSCAD compile call; excludes output-readiness wait and transfer"
# The detailed matrix files of a line-constants solve, each with its name in a run's outputs.
const DETAILED_OUTPUTS = ("_zm.out" => "result_zm.out", "_zp.out" => "result_zp.out",
    "_ym.out" => "result_ym.out", "_yp.out" => "result_yp.out")

_string(value) = pyconvert(String, pystr(value))
_components(value) = pyconvert(Vector{Py}, value)
_parameters(component) = pyconvert(Dict{String, Any}, component.parameters())

_same_value(observed::Number, requested::Number) = observed == requested
_same_value(observed, requested) = string(observed) == string(requested)

function _set!(component, name::AbstractString, value; readback = value)
    haskey(_parameters(component), name) || throw(KeyError(
        "PSCAD component $(_string(component.defn_name)) has no field $name",
    ))
    component.parameters(; Symbol(name) => value)
    observed = _parameters(component)[name]
    _same_value(observed, readback) || throw(ArgumentError(
        "PSCAD field $name rejected $(repr(value)); expected readback " *
        "$(repr(readback)), found $(repr(observed))",
    ))
    return value
end

# A conductor-elimination field that retains its conductor.
_disabled(value::Number) = iszero(value)
_disabled(value) = uppercase(strip(string(value))) in
                   ("0", "NO", "FALSE", "DISABLE", "DISABLED", "RETAIN", "NONE")

# A Python error ends the check of the forward. Any other error, an interrupt included,
# propagates.
_python_error(::PyException) = nothing
_python_error(_) = rethrow()

# Ends a child process of the run by its handle: SIGTERM, then SIGKILL (9) when the process
# has not exited `timeout` [s] later.
function _end(process, timeout::Real)
    process_running(process) && kill(process)
    timedwait(() -> process_exited(process), timeout) === :ok || kill(process, 9)
    wait(process)
    return nothing
end

# Whether PSCAD listens on `config.port` on the station, by the station's own listing of its
# TCP sockets. Through the forward, the delays of a connection to a dead port and to a
# listening one overlap.
function _listening(config::RemoteConfig)
    output, errors = IOBuffer(), IOBuffer()
    listing = run(pipeline(ignorestatus(Cmd([config.ssh..., config.host, "netstat -ano -p TCP"]));
        stdin = devnull, stdout = output, stderr = errors); wait = false)
    timedwait(() -> process_exited(listing), config.timeout) === :ok ||
        _end(listing, config.timeout)
    wait(listing)
    success(listing) || throw(ErrorException("the listing of the TCP sockets of " *
        "$(config.host) failed with exit code $(listing.exitcode) and signal " *
        "$(listing.termsignal): $(String(take!(errors)))"))
    addresses = ("127.0.0.1:$(config.port)", "0.0.0.0:$(config.port)")
    return any(eachline(IOBuffer(take!(output)))) do line
        fields = split(line)
        length(fields) >= 4 && fields[1] == "TCP" && fields[2] in addresses &&
            fields[4] == "LISTENING"
    end
end

# The forward of PSCAD's automation port on the station to the caller's `port`, open once the
# caller's loopback accepts a connection on `port`.
function _forward(config::RemoteConfig, port::Integer)
    errors = IOBuffer()
    forward = run(pipeline(Cmd([config.ssh..., "-N", "-o", "ExitOnForwardFailure=yes", "-L",
            "127.0.0.1:$port:127.0.0.1:$(config.port)", config.host]);
        stdin = devnull, stdout = devnull, stderr = errors); wait = false)
    socket = pyimport("socket")
    try
        deadline = time() + config.timeout
        while true
            accepted = try
                socket.create_connection(("127.0.0.1", port), 1).close()
                true
            catch error
                _python_error(error)
                false
            end
            accepted && break
            time() < deadline || _end(forward, config.timeout)
            process_exited(forward) && (wait(forward); throw(ErrorException(
                "the forward of port $(config.port) of $(config.host) did not open within " *
                "$(config.timeout) seconds: $(String(take!(errors)))")))
            sleep(0.1)
        end
    catch
        _end(forward, config.timeout)
        rethrow()
    end
    return forward
end

# PSCAD started over SSH in the station's session 0, with its working directory at
# `shared_root`, once the station lists it as listening on `config.port`. The process ends
# with its SSH command.
function _start(config::RemoteConfig)
    errors = IOBuffer()
    command = "cmd /c \"cd /d $(config.shared_root) && \"$(config.pscad_executable)\" " *
              "/startup:au /port:$(config.port) /silence:true /splash:false\""
    starting = run(pipeline(Cmd([config.ssh..., config.host, command]);
        stdin = devnull, stdout = devnull, stderr = errors); wait = false)
    try
        deadline = time() + config.timeout
        while !_listening(config)
            time() < deadline || _end(starting, config.timeout)
            process_exited(starting) && (wait(starting); throw(ErrorException(
                "PSCAD did not start on $(config.host) with $(config.pscad_executable) in " *
                "$(config.shared_root) and serve automation on port $(config.port) within " *
                "$(config.timeout) seconds: $(String(take!(errors)))")))
            sleep(2)
        end
    catch
        _end(starting, config.timeout)
        rethrow()
    end
    return starting
end

# The application of the PSCAD that serves automation through the forward on the caller's
# `port`, licensed and of the required version.
function _connect(config::RemoteConfig, port::Integer)
    address = "$(config.host):$(config.port)"
    app = pyimport("mhi.pscad").connect(host = "127.0.0.1", port = port, timeout = 30)
    pyconvert(Bool, app.licensed()) || throw(ErrorException(
        "PSCAD at $address holds no license"))
    version = _string(app.version)
    version == PSCAD_VERSION || throw(ErrorException(
        "PSCAD at $address is version $version; native runs require PSCAD $PSCAD_VERSION"))
    return app
end

# Calls `f(port)` with the forward of PSCAD's automation port on a free `port` of the caller's
# loopback. The forward ends on every exit.
function _forwarded(f, config::RemoteConfig)
    probe = pyimport("socket").socket()
    port = try
        probe.bind(("127.0.0.1", 0))
        pyconvert(Int, probe.getsockname()[1])
    finally
        probe.close()
    end
    forward = _forward(config, port)
    try
        return f(port)
    finally
        _end(forward, config.timeout)
    end
end

# Calls `f(app)` on the PSCAD of the station. A PSCAD that the station lists as listening on
# `config.port` is borrowed: it serves automation in the station's desktop session, and a modal
# dialog of it blocks every automation call until a person answers it. Its `silence` flag is
# set, and returns to its previous value on every exit, an error or an interrupt included. It
# keeps running when the run ends. Otherwise the run starts PSCAD with its dialogs suppressed
# and quits it on every exit. That PSCAD's SSH command ends once PSCAD exits, or at the latest
# after `config.timeout`.
function _connected(f, config::RemoteConfig)
    installed = _string(pyimport("importlib.metadata").version("mhi.pscad"))
    installed == MHI_PSCAD_VERSION || throw(ErrorException(
        "the PSCAD backend requires mhi.pscad $MHI_PSCAD_VERSION; found $installed"))
    if _listening(config)
        return _forwarded(config) do port
            app = _connect(config, port)
            silence = pyconvert(Bool, app.silence)
            app.silence = true
            try
                return f(app)
            finally
                app.silence = silence
            end
        end
    end
    starting = _start(config)
    try
        return _forwarded(config) do port
            app = _connect(config, port)
            try
                return f(app)
            finally
                try
                    app.quit()
                catch error
                    @warn "PSCAD did not quit" exception = (error, catch_backtrace())
                end
                timedwait(() -> process_exited(starting), config.timeout)
            end
        end
    finally
        _end(starting, config.timeout)
    end
end

# A PSCAD in session 0 cannot read its settings through `app.settings()`. The remote method
# `app._settings` reads the selected line-constants program verbatim.
function _identity(app)
    metadata = pyimport("importlib.metadata")
    return Dict{String, String}(
        "schema" => "2",
        "version" => _string(app.version),
        "automation_version" => _string(metadata.version("mhi.pscad")),
        "common_version" => _string(metadata.version("mhi.common")),
        "line_constants" => _string(app._settings(pydict())["file_lcp"]))
end

identify(config::RemoteConfig) = _connected(app -> validate(_identity(app), config), config)

# Each step of a run writes a line of the run's log and a debug record of the PSCAD module,
# which a run's verbosity shows at level 2 of `PSCAD`.
function _report(console::IO, message::AbstractString)
    println(console, "[", round(time(); digits = 3), "] ", message)
    flush(console)
    @debug message _module = PSCAD
    return nothing
end

# An interrupt stays an interrupt. Any other failure becomes the run's error.
_failure(::InterruptException, _) = rethrow()
_failure(_, message) = throw(ErrorException(message))

# The station's path of a file under the shared directory.
_remote_path(root::AbstractString, parts::AbstractString...) =
    join((rstrip(root, ('\\', '/')), parts...), '\\')

# The rows whose first two fields are numbers: the frequency samples of a detailed file.
function _data_rows(path::AbstractString)
    return count(eachline(path)) do line
        fields = split(line)
        length(fields) >= 2 &&
            all(field -> isfinite(something(tryparse(Float64, field), NaN)), fields[1:2])
    end
end

"""
Wait at most `timeout` \\[s\\] for `path` to hold `rows` data rows, unchanged between two
checks `poll_interval` \\[s\\] apart: PSCAD writes the detailed files while it solves.
"""
function _wait_output(path::AbstractString, rows::Integer; timeout::Real,
        poll_interval::Real = 0.1)
    deadline = time() + timeout
    # The size and the data rows of the file, both zero before PSCAD creates it.
    previous = (0, 0)
    while time() < deadline
        current = isfile(path) ? (filesize(path), _data_rows(path)) : (0, 0)
        last(current) <= rows || throw(ArgumentError(
            "PSCAD output $path contains $(last(current)) rows; expected $rows"))
        current == previous && last(current) == rows && return path
        previous = current
        sleep(poll_interval)
    end
    throw(ArgumentError("PSCAD output $path did not become complete within $timeout " *
                        "seconds; last observed $(last(previous)) of $rows rows. PSCAD writes " *
                        "the detailed files into its working directory, which must be the " *
                        "configuration's `shared_root`. A borrowed PSCAD must have been " *
                        "started there."))
end

function run_remote_pscad(config::RemoteConfig, project::AbstractString,
        output::AbstractString, settings::AbstractDict, frequencies::AbstractVector;
        line_name::AbstractString, expected_identity = nothing)
    isdir(output) && !isempty(readdir(output)) && throw(ArgumentError(
        "PSCAD output directory is not empty; select a new run directory: $output"))
    parts = splitpath(relpath(abspath(project), config.local_root))
    first(parts) == ".." && throw(ArgumentError(
        "PSCAD project must be inside remote.local_root: $project"))
    return _connected(config) do app
        _run(app, config, project, _remote_path(config.shared_root, parts...), output,
            settings, frequencies, line_name, expected_identity)
    end
end

# The run on a connected application. It checks the identity before it loads the project from
# the shared directory. After the settings, it solves the line constants and collects the
# detailed files. Every exit unloads the project that the run loaded.
function _run(app, config::RemoteConfig, project, station, output, settings, frequencies,
        line_name, expected_identity)
    mkpath(output)
    console_path = joinpath(output, "pscad-console.txt")
    console = open(console_path, "w")
    try
        _report(console, "Connected to PSCAD $PSCAD_VERSION at $(config.host):$(config.port)")
        identity = validate(_identity(app), config)
        something(expected_identity, identity) == identity || throw(ErrorException(
            "PSCAD installation does not match the expected solver identity"))
        name = first(splitext(basename(project)))
        # Whether PSCAD holds a project of the run's name, and PSCAD's messages about a project.
        held() = any(entry -> _string(entry["name"]) == name, app.projects())
        messages(loaded) = join((_string(message) for message in loaded.messages()), '\n')
        held() && throw(ErrorException(
            "PSCAD at $(config.host):$(config.port) already holds a project named $name; " *
            "unload it in PSCAD before the run"))
        phase = "project load"
        try
            _report(console, "Loading project $name from $station")
            app._load(station)
            loaded = app.project(name)
            lines = vcat(_components(loaded.find_all("TLine")),
                _components(loaded.find_all("Cable")))
            length(lines) == 1 || throw(ArgumentError(
                "generated PSCAD project must expose exactly one line-data row; " *
                "found $(length(lines))"))
            line = only(lines)
            _set!(line, "Name", line_name)
            _report(console,
                "Found line-data row $(_string(line.defn_name)); line name is $line_name")
            # The rows of the line's canvas, by definition, in one pass over its components.
            rows = Dict(definition => Py[] for definition in
                        ("master:line_frephase_options", "master:line_ground", "master:cable_coax"))
            for component in _components(loaded.canvas(
                    last(split(_string(line.defn_name), ':'; limit = 2))).components())
                definition = lowercase(_string(component.defn_name))
                haskey(rows, definition) && push!(rows[definition], component)
            end
            phase = "input configuration"
            for definition in ("master:line_frephase_options", "master:line_ground")
                length(rows[definition]) == 1 || throw(ArgumentError(
                    "PSCAD project exposes $(length(rows[definition])) components named " *
                    definition))
            end
            frequency = only(rows["master:line_frephase_options"])
            ground = only(rows["master:line_ground"])
            _set!(frequency, "FS", Float64(first(frequencies)))
            _set!(frequency, "FE", Float64(last(frequencies)))
            _set!(frequency, "Numf", length(frequencies) - 1)
            _set!(frequency, "Output", 1; readback = "YES")
            readbacks = Dict{String, Any}()
            for (row, component) in
                (("ground", ground), ("frequency", frequency), ("configuration", line))
                observed = Dict{String, Any}()
                for (field, setting) in settings[row]
                    _set!(component, field, setting["value"]; readback = setting["readback"])
                    observed[field] = _parameters(component)[field]
                end
                readbacks[row] = observed
            end
            open(io -> TOML.print(io, readbacks; sorted = true),
                joinpath(output, "native-settings.toml"), "w")
            _report(console, "Applied frequency range $(first(frequencies)) Hz to " *
                             "$(last(frequencies)) Hz with $(length(frequencies)) samples")
            # The generated PSCAD project retains every requested cable port.
            cables = rows["master:cable_coax"]
            isempty(cables) && throw(ArgumentError(
                "generated PSCAD project exposes no master:Cable_Coax components"))
            for cable in cables, (field, value) in _parameters(cable)
                occursin(r"^elim\d+$", field) && !_disabled(value) && throw(ArgumentError(
                    "PSCAD field $field requests conductor elimination; " *
                    "requested ports must be retained"))
            end
            _report(console, "Verified that all cable terminals are retained")
            loaded.save()
            println(console, "PSCAD messages before computation:\n", messages(loaded))
            phase = "line-constants computation"
            _report(console, "Starting PSCAD line-constants computation")
            started = time_ns()
            line.compile()
            elapsed = (time_ns() - started) * 1e-9
            _report(console,
                "PSCAD compile call returned in $(round(elapsed; digits = 3)) seconds")
            println(console, "PSCAD messages after the compile call:\n", messages(loaded))
            write(joinpath(output, "timing.txt"), string(elapsed))
            phase = "detailed-output completion"
            samples = length(frequencies)
            _report(console, "Waiting for $samples rows in each detailed PSCAD matrix file")
            for (suffix, _) in DETAILED_OUTPUTS
                _wait_output(joinpath(config.local_root, line_name * suffix), samples;
                    timeout = config.timeout)
            end
            # The line's files are its detailed files and the line-constants program's log.
            destinations = Dict(line_name * suffix => name for (suffix, name) in DETAILED_OUTPUTS)
            for file in readdir(config.local_root)
                startswith(file, line_name * "_") || startswith(file, line_name * ".") ||
                    continue
                mv(joinpath(config.local_root, file),
                    joinpath(output, get(destinations, file, file)))
            end
            _report(console, "Collected detailed PSCAD Z and Y outputs")
            open(io -> TOML.print(io, identity; sorted = true),
                joinpath(output, "solver.toml"), "w")
            @debug("PSCAD frequency scan completed", _module = PSCAD,
                compile_call_seconds = elapsed, timing_scope = PSCAD_TIMING_SCOPE)
            return (elapsed_seconds = elapsed, elapsed_scope = PSCAD_TIMING_SCOPE,
                console_path, output_dir = String(output))
        catch error
            # PSCAD's messages and the project's output. A failed read takes their place, and
            # the failure stays the error.
            diagnostics = try
                held_project = app.project(name)
                "PSCAD messages:\n$(messages(held_project))\n" *
                "PSCAD project output:\n$(strip(_string(held_project.output())))"
            catch reading
                "PSCAD diagnostics are unavailable: $(sprint(showerror, reading))"
            end
            message = "PSCAD run failed during $phase: $(sprint(showerror, error))\n" *
                      "$diagnostics\nFull PSCAD log: $console_path"
            println(console, message)
            _failure(error, message)
        finally
            try
                held() && app.project(name).unload()
            catch error
                println(console, "PSCAD project unload failed: ", sprint(showerror, error))
            end
        end
    finally
        close(console)
    end
end

end # module LineCableModelsPSCADExt
