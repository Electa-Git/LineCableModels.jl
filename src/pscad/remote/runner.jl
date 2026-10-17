using PythonCall
using TOML

include(joinpath(@__DIR__, "files.jl"))

const MHI_PSCAD_VERSION = "3.1.2"
const OUTPUT_SUFFIXES = (
    "_zm.out" => "result_zm.out",
    "_zp.out" => "result_zp.out",
    "_ym.out" => "result_ym.out",
    "_yp.out" => "result_yp.out"
)
_string(value) = pyconvert(String, pybuiltins.str(value))
_components(value) = pyconvert(Vector{Py}, pybuiltins.list(value))

function _schematic_origins(path)
    root = pyimport("xml.etree.ElementTree").parse(String(path)).getroot()
    objects = Dict{Tuple{String, String}, NamedTuple{(:x, :y, :kind), Tuple{Int, Int, String}}}()
    for definition in root.findall("./definitions/Definition")
        name = _string(definition.get("name"))
        for canvas in definition.findall("schematic"), node in canvas.iter()
            pyis(node.get("id"), pybuiltins.None) && continue
            pyis(node.get("x"), pybuiltins.None) && continue
            objects[(name, _string(node.get("id")))] = (
                x = parse(Int, _string(node.get("x"))),
                y = parse(Int, _string(node.get("y"))),
                kind = _string(node.get("defn", node.get("classid", node.tag))))
        end
    end
    return objects
end

function _load_relocations(before, after)
    moved = String[]
    for key in sort!(collect(keys(before)))
        haskey(after, key) || continue
        a, b = before[key], after[key]
        (a.x, a.y) == (b.x, b.y) && continue
        push!(moved, "$(key[1]) id=$(key[2]) $(a.kind): $((a.x, a.y)) -> $((b.x, b.y))")
    end
    return moved
end

function _definition(component)
    return lowercase(_string(component.defn_name))
end

function _single_component(components, definition)
    matches = filter(component -> _definition(component) == lowercase(definition), components)
    length(matches) == 1 || throw(ArgumentError(
        "PSCAD project exposes $(length(matches)) components named $definition",
    ))
    return only(matches)
end

function _parameters(component)
    raw = pyconvert(Dict, component.parameters())
    return Dict(string(key) => value for (key, value) in raw)
end

_same_value(observed::Number, requested::Number) = observed == requested
_same_value(observed, requested) = string(observed) == string(requested)

function _set!(component, name::AbstractString, value; readback = value)
    haskey(_parameters(component), name) || throw(KeyError(
        "PSCAD component $(_string(component.defn_name)) has no field $name",
    ))
    component.parameters(; Dict(Symbol(name) => value)...)
    observed = _parameters(component)[name]
    _same_value(observed, readback) || throw(ArgumentError(
        "PSCAD field $name rejected $(repr(value)); expected readback " *
        "$(repr(readback)), found $(repr(observed))",
    ))
    return value
end

function _disabled(value)
    value isa Number && return iszero(value)
    return uppercase(strip(string(value))) in (
        "0", "NO", "FALSE", "DISABLE", "DISABLED", "RETAIN", "NONE")
end

function _messages(project)
    messages = try
        _components(project.messages())
    catch error
        return "Unable to read PSCAD messages: $(sprint(showerror, error))\n"
    end
    return join((_string(message) for message in messages), '\n') * "\n"
end

function _report(console::IO, verbosity::Int, level::Int, message::AbstractString)
    line = "[$(round(time(); digits = 3))] $(String(message))"
    println(console, line)
    flush(console)
    if verbosity >= level
        println(stdout, line)
        flush(stdout)
    end
    return nothing
end

function _project_diagnostics(project)
    parts=String["PSCAD messages:\n$(strip(_messages(project)))"]
    try
        output=strip(_string(project.output()))
        isempty(output) || push!(parts, "PSCAD project output:\n$output")
    catch error
        push!(parts, "Unable to read PSCAD project output: $(sprint(showerror, error))")
    end
    return join(parts, '\n')
end

function _record_diagnostics(
        console::IO,
        project,
        verbosity::Int,
        heading::AbstractString;
        error_stream::Bool = false
)
    text="$(String(heading))\n$(_project_diagnostics(project))"
    println(console, text)
    flush(console)
    if verbosity >= 2
        stream=error_stream ? stderr : stdout
        println(stream, text)
        flush(stream)
    end
    return text
end

function _build_messages(project, path, case_directory)
    messages = _components(project.messages())
    open(path, "w") do io
        for message in messages
            println(io,
                join(
                    (_string(message.status), _string(message.label),
                        _string(message.name), replace(_string(message.text), '\n' => ' ', '\r' =>
                            ' ')),
                    '\t'))
        end
    end
    if any(message -> lowercase(_string(message.status)) == "error", messages)
        output = dirname(path)
        diagnostics = String[]
        logs = setdiff(_output_matches([case_directory], ".log"),
            _output_matches([output], ".log"))
        for source in logs
            destination = joinpath(output, relpath(source, case_directory))
            try
                mkpath(dirname(destination))
                cp(source, destination; force = true)
                push!(diagnostics, "Last lines of $(relpath(source, case_directory)): " *
                    join(last(readlines(destination), 12), " | "))
            catch log_error
                push!(diagnostics, "Unable to collect $source: $(sprint(showerror, log_error))")
            end
        end
        error("PSCAD build messages contain error status; see $path\n" *
            join(diagnostics, '\n'))
    end
    return messages
end

function main(arguments)
    length(arguments) == 14 || throw(ArgumentError(
        "runner expects project, output, project name, output stem, formulation, " *
        "FS, FE, Numf, PSCAD version, " *
        "verbosity, mode, port, ownership, and scratch root",
    ))
    project_path, output, project_name, output_stem, formulation,
    fs, fe, numf, pscad_version,
    verbosity_text = arguments[1:10]
    mode, port_text, owned_text, remote_root = arguments[11:14]
    port = parse(Int, port_text)
    (port == 0 || port in 1024:65535) || throw(ArgumentError("invalid PSCAD port"))
    owned = parse(Bool, owned_text)
    port == 0 && !owned && throw(ArgumentError("a private PSCAD launch must be owned"))
    mode in ("constants", "transient") || throw(ArgumentError("invalid PSCAD mode $mode"))
    verbosity = parse(Int, verbosity_text)
    verbosity in 0:2 || throw(ArgumentError("verbosity must be 0, 1, or 2"))
    pscad_version == "5.1.0" || throw(ArgumentError("PSCAD 5.1.0 is required"))
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_]{0,19}$", output_stem) ||
        throw(ArgumentError("invalid PSCAD output stem $output_stem"))
    mkpath(output)
    app = nothing
    previous_silence = nothing
    project = nothing
    load_attempted = false
    phase = "initialization"
    input_path = joinpath(dirname(project_path), "computation.toml")
    input = TOML.parsefile(input_path)
    get(input, "schema_version", nothing) == 4 || throw(ArgumentError(
        "PSCAD runner requires the complete version-4 numerical input record"))
    identify = pyimport("runpy").run_path(joinpath(@__DIR__, "identity.py"))["identify"]
    console = open(joinpath(output, "pscad-console.txt"), "w")
    build_path = joinpath(output, "build-messages.txt")
    write(build_path, "")
    try
        metadata = pyimport("importlib.metadata")
        observed = _string(metadata.version("mhi.pscad"))
        observed == MHI_PSCAD_VERSION || throw(ArgumentError(
            "mhi.pscad $MHI_PSCAD_VERSION is required; found $observed",
        ))
        pscad = pyimport("mhi.pscad")
        phase = "PSCAD acquisition"
        if port == 0
            _report(console, verbosity, 1, "Launching PSCAD 5.1.0")
            app = cd(abspath(output)) do
                pscad.launch(; version = pscad_version, x64 = true, minimize = true,
                    splash = false, silence = true, timeout = 60)
            end
        else
            _report(console, verbosity, 1, "Connecting to PSCAD on port $port (owned=$owned)")
            connect_started = time_ns()
            app = pscad.connect(; host = "127.0.0.1", port, timeout = 30)
            connect_seconds = (time_ns() - connect_started) / 1e9
            _report(console, verbosity, 1,
                "PSCAD connect took $(round(connect_seconds; digits = 3)) s")
        end
        phase = "PSCAD dialog suppression"
        previous_silence = pyconvert(Bool, app.silence)
        app.silence = true
        _report(console, verbosity, 1, "silence: previous=$previous_silence, now=true")
        if mode == "transient"
            phase = "compiler detection"
            application_settings = try
                app.settings()
            catch capability_error
                occursin("Unable to retrieve detected software", sprint(showerror, capability_error)) ||
                    rethrow()
                error("PSCAD cannot detect compilers in this Windows session (non-interactive, e.g. SSH session 0); set pscad_port in the station configuration so runs use a PSCAD in the logged-on desktop session")
            end
            _report(console, verbosity, 1, "fortran_version=$(_string(application_settings["fortran_version"]))")
        end
        phase = "PSCAD licence and identity"
        pyconvert(Bool, app.licensed()) || error("PSCAD refused the configured license")
        _string(app.version) == pscad_version ||
            error("PSCAD launched an unexpected version")
        initial_identity = pyconvert(Dict{String, String}, identify(pscad_version, app))
        expected_identity = get(input, "expected_solver", nothing)
        expected_identity === nothing || initial_identity == expected_identity ||
            error("PSCAD installation does not match the expected solver identity")
        phase = "project load"
        for entry in _components(app.projects())
            _string(entry["name"]) == project_name || continue
            existing = app.project(project_name)
            filename = _string(existing.filename)
            root = lowercase(replace(normpath(abspath(remote_root)), '\\' => '/'))
            path = lowercase(replace(normpath(abspath(filename)), '\\' => '/'))
            startswith(path, rstrip(root, '/') * "/") ||
                error("Refusing to replace loaded PSCAD project $project_name from $filename; it is outside the station scratch root $remote_root")
            existing.unload()
            _report(console, verbosity, 1, "Unloaded stale PSCAD project $project_name from $filename")
        end
        _report(console, verbosity, 1, "Loading generated project $project_name")
        load_attempted = true
        origins = _schematic_origins(project_path)
        app.load(abspath(project_path))
        project = app.project(project_name)
        lines = vcat(
            _components(project.find_all("TLine")),
            _components(project.find_all("Cable"))
        )
        length(lines) == 1 || throw(ArgumentError(
            "generated PSCAD project must expose exactly one line-data row; found $(length(lines))",
        ))
        line = only(lines)
        mode == "constants" && _set!(line, "Name", output_stem)
        _report(
            console,
            verbosity,
            2,
            "Found line-data row $(_string(line.defn_name)); output stem is $output_stem"
        )
        canvas = try
            project.canvas(last(split(_string(line.defn_name), ':'; limit = 2)))
        catch
            line.canvas()
        end
        components = _components(canvas.components())
        phase = "input configuration"
        frequency = mode == "constants" ?
            _single_component(components, "master:line_frephase_options") :
            only(filter(component -> _definition(component) in
                ("master:line_frephase_options", "master:line_fredep_options"), components))
        ground = _single_component(components, "master:line_ground")
        if mode == "constants"
            _set!(frequency, "FS", parse(Float64, fs))
            _set!(frequency, "FE", parse(Float64, fe))
            _set!(frequency, "Numf", parse(Int, numf))
            _set!(frequency, "Output", 1; readback = "YES")
        end
        readbacks = Dict{String, Any}()
        for (name, component) in
            (("ground", ground), ("frequency", frequency), ("configuration", line))
            controls = input["native_settings"][name]
            observed = Dict{String, Any}()
            for (field, setting) in controls
                _set!(component, field, setting["value"]; readback = setting["readback"])
                observed[field] = _parameters(component)[field]
            end
            readbacks[name] = observed
        end
        if mode == "transient"
            _set!(project, "output_filename", output_stem * ".out")
            settings = _parameters(project)
            readbacks["project"] = Dict(field => settings[field] for field in
                ("time_duration", "time_step", "sample_step", "PlotType", "output_filename"))
            open(joinpath(output, "execution.toml"), "w") do io
                TOML.print(io, Dict("mode" => mode, "stem" => output_stem); sorted = true)
            end
        end
        open(joinpath(output, "native-settings.toml"), "w") do io
            TOML.print(io, readbacks; sorted = true)
        end
        _report(
            console,
            verbosity,
            2,
            "Applied frequency range $fs Hz to $fe Hz and formulation $formulation"
        )
        # The generated PSCAD project retains every requested cable port.
        cables = filter(
            component -> _definition(component) == "master:cable_coax",
            components
        )
        isempty(cables) && throw(ArgumentError(
            "generated PSCAD project exposes no master:Cable_Coax components",
        ))
        for cable in cables
            parameters = _parameters(cable)
            elimination_fields = sort!(filter(
                field -> occursin(r"^elim\d+$", field),
                collect(keys(parameters))
            ))
            for field in elimination_fields
                _disabled(parameters[field]) || throw(ArgumentError(
                    "PSCAD field $field requests conductor elimination; requested ports must be retained",
                ))
            end
        end
        _report(console, verbosity, 2, "Verified that all cable terminals are retained")
        project.save()
        phase = "load geometry"
        moved = _load_relocations(origins, _schematic_origins(project_path))
        if !isempty(moved)
            relocation_path = joinpath(output, "load-relocations.txt")
            write(relocation_path, join(moved, '\n') * "\n")
            error("PSCAD moved $(length(moved)) schematic objects on load; see $relocation_path\n" *
                join(first(moved, 20), '\n'))
        end

        _record_diagnostics(console, project, verbosity, "PSCAD diagnostics before computation")
        if mode == "transient"
            phase = "project run"
            _report(console, verbosity, 1, "Starting PSCAD project.run()")
            started = time_ns()
            # PowerImpedance.jl `src/PSCAD/remote/runner.jl:235-259`.
            project.run()
            elapsed = (time_ns() - started) * 1e-9
            write(joinpath(output, "timing.txt"), string(elapsed))
            phase = "build messages"
            _record_diagnostics(console, project, verbosity, "PSCAD diagnostics after run call")
            _build_messages(project, build_path, dirname(project_path))
            roots = [abspath(output), abspath(dirname(project_path)),
                abspath(_string(project.temp_folder))]
            phase = "transient-output completion"
            settings = readbacks["project"]
            expected_rows = round(Int, parse(Float64, string(settings["time_duration"])) /
                (parse(Float64, string(settings["sample_step"])) * 1e-6)) + 1
            files = _wait_outputs(roots, output_stem, expected_rows)
            for source in (files.inf, files.outs...)
                destination = joinpath(output, basename(source))
                abspath(source) == abspath(destination) || cp(source, destination; force = true)
            end
            _report(console, verbosity, 1, "Collected PSCAD transient outputs")
        else
            phase = "line-constants computation"
            _report(console, verbosity, 1, "Starting PSCAD line-constants computation")
            started = time_ns()
            line.compile()
            elapsed = (time_ns() - started)*1e-9
            _report(
                console,
                verbosity,
                1,
                "PSCAD compile call returned in $(round(elapsed; digits = 3)) seconds"
            )
            _record_diagnostics(console, project, verbosity, "PSCAD diagnostics after compile call")
            phase = "build messages"
            _build_messages(project, build_path, dirname(project_path))
            write(joinpath(output, "timing.txt"), string(elapsed))
            roots = [
                abspath(output),
                abspath(dirname(project_path)),
                abspath(_string(project.temp_folder))
            ]
            phase = "detailed-output completion"
            expected_rows=parse(Int, numf) + 1
            _report(
                console,
                verbosity,
                1,
                "Waiting for $expected_rows rows in each detailed PSCAD matrix file"
            )
            for (suffix, name) in OUTPUT_SUFFIXES
                source=_wait_output(roots, suffix, expected_rows)
                destination=joinpath(output, name)
                cp(source, destination; force = true)
                _data_rows(destination) == expected_rows || throw(ArgumentError(
                    "copied PSCAD output $destination is incomplete",
                ))
                _report(console, verbosity, 2, "Validated $name with $expected_rows rows")
            end
            _report(console, verbosity, 1, "Collected detailed PSCAD Z and Y outputs")
        end
        observed = pyconvert(Dict{String, String}, identify(pscad_version, app))
        observed == initial_identity ||
            error("PSCAD installation changed during computation")
        open(joinpath(output, "solver.toml"), "w") do io
            TOML.print(io, observed; sorted = true)
        end
    catch error
        project === nothing || _record_diagnostics(
            console,
            project,
            verbosity,
            "PSCAD diagnostics at failure";
            error_stream = true
        )
        message="PSCAD runner failed during $phase: $(sprint(showerror, error))"
        if project !== nothing
            try
                _build_messages(project, build_path, dirname(project_path))
            catch build_error
                diagnostic = sprint(showerror, build_error)
                occursin(diagnostic, message) || (message *= "\n$diagnostic")
            end
        end
        println(console, message)
        flush(console)
        throw(ErrorException(message))
    finally
        if project === nothing && load_attempted
            try
                loaded = app.project(project_name)
                normpath(abspath(_string(loaded.filename))) ==
                normpath(abspath(project_path)) &&
                    (project = loaded)
            catch error
                println(console, "PSCAD project lookup during cleanup failed: ", sprint(showerror, error))
            end
        end
        project === nothing || try
            project.unload()
        catch error
            println(console, "PSCAD project unload failed: ", sprint(showerror, error))
        end
        if app !== nothing && !owned && previous_silence !== nothing
            try
                app.silence = previous_silence
            catch error
                println(console, "PSCAD silence restore failed: ", sprint(showerror, error))
            end
        end
        app === nothing || !owned ||
            try
                app.quit()
            catch error
                println(console, "PSCAD application close failed: ", sprint(showerror, error))
            end
        flush(console)
        close(console)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main(ARGS)
    catch error
        println(stderr, sprint(showerror, error))
        exit(1)
    end
end
