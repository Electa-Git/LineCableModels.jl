using PythonCall

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
    return uppercase(strip(string(value))) in
           ("0", "NO", "FALSE", "DISABLE", "DISABLED", "RETAIN", "NONE")
end

function _verify_retained_ports(components)
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
                "PSCAD field $field requests conductor elimination; gauntlet ports must be retained",
            ))
        end
    end
    return nothing
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

function main(arguments)
    length(arguments) == 13 || throw(ArgumentError(
        "runner expects project, output, project name, output stem, formulation, " *
        "earth field, earth value, earth readback, FS, FE, Numf, PSCAD version, " *
        "and verbosity",
    ))
    project_path, output, project_name, output_stem, formulation, earth_field,
    earth_value, earth_readback, fs, fe, numf, pscad_version,
    verbosity_text = arguments
    verbosity = parse(Int, verbosity_text)
    verbosity in 0:2 || throw(ArgumentError("verbosity must be 0, 1, or 2"))
    pscad_version == "5.1.0" || throw(ArgumentError("PSCAD 5.1.0 is required"))
    earth_field in ("EarthForm", "EarthForm2", "EarthForm3") || throw(ArgumentError(
        "unsupported PSCAD earth formulation field $earth_field",
    ))
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_]{0,19}$", output_stem) ||
        throw(ArgumentError("invalid PSCAD output stem $output_stem"))
    mkpath(output)
    app = nothing
    project = nothing
    phase = "initialization"
    console = open(joinpath(output, "pscad-console.txt"), "w")
    try
        metadata = pyimport("importlib.metadata")
        observed = _string(metadata.version("mhi.pscad"))
        observed == MHI_PSCAD_VERSION || throw(ArgumentError(
            "mhi.pscad $MHI_PSCAD_VERSION is required; found $observed",
        ))
        pscad = pyimport("mhi.pscad")
        phase = "PSCAD launch"
        _report(console, verbosity, 1, "Launching PSCAD 5.1.0")
        app = cd(abspath(output)) do
            pscad.launch(;
                version = pscad_version,
                x64 = true,
                minimize = true,
                splash = false,
                silence = true,
                timeout = 60
            )
        end
        pyconvert(Bool, app.licensed()) || error("PSCAD refused the configured license")
        phase = "project load"
        _report(console, verbosity, 1, "Loading generated project $project_name")
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
        _set!(line, "Name", output_stem)
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
        frequency = _single_component(components, "master:line_frephase_options")
        ground = _single_component(components, "master:line_ground")
        _set!(frequency, "FS", parse(Float64, fs))
        _set!(frequency, "FE", parse(Float64, fe))
        _set!(frequency, "Numf", parse(Int, numf))
        _set!(frequency, "FDIS", 3; readback = "LOG_LINEAR")
        _set!(frequency, "Output", 1; readback = "YES")
        _set!(
            ground,
            earth_field,
            parse(Int, earth_value);
            readback = earth_readback
        )
        _report(
            console,
            verbosity,
            2,
            "Applied frequency range $fs Hz to $fe Hz and formulation $formulation"
        )
        _verify_retained_ports(components)
        _report(console, verbosity, 2, "Verified that all cable terminals are retained")
        project.save()
        _record_diagnostics(console, project, verbosity, "PSCAD diagnostics before calculation")
        phase = "line-constants calculation"
        _report(console, verbosity, 1, "Starting PSCAD line-constants calculation")
        started = time()
        line.compile()
        elapsed = time() - started
        _report(
            console,
            verbosity,
            1,
            "PSCAD line-constants calculation finished in $(round(elapsed; digits = 3)) seconds"
        )
        _record_diagnostics(console, project, verbosity, "PSCAD diagnostics after calculation")
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
    catch error
        project === nothing || _record_diagnostics(
            console,
            project,
            verbosity,
            "PSCAD diagnostics at failure";
            error_stream = true
        )
        message="PSCAD runner failed during $phase: $(sprint(showerror, error))"
        println(console, message)
        flush(console)
        throw(ErrorException(message))
    finally
        project === nothing || try
            project.unload()
        catch error
            println(console, "PSCAD project unload failed: ", sprint(showerror, error))
        end
        app === nothing || try
            app.quit()
        catch error
            println(console, "PSCAD application close failed: ", sprint(showerror, error))
        end
        flush(console)
        close(console)
    end
    return nothing
end

try
    main(ARGS)
catch error
    println(stderr, sprint(showerror, error))
    exit(1)
end
