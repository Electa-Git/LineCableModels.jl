using PythonCall

const MHI_PSCAD_VERSION = "3.1.2"
const OUTPUT_SUFFIXES = (
    "_zm.out" => "result_zm.out",
    "_zp.out" => "result_zp.out",
    "_ym.out" => "result_ym.out",
    "_yp.out" => "result_yp.out"
)
const FORMULATIONS = Dict(
    "overhead_deri_carson_lossless" =>
        (field = "EarthForm2", value = 0, readback = "DERISEMLYEN"),
    "underground_wedepohl_pollaczek_lossless" =>
        (field = "EarthForm", value = 0, readback = "WEDEPOHL")
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

function _find_output(roots, suffix)
    matches = String[]
    for root in unique(filter(isdir, roots))
        for (directory, _, files) in walkdir(root), file in files

            endswith(lowercase(file), suffix) &&
                push!(matches, realpath(joinpath(directory, file)))
        end
    end
    unique!(matches)
    length(matches) == 1 || throw(ArgumentError(
        "PSCAD emitted $(length(matches)) files ending in $suffix",
    ))
    return only(matches)
end

function main(arguments)
    length(arguments) == 9 || throw(ArgumentError(
        "runner expects project, output, project name, formulation, FS, FE, Numf, PSCAD version, and verbosity",
    ))
    project_path, output, project_name, formulation, fs, fe,
    numf, pscad_version,
    verbosity_text = arguments
    verbosity = parse(Int, verbosity_text)
    verbosity in 0:2 || throw(ArgumentError("verbosity must be 0, 1, or 2"))
    pscad_version == "5.1.0" || throw(ArgumentError("PSCAD 5.1.0 is required"))
    haskey(FORMULATIONS, formulation) || throw(ArgumentError(
        "unsupported PSCAD formulation $formulation",
    ))
    mkpath(output)
    app = nothing
    project = nothing
    console = open(joinpath(output, "pscad-console.txt"), "w")
    try
        metadata = pyimport("importlib.metadata")
        observed = _string(metadata.version("mhi.pscad"))
        observed == MHI_PSCAD_VERSION || throw(ArgumentError(
            "mhi.pscad $MHI_PSCAD_VERSION is required; found $observed",
        ))
        pscad = pyimport("mhi.pscad")
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
        _report(
            console,
            verbosity,
            2,
            "Found line-data row $(_string(line.defn_name))"
        )
        canvas = try
            project.canvas(last(split(_string(line.defn_name), ':'; limit = 2)))
        catch
            line.canvas()
        end
        components = _components(canvas.components())
        frequency = _single_component(components, "master:line_frephase_options")
        ground = _single_component(components, "master:line_ground")
        _set!(frequency, "FS", parse(Float64, fs))
        _set!(frequency, "FE", parse(Float64, fe))
        _set!(frequency, "Numf", parse(Int, numf))
        _set!(frequency, "FDIS", 3; readback = "LOG_LINEAR")
        _set!(frequency, "Output", 1; readback = "YES")
        specification = FORMULATIONS[formulation]
        _set!(
            ground,
            specification.field,
            specification.value;
            readback = specification.readback
        )
        _report(
            console,
            verbosity,
            2,
            "Applied frequency range $fs Hz to $fe Hz and formulation $formulation"
        )
        project.save()
        if verbosity >= 2
            load_messages = _messages(project)
            println(console, load_messages)
            flush(console)
            println(stdout, load_messages)
            flush(stdout)
        end
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
        messages = _messages(project)
        println(console, messages)
        flush(console)
        if verbosity >= 2
            println(stdout, messages)
            flush(stdout)
        end
        try
            project_output = _string(project.output())
            println(console, project_output)
            flush(console)
            if verbosity >= 2
                println(stdout, project_output)
                flush(stdout)
            end
        catch error
            println(console, "Unable to read PSCAD project output: ", sprint(showerror, error))
        end
        write(joinpath(output, "timing.txt"), string(elapsed))
        roots = [
            abspath(output),
            abspath(dirname(project_path)),
            abspath(_string(project.temp_folder))
        ]
        for (suffix, name) in OUTPUT_SUFFIXES
            cp(_find_output(roots, suffix), joinpath(output, name); force = true)
        end
        _report(console, verbosity, 1, "Collected detailed PSCAD Z and Y outputs")
    catch error
        println(console, "PSCAD runner failed: ", sprint(showerror, error, catch_backtrace()))
        flush(console)
        rethrow()
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

main(ARGS)
