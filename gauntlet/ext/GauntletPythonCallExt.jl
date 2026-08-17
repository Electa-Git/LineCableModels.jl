module GauntletPythonCallExt

using Gauntlet
using JSON3
using PythonCall
using TOML

import Gauntlet: automate

const G = Gauntlet

_string(value) = pyconvert(String, pybuiltins.str(value))
_bool(value) = pyconvert(Bool, value)
_int(value) = pyconvert(Int, value)

function _plain(value::Py)
    pyisnone(value) && return nothing
    pyisinstance(value, pybuiltins.bool) && return pyconvert(Bool, value)
    pyisinstance(value, pybuiltins.int) && return pyconvert(Int, value)
    pyisinstance(value, pybuiltins.float) && return pyconvert(Float64, value)
    pyisinstance(value, pybuiltins.str) && return pyconvert(String, value)
    pyisinstance(value, pybuiltins.dict) && return Dict(
        _string(key) => _plain(item) for (key, item) in PyDict{Py, Py}(value)
    )
    if pyisinstance(value, pybuiltins.list) || pyisinstance(value, pybuiltins.tuple) ||
       pyisinstance(value, pybuiltins.set) || pyisinstance(value, pybuiltins.range)
        return [_plain(item) for item in PyList{Py}(pybuiltins.list(value))]
    end
    return _string(value)
end

_values(component) = _plain(component.parameters())

function _attribute(value, name::AbstractString, default = nothing)
    pyhasattr(value, name) || return default
    return _plain(pygetattr(value, name))
end

function _components(value)
    return pyconvert(Vector{Py}, value)
end

function _lines(project)
    vcat(
        _components(project.find_all("TLine")), _components(project.find_all("Cable"))
    )
end

function _canvas(project, line)
    try
        return line.canvas()
    catch
        return project.canvas(last(split(_string(line.defn_name), ':'; limit = 2)))
    end
end

function _choices(component, name)
    raw = component.range(name)
    return sort!([_plain(item) for item in PyList{Py}(pybuiltins.list(raw))]; by = string)
end

function _choose(available, tokens...)
    for token in tokens, item in available

        occursin(lowercase(token), lowercase(string(item))) && return item
    end
    throw(ArgumentError("none of $(tokens) appears in $(available)"))
end

function _snapshot(component; ranges::Bool = false)
    parameters = _values(component)
    row = Dict{String, Any}(
        "iid" => _attribute(component, "iid"),
        "type" => _string(pybuiltins.type(component).__name__),
        "definition" => _attribute(component, "defn_name"),
        "parameters" => parameters
    )
    if ranges
        selectors = Set(("EarthForm", "EarthForm2", "EarthForm3", "GrRho",
            "Output", "Numf", "FDIS", "CPASS", "NFP", "DCenab", "DCCOR",
            "Interp1", "ElimGW", "LC", "elimp"))
        union!(selectors, filter(name -> startswith(lowercase(name), "elim"), keys(parameters)))
        values = Dict{String, Any}()
        for name in intersect(keys(parameters), selectors)
            try
                values[name] = _choices(component, name)
            catch
            end
        end
        row["ranges"] = values
    end
    return row
end

function _line_snapshot(project, line; ranges::Bool = false)
    return Dict(
        "outer" => _snapshot(line; ranges),
        "internal_canvas" => [_snapshot(component; ranges)
         for component in _components(_canvas(project, line).components())]
    )
end

function _unique_lines(project)
    rows = G.LineKey[]
    seen = Set{String}()
    for line in _lines(project)
        definition = _string(line.defn_name)
        lowercase(definition) in seen && continue
        push!(seen, lowercase(definition))
        parameters = _values(line)
        name = string(get(parameters, "Name", get(parameters, "name", last(split(definition, ':')))))
        push!(rows, G.LineKey(_int(line.iid), definition,
            _string(pybuiltins.type(line).__name__), name))
    end
    return rows
end

function _find_line(project, key::G.LineKey)
    lines = _lines(project)
    match = findfirst(line -> _int(line.iid) == key.iid, lines)
    match === nothing || return lines[match]
    match = findfirst(
        line -> lowercase(_string(line.defn_name)) ==
                lowercase(key.definition), lines)
    match === nothing && throw(KeyError(key.definition))
    return lines[match]
end

function _relevant(project, line)
    components = _components(_canvas(project, line).components())
    parameters = _values.(components)
    ground = findfirst(eachindex(components)) do index
        definition = lowercase(_string(components[index].defn_name))
        endswith(definition, ":line_ground") ||
            all(in(keys(parameters[index])), ("EarthForm", "EarthForm2", "GRRES"))
    end
    frequency = findfirst(eachindex(components)) do index
        definition = lowercase(_string(components[index].defn_name))
        endswith(definition, ":line_frephase_options") ||
            all(in(keys(parameters[index])), ("FS", "FE", "Numf", "Output"))
    end
    definitions = lowercase.(_string.(getproperty.(components, :defn_name)))
    aerial = any(name -> occursin("line_tower", name), definitions) ||
             any(parameters) do values
        haskey(values, "ElimGW") && haskey(values, "NG")
    end
    underground = any(definitions) do name
        occursin("cable_", name) || occursin("ground_plane", name)
    end || any(parameters) do values
        haskey(values, "LL") || haskey(values, "cnum") ||
        any(name -> occursin(r"^LL\d+$", name), keys(values))
    end
    return (;
        components,
        parameters,
        ground = ground === nothing ? nothing : components[ground],
        frequency = frequency === nothing ? nothing : components[frequency],
        aerial,
        underground
    )
end

function _conductor_layers(layout)
    layout isa Number && return Int(layout) ÷ 2 + 1
    return max(1, length(collect(eachmatch(r"(?:^|_)C\d+"i, string(layout)))))
end

function _eliminated(value, layers, parameters, suffix = "")
    text = lowercase(string(value))
    occursin("all", text) && occursin("concentric", text) && return layers
    occursin("outermost", text) && return min(1, layers)
    if occursin("spec", text)
        return count(1:layers) do layer
            value = get(parameters, "elim$(layer)$(suffix)", nothing)
            value === true || occursin("elim", lowercase(string(value)))
        end
    end
    return 0
end

function _requested(setting, available)
    setting === nothing && return available
    setting == "available" && return available
    requested = setting isa AbstractVector ? setting : [setting]
    lookup = Dict(lowercase(string(item)) => item for item in available)
    all(item -> haskey(lookup, lowercase(string(item))), requested) || throw(ArgumentError(
        "requested PSCAD choice is absent from $(available)",
    ))
    return [lookup[lowercase(string(item))] for item in requested]
end

function _formulation_axes(info, config)
    settings = config["axes"]["earth_formulation"]
    get(settings, "enabled", true) || return Pair{String, Vector}[]
    info.ground === nothing && return Pair{String, Vector}[]
    axes = Pair{String, Vector}[]
    info.aerial && push!(axes,
        "EarthForm2" => _requested(get(settings, "aerial", nothing),
            _choices(info.ground, "EarthForm2")))
    info.underground && push!(axes,
        "EarthForm" => _requested(get(settings, "underground", nothing),
            _choices(info.ground, "EarthForm")))
    info.aerial && info.underground &&
        push!(axes,
            "EarthForm3" => _requested(get(settings, "mixed", nothing),
                _choices(info.ground, "EarthForm3")))
    return axes
end

function _products(axes)
    isempty(axes) && return [Dict{String, Any}()]
    result = [Dict{String, Any}()]
    for (name, choices) in axes
        result = [merge(row, Dict(name => choice)) for row in result for choice in choices]
    end
    return result
end

function _has_reduction(info)
    return any(info.parameters) do values
        haskey(values, "ElimGW") && Int(get(values, "NG", 0)) > 0 ||
            haskey(values, "LC") ||
            haskey(values, "elimp") && Bool(get(values, "poins", true)) ||
            any(name -> occursin(r"^LC\d+$"i, name) || startswith(lowercase(name), "elim"), keys(values))
    end
end

function _variants(info, config)
    get(config["frequency_dependent_phase"], "enabled", true) &&
        info.frequency === nothing &&
        return G.Variant[]
    rho = config["axes"]["earth_resistivity"]
    resistivities = get(rho, "enabled", true) && info.ground !== nothing ?
                    Union{Nothing, Float64}[Float64(value) for value in rho["values_ohm_m"]] :
                    Union{Nothing, Float64}[nothing]
    formulations = _products(_formulation_axes(info, config))
    reduction = config["axes"]["kron_reduction"]
    reductions = get(reduction, "enabled", true) && _has_reduction(info) ?
                 String.(get(reduction, "states", ["retained", "reduced"])) :
                 ["not_applicable"]
    return [G.Variant(rho, formulation, state)
            for rho in resistivities, formulation in formulations, state in reductions][:]
end

function _set_parameters(component, values)
    isempty(values) || component.parameters(;
        (Symbol(name) => value for (name, value) in values)...)
end

function _frequency!(component, config)
    settings = config["frequency_dependent_phase"]
    parameters = _values(component)
    passivity = get(settings, "passivity", Dict{String, Any}())
    requested = Dict{String, Any}(
        "Output" => get(settings, "detailed_output", true) ? "YES" : "NO",
        "FS" => Float64(settings["start_hz"]),
        "FE" => Float64(settings["end_hz"]),
        "Numf" => Int(settings["samples"]),
        "FDIS" => get(settings, "distribution", "LOG_LINEAR"),
        "Interp1" => get(settings, "interpolate_travel_time", true) ? "YES" : "NO",
        "DCenab" => get(settings, "dc_correction", true) ? "ENABLED" : "DISABLED",
        "DCCOR" => get(settings, "dc_correction_method", "FUNCTIONAL_FORM"),
        "CPASS" => get(passivity, "mode", "DISABLE"),
        "NFP" => Int(get(passivity, "samples", 1000)),
        "FSP" => Float64(get(passivity, "start_hz", settings["start_hz"])),
        "FEP" => Float64(get(passivity, "end_hz", settings["end_hz"]))
    )
    mutation = Dict(name => value
    for (name, value) in requested if haskey(parameters, name))
    _set_parameters(component, mutation)
    return mutation
end

function _reduction!(line, info, state)
    state == "not_applicable" && return Dict{String, Any}[]
    state in ("retained", "reduced") || throw(ArgumentError("unknown reduction $state"))
    changed = Dict{String, Any}[]
    dimension_delta = 0
    for (component, parameters) in zip(info.components, info.parameters)
        mutation = Dict{String, Any}()
        if haskey(parameters, "ElimGW") && Int(get(parameters, "NG", 0)) > 0
            old = occursin("enable", lowercase(string(parameters["ElimGW"]))) ?
                  Int(parameters["NG"]) : 0
            new = state == "reduced" ? Int(parameters["NG"]) : 0
            dimension_delta += old - new
            mutation["ElimGW"] = _choose(_choices(component, "ElimGW"),
                state == "retained" ? "DISABLED" : "ENABLED")
        end
        if haskey(parameters, "LC")
            layers = _conductor_layers(get(parameters, "LL", 1)) - 1
            old = _eliminated(parameters["LC"], layers, parameters)
            new = state == "reduced" ? layers : 0
            dimension_delta += old - new
            tokens = state == "retained" ? ("NONE", "RETAIN") :
                     ("ALL_CONCENTRIC", "OUTERMOST", "ELIM")
            mutation["LC"] = _choose(_choices(component, "LC"), tokens...)
        elseif haskey(parameters, "elimp")
            if Bool(get(parameters, "poins", true))
                old = Bool(parameters["elimp"]) ? 1 : 0
                new = state == "reduced" ? 1 : 0
                dimension_delta += old - new
                mutation["elimp"] = state == "reduced"
            end
            for index in 1:Int(get(parameters, "cnum", 0))
                lc = "LC$index"
                haskey(parameters, lc) || continue
                layers = _conductor_layers(get(parameters, "LL$index", 1)) - 1
                old = _eliminated(parameters[lc], layers, parameters, string(index))
                new = state == "reduced" ? layers : 0
                dimension_delta += old - new
                tokens = state == "retained" ? ("NONE", "RETAIN") :
                         ("ALL_CONCENTRIC", "OUTERMOST", "ELIM")
                mutation[lc] = _choose(_choices(component, lc), tokens...)
            end
        elseif !haskey(parameters, "ElimGW")
            for name in keys(parameters)
                startswith(lowercase(name), "elim") || continue
                choices = try
                    _choices(component, name)
                catch
                    continue
                end
                tokens = state == "retained" ? ("RETAIN", "DISABLE") :
                         ("ELIM", "ENABLE")
                mutation[name] = _choose(choices, tokens...)
            end
        end
        _set_parameters(component, mutation)
        isempty(mutation) || push!(changed,
            Dict(
                "component" => _string(component.defn_name),
                "iid" => _int(component.iid),
                "parameters" => mutation
            ))
    end
    if dimension_delta != 0
        outer = _values(line)
        old = Int(get(outer, "Dim", 0))
        if old > 0
            new = old + dimension_delta
            new > 0 ||
                throw(DomainError(new, "reduction produced an invalid line dimension"))
            _set_parameters(line, Dict("Dim" => new))
            push!(changed,
                Dict("component" => _string(line.defn_name),
                    "iid" => _int(line.iid), "parameters" => Dict("Dim" => new)))
        end
    end
    return changed
end

function _apply!(project, line, variant, config)
    info = _relevant(project, line)
    frequency = get(config["frequency_dependent_phase"], "enabled", true) ?
                _frequency!(info.frequency, config) : Dict{String, Any}()
    ground = copy(variant.formulations)
    if info.ground !== nothing && variant.resistivity_ohm_m !== nothing
        parameters = _values(info.ground)
        haskey(parameters, "GrRho") && (ground["GrRho"] = _choose(
            _choices(info.ground, "GrRho"), "CONSTANT_RESISTIVITY"))
        ground["GRRES"] = variant.resistivity_ohm_m
    end
    info.ground === nothing || _set_parameters(info.ground, ground)
    return Dict("frequency" => frequency, "ground" => ground,
        "reduction" => _reduction!(line, info, variant.reduction))
end

function _config(path)
    defaults = TOML.parsefile(joinpath(G._HARVEST_ROOT, "benchmark.toml"))
    path === nothing && return defaults
    function combine(base, override)
        result = deepcopy(base)
        for (key, value) in override
            result[key] = value isa AbstractDict &&
                          get(result, key, nothing) isa AbstractDict ?
                          combine(result[key], value) : deepcopy(value)
        end
        return result
    end
    return combine(defaults, TOML.parsefile(path))
end

function _session(config)
    settings = config["pscad"]
    module_value = pyimport("mhi.pscad")
    app = module_value.launch(;
        version = String(settings["version"]),
        x64 = Bool(get(settings, "x64", true)),
        minimize = Bool(get(settings, "minimize", true)),
        splash = false,
        timeout = Int(get(settings, "launch_timeout_seconds", 60))
    )
    if !_bool(app.licensed())
        Bool(get(settings, "acquire_certificate", true)) || error(
            "PSCAD has no active certificate and acquisition is disabled",
        )
        certificates = collect(values(PyDict{Py, Py}(app.get_available_certificates(;
            refresh = true))))
        available = filter(value -> _bool(value.available), certificates)
        isempty(available) && error("no PSCAD certificate is available")
        result = _int(app.get_certificate(first(available)))
        result == 0 && _bool(app.licensed()) || error(
            "PSCAD certificate acquisition failed with code $result",
        )
    end
    return app, _string(app.get_current_certificate())
end

function _load(app, source)
    app.load(abspath(source))
    return app.project(G._project_identity(source).name)
end

function _inspect(source, config)
    app, certificate = _session(config)
    project = nothing
    try
        project = _load(app, source)
        definitions = Dict{String, Any}[]
        for key in _unique_lines(project)
            line = _find_line(project, key)
            info = _relevant(project, line)
            variants = _variants(info, config)
            push!(definitions,
                Dict(
                    "key" => Dict("iid" => key.iid, "definition" => key.definition,
                        "component_type" => key.component_type, "name" => key.name),
                    "geometry_class" => Dict("aerial" => info.aerial,
                        "underground" => info.underground,
                        "mixed" => info.aerial && info.underground),
                    "variant_count" => length(variants),
                    "variants" => G._variant_dict.(variants),
                    "model" => _line_snapshot(project, line; ranges = true)
                ))
        end
        return Dict(
            "source" => abspath(source), "source_sha256" => G._sha256(source),
            "pscad_version" => config["pscad"]["version"],
            "certificate" => certificate, "definitions" => definitions,
            "total_variants" => sum(row["variant_count"] for row in definitions)
        )
    finally
        project === nothing || try
            project.unload()
        catch
        end
        try
            app.quit()
        catch
        end
    end
end

function _copy_project(source, destination)
    ispath(destination) && rm(destination; recursive = true, force = true)
    cp(dirname(source), destination; force = true)
    for (directory, directories, names) in walkdir(destination; topdown = false)
        for name in directories
            occursin(r"\.(?:gf|if|ef|gfortran|intel|clang|vc)\w*$"i, name) &&
                rm(joinpath(directory, name); recursive = true, force = true)
        end
        for name in names
            extension = lowercase(splitext(name)[2])
            if extension == ".pscx" && lowercase(name) != lowercase(basename(source)) ||
               extension == ".pswx"
                rm(joinpath(directory, name); force = true)
            end
        end
    end
    return joinpath(destination, basename(source))
end

function _run(source, config, options)
    output = abspath(String(get(options, "output", config["output"]["root"])))
    mkpath(output)
    app, certificate = _session(config)
    summary = Dict{String, Any}(
        "source" => abspath(source), "source_sha256" => G._sha256(source),
        "pscad_version" => String(config["pscad"]["version"]),
        "certificate" => certificate, "records" => Dict{String, Any}[],
        "excluded_definitions" => Dict{String, Any}[]
    )
    group = joinpath(output,
        G._slug(splitext(basename(source))[1], 38) * "_" * summary["source_sha256"][1:8])
    donor = nothing
    try
        donor = _load(app, source)
        keys = _unique_lines(donor)
        selector = get(options, "line", nothing)
        selector === nothing || filter!(
            key -> begin
                needle = lowercase(String(selector))
                occursin(needle, lowercase(key.definition)) ||
                    occursin(needle, lowercase(key.name)) ||
                    needle == string(key.iid)
            end,
            keys)
        isempty(keys) && throw(ArgumentError("no line/cable definition matched"))
        plan = Tuple{G.LineKey, G.Variant}[]
        models = Dict{String, Any}()
        for key in keys
            line = _find_line(donor, key)
            models[key.definition] = _line_snapshot(donor, line; ranges = true)
            variants = _variants(_relevant(donor, line), config)
            isempty(variants) && push!(summary["excluded_definitions"],
                Dict(
                    "definition" => key.definition,
                    "reason" => "not a frequency-dependent phase model"
                ))
            append!(plan, ((key, variant) for variant in variants))
        end
        donor.unload()
        donor = nothing
        haskey(options, "limit") &&
            resize!(plan, min(length(plan), parse(Int, options["limit"])))
        summary["planned_variants"] = length(plan)
        if haskey(options, "dry-run")
            summary["plan"] = [Dict("definition" => key.definition,
                                   "variant" => G._variant_dict(variant))
                               for (key, variant) in plan]
            return summary
        end
        for (position, (key, variant)) in enumerate(plan)
            specification = Dict(
                "source_sha256" => summary["source_sha256"],
                "line_definition" => key.definition,
                "variant" => G._variant_dict(variant),
                "frequency_dependent_phase" => config["frequency_dependent_phase"],
                "pscad_version" => summary["pscad_version"]
            )
            spec = G._sha256_json(specification)
            record = joinpath(group, G._slug(last(split(key.definition, ':')), 30),
                "v_$(spec[1:16])")
            manifest_path = joinpath(record, "manifest.json")
            if !haskey(options, "no-resume") && isfile(manifest_path)
                existing = JSON3.read(read(manifest_path, String), Dict{String, Any})
                if get(existing, "specification_sha256", "") == spec &&
                   get(existing, "status", "") in ("success", "solver_rejected")
                    push!(summary["records"],
                        Dict("path" => record,
                            "status" => "resumed_$(existing["status"])"))
                    continue
                end
            end
            @info "Running PSCAD reference" case=position total=length(plan) definition=key.definition
            mkpath(record)
            variant_source = _copy_project(source, joinpath(record, "project"))
            manifest = Dict{String, Any}(
                "schema_version" => 1, "status" => "running",
                "source" => abspath(source), "source_sha256" => summary["source_sha256"],
                "source_label" => "PSCAD reference",
                "line" => Dict("iid" => key.iid, "definition" => key.definition,
                    "component_type" => key.component_type, "name" => key.name),
                "variant" => G._variant_dict(variant),
                "specification_sha256" => spec,
                "pscad_version" => summary["pscad_version"],
                "certificate" => certificate, "donor_model" => models[key.definition],
                "mutations" => Dict{String, Any}(), "artifacts" => Dict{String, Any}[]
            )
            G._write_pscad_json(manifest_path, manifest)
            project = nothing
            started = time()
            try
                project = _load(app, variant_source)
                line = _find_line(project, key)
                manifest["mutations"] = _apply!(project, line, variant, config)
                manifest["declared_model"] = _line_snapshot(project, line; ranges = true)
                project.save()
                temp = abspath(_string(project.temp_folder))
                process_root = pwd()
                project_root = dirname(variant_source)
                before = G._snapshot(temp)
                process_before = G._snapshot(process_root; recursive = false)
                project_before = G._snapshot(project_root; recursive = false)
                line.compile()
                after = G._wait_for_artifacts(temp, before;
                    timeout = Float64(get(config["pscad"], "artifact_timeout_seconds", 300)),
                    quiet = Float64(get(config["pscad"], "artifact_quiet_seconds", 2)),
                    process_root, process_before, project_root, project_before)
                artifacts = joinpath(record, "artifacts")
                append!(manifest["artifacts"],
                    G._copy_changed(temp,
                        G._changed(before, after.temp), artifacts, :temp))
                append!(manifest["artifacts"],
                    G._copy_changed(process_root,
                        G._changed(process_before, after.process), artifacts, :process_cwd))
                append!(manifest["artifacts"],
                    G._copy_changed(project_root,
                        G._changed(project_before, after.project), artifacts, :project_cwd))
                inputs = [joinpath(record, row["path"])
                          for row in manifest["artifacts"]
                          if row["kind"] == "normalized_lcp_input"]
                constants = count(row -> row["kind"] == "emtdc_constants", manifest["artifacts"])
                results = count(
                    row -> row["kind"] in
                           ("human_readable_matrices", "detailed_frequency_output"),
                    manifest["artifacts"])
                isempty(inputs) ||
                    (manifest["case_id"] = G._case_id(summary["pscad_version"], inputs))
                rejection = G._solver_rejection(record, manifest["artifacts"])
                if !isempty(inputs) && rejection !== nothing
                    manifest["status"] = "solver_rejected"
                    manifest["solver_diagnostic"] = rejection
                elseif !isempty(inputs) && constants > 0 && results > 0
                    manifest["status"] = "success"
                else
                    error("incomplete LCP artifact set: inputs=$(length(inputs)), " *
                          "constants=$constants, results=$results")
                end
            catch error
                error isa InterruptException && rethrow()
                manifest["status"] = "failed"
                manifest["error"] = Dict("type" => string(typeof(error)),
                    "message" => sprint(showerror, error))
            finally
                manifest["elapsed_seconds"] = time() - started
                G._write_pscad_json(manifest_path, manifest)
                project === nothing || try
                    project.unload()
                catch
                end
            end
            push!(summary["records"],
                Dict("path" => record,
                    "status" => manifest["status"], "case_id" =>
                        get(manifest, "case_id", nothing)))
        end
        return summary
    finally
        donor === nothing || try
            donor.unload()
        catch
        end
        summary["successful"] = count(
            row -> row["status"] in
                   ("success", "resumed_success"), summary["records"])
        summary["failed"] = count(row -> row["status"] == "failed", summary["records"])
        G._write_pscad_json(joinpath(group, "run_summary.json"), summary)
        try
            app.quit()
        catch
        end
    end
end

function automate(::Gauntlet.PSCAD, ::Val{:inspect}, options)
    return _inspect(String(options["source"]), _config(get(options, "config", nothing)))
end

function automate(::Gauntlet.PSCAD, ::Val{:case}, options)
    return _run(String(options["source"]), _config(get(options, "config", nothing)), options)
end

function automate(::Gauntlet.PSCAD, ::Val{:batch}, options)
    root = abspath(String(options["source"]))
    projects = String[]
    for (directory, _, names) in walkdir(root), name in names

        lowercase(splitext(name)[2]) == ".pscx" &&
            push!(projects, joinpath(directory, name))
    end
    return Dict("root" => root,
        "runs" => [_run(source, _config(get(options, "config", nothing)), options)
         for source in sort!(projects)])
end

end
