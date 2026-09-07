const REPOSITORY_ROOT = pkgdir(LineCableModels)

const FLATTEN_IMPLEMENTATION_PATHS = (
    "src/datamodel/flatten.jl",
    "src/datamodel/baseparams/geometry.jl",
    "src/datamodel/baseparams/resistance.jl",
    "src/datamodel/baseparams/inductance.jl",
    "src/datamodel/baseparams/dielectrics.jl",
    "src/datamodel/placement/paths.jl",
    "src/datamodel/geometry/primitives.jl",
    "src/datamodel/geometry/sector.jl",
    "src/datamodel/geometry/ellipse.jl",
    "src/datamodel/design/assembly.jl",
    "src/datamodel/placement/bounded.jl",
    "src/materials/material.jl"
)

const COAXIAL_IMPLEMENTATION_PATHS = (
    FLATTEN_IMPLEMENTATION_PATHS...,
    "src/engine/formulations.jl",
    "src/engine/blueprint.jl",
    "src/engine/input.jl",
    "src/engine/impedance.jl",
    "src/engine/admittance.jl",
    "src/engine/earthreturn.jl",
    "src/engine/earthkernels.jl",
    "src/engine/matrixops.jl",
    "src/engine/reduction.jl",
    "src/engine/lineparameters.jl",
    "src/engine/options.jl",
    "src/engine/problems.jl"
)

include("fingerprints.jl")

_selection_value(value::Union{Nothing, Bool, Number, AbstractString, Symbol}) = value
_selection_value(value::Type) = string(value)
_selection_value(value::NamedTuple) = map(_selection_value, value)
_selection_value(value::Tuple) = map(_selection_value, value)
_selection_value(value::AbstractVector) = _selection_value.(value)
function _selection_value(value)
    ismutabletype(typeof(value)) && throw(ArgumentError(
        "formula provenance cannot fingerprint mutable $(typeof(value)); provide immutable numerical route inputs"))
    names = fieldnames(typeof(value))
    return (type=string(typeof(value)), fields=NamedTuple{names}(
        map(name -> _selection_value(getfield(value, name)), names)))
end


"Return a digest of the materialised numerical declarations in one problem."
function numerical_input_sha256(problem::LineCableModels.Engine.LineParametersProblem)
    # JSON transport retains the declaration, not the resolved physical tree.
    # Both are inputs to reuse: a placement change must not resurrect results
    # computed from the same declaration by an older resolver.
    resolved = map(problem.system.designs) do design
        _selection_value((geometry=design.geometry,
            terminal_order=design.terminal_order, terminal_map=design.terminal_map))
    end
    return semantic_sha256((schema_version=2,
        declaration=LineCableModels.ImportExport.serialize_value(problem), resolved))
end

function repository_provenance()
    commit = readchomp(`git -C $REPOSITORY_ROOT rev-parse HEAD`)
    dirty = !isempty(readchomp(`git -C $REPOSITORY_ROOT status --porcelain`))
    return (; commit, dirty)
end

function git_blob_record(relative::AbstractString)
    path = joinpath(REPOSITORY_ROOT, relative)
    isfile(path) || throw(ArgumentError("implementation source is missing: $relative"))
    blob = readchomp(`git -C $REPOSITORY_ROOT hash-object $relative`)
    return (path = String(relative), blob)
end

function _formula_paths(family::AbstractString, formula)
    identifier = LineCableModels.formula_id(formula)
    relative = joinpath(
        "src", "engine", family, "formulas",
        lowercase(string(identifier)) * ".jl"
    )
    isfile(joinpath(REPOSITORY_ROOT, relative)) || return String[]
    directory = dirname(joinpath(REPOSITORY_ROOT, relative))
    candidates = sort!(filter(endswith(".jl"), readdir(directory)))
    selected = Set((basename(relative),))
    pending = [basename(relative)]
    while !isempty(pending)
        source = lowercase(read(joinpath(directory, popfirst!(pending)), String))
        for candidate in candidates
            candidate in selected && continue
            stem = lowercase(splitext(candidate)[1])
            occursin(stem, source) || continue
            push!(selected, candidate)
            push!(pending, candidate)
        end
    end
    return [joinpath("src", "engine", family, "formulas", file)
            for file in sort!(collect(selected))]
end

function _ehem_path(sequence)
    sequence === nothing && return nothing
    rule = LineCableModels.EarthProps.EHEM.rule(sequence)
    rule isa LineCableModels.EarthProps.EHEM.Layer &&
        return "src/earthprops/ehem/layer.jl"
    identifier = LineCableModels.formula_id(rule)
    return joinpath(
        "src", "earthprops", "ehem", "formulas",
        lowercase(string(identifier)) * ".jl"
    )
end

function _fd_path(formula)
    formula === nothing && return nothing
    identifier = LineCableModels.formula_id(formula)
    return joinpath(
        "src", "earthprops", "fd", "formulas",
        lowercase(string(identifier)) * ".jl"
    )
end

function _selection_record(value)
    value === nothing && return nothing
    identifier = applicable(LineCableModels.formula_id, value) ?
                 LineCableModels.formula_id(value) : string(typeof(value))
    assumptions = hasproperty(value, :assumptions) ?
                  _selection_value(value.assumptions) : (;)
    routes = if hasproperty(value, :routes)
        map(_selection_value, value.routes)
    elseif hasproperty(value, :route)
        _selection_value(value.route)
    else
        nothing
    end
    return (; identifier, assumptions, routes)
end

function formulation_record(formulation::LineCableModels.Engine.LineParametersFormulation)
    methods = formulation.methods
    equivalent = methods.equivalent_earth
    equivalent_record = if equivalent === nothing
        nothing
    else
        rule = LineCableModels.EarthProps.EHEM.rule(equivalent)
        rule isa LineCableModels.EarthProps.EHEM.Layer ?
        (order = nameof(typeof(equivalent)), rule = :Layer, layer = rule.layer) :
        (order = nameof(typeof(equivalent)), rule = _selection_record(rule))
    end
    return (
        schema_version = 1,
        backend = :coaxial,
        internal_impedance = _selection_record(methods.internal_impedance),
        insulation_impedance = _selection_record(methods.insulation_impedance),
        earth_impedance = _selection_record(methods.earth_impedance),
        insulation_admittance = _selection_record(methods.insulation_admittance),
        semicon_admittance = _selection_record(methods.semicon_admittance),
        earth_admittance = _selection_record(methods.earth_admittance),
        earth_properties = _selection_record(methods.earth_properties),
        equivalent_earth = equivalent_record,
        pipe_impedance = _selection_record(methods.pipe_impedance),
        options = formulation.options
    )
end

"Return the selected source blobs and declarative fingerprint of one formulation."
function implementation_record(
        formulation::LineCableModels.Engine.LineParametersFormulation;
        external_sources = ()
)
    methods = formulation.methods
    paths = String[COAXIAL_IMPLEMENTATION_PATHS...]
    append!(paths, (
        "src/engine/internalimpedance/interface.jl",
        "src/engine/insulationimpedance/interface.jl",
        "src/engine/insulationadmittance/interface.jl",
        "src/engine/semiconadmittance/interface.jl",
        "src/engine/earthimpedance/interface.jl",
        "src/engine/earthimpedance/homogeneous.jl",
        "src/engine/earthadmittance/interface.jl",
        "src/engine/earthadmittance/homogeneous.jl",
        "src/engine/pipeimpedance/interface.jl",
        "src/earthprops/ehem/interface.jl",
        "src/earthprops/fd/interface.jl"
    ))
    for (family, formula) in (
        ("internalimpedance", methods.internal_impedance),
        ("insulationimpedance", methods.insulation_impedance),
        ("insulationadmittance", methods.insulation_admittance),
        ("semiconadmittance", methods.semicon_admittance),
        ("earthimpedance", methods.earth_impedance),
        ("earthadmittance", methods.earth_admittance),
        ("pipeimpedance", methods.pipe_impedance)
    )
        append!(paths, _formula_paths(family, formula))
    end
    for path in (
        _fd_path(methods.earth_properties),
        _ehem_path(methods.equivalent_earth)
    )
        path === nothing || push!(paths, path)
    end
    append!(paths, String.(external_sources))
    unique!(sort!(paths))
    selection = formulation_record(formulation)
    return (
        selection,
        selection_sha256 = semantic_sha256(selection),
        blobs = git_blob_record.(paths)
    )
end
