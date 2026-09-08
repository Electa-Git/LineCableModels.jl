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
    "src/materials/material.jl",
    "src/materials/radialdielectric.jl"
)

const COAXIAL_IMPLEMENTATION_PATHS = (
    FLATTEN_IMPLEMENTATION_PATHS...,
    "src/materials/temperaturedependent/TemperatureDependent.jl",
    "src/materials/temperaturedependent/interface.jl",
    "src/materials/temperaturedependent/formulas/default.jl",
    "src/engine/formulations.jl",
    "src/engine/blueprint.jl",
    "src/engine/input.jl",
    "src/engine/impedance.jl",
    "src/engine/admittance.jl",
    "src/engine/earthreturn.jl",
    "src/engine/earthkernels.jl",
    "src/engine/integration.jl",
    "src/formulas.jl",
    "src/grammar/formulas.jl",
    "src/engine/matrixops.jl",
    "src/engine/reduction.jl",
    "src/engine/lineparameters.jl",
    "src/engine/options.jl",
    "src/engine/problems.jl"
)

include("fingerprints.jl")

_selection_value(value::Union{Nothing, Bool, Number, AbstractString, Symbol}) = value
function _selection_value(value::Type)
    sprint(show, value; context = (:module=>nothing, :compact=>false))
end
_selection_value(value::NamedTuple) = map(_selection_value, value)
_selection_value(value::Tuple) = map(_selection_value, value)
_selection_value(value::AbstractVector) = _selection_value.(value)
function _selection_value(value)
    ismutabletype(typeof(value)) && throw(ArgumentError(
        "formula provenance cannot fingerprint mutable $(typeof(value)); provide immutable numerical route inputs"))
    names = fieldnames(typeof(value))
    return (type = _selection_value(typeof(value)),
        fields = NamedTuple{names}(
            map(name -> _selection_value(getfield(value, name)), names)))
end

"Return a digest of the materialised numerical declarations in one problem."
function numerical_input_sha256(problem::LineCableModels.Engine.LineParametersProblem)
    # JSON transport retains the declaration, not the resolved physical tree.
    # Both are inputs to reuse: a placement change must not resurrect results
    # computed from the same declaration by an older resolver.
    resolved = map(problem.system.designs) do design
        _selection_value((geometry = design.geometry,
            terminal_order = design.terminal_order, terminal_map = design.terminal_map))
    end
    return semantic_sha256((schema_version = 2,
        declaration = LineCableModels.ImportExport.serialize_value(problem), resolved))
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
    return [relative]
end

function _ehem_path(sequence)
    sequence === nothing && return nothing
    rule = LineCableModels.Earth.EquivalentHomogeneous.rule(sequence)
    identifier = LineCableModels.formula_id(rule)
    return joinpath(
        "src", "earth", "equivalenthomogeneous", "formulas",
        lowercase(string(identifier)) * ".jl"
    )
end

function _fd_path(formula)
    formula === nothing && return nothing
    identifier = LineCableModels.formula_id(formula)
    return joinpath(
        "src", "earth", "frequencydependent", "formulas",
        lowercase(string(identifier)) * ".jl"
    )
end

function _selection_record(value)
    value === nothing && return nothing
    identifier = applicable(LineCableModels.formula_id, value) ?
                 LineCableModels.formula_id(value) : string(typeof(value))
    assumptions = hasproperty(value, :assumptions) ?
                  _selection_value(value.assumptions) : (;)
    binding = hasproperty(value, :binding) ? _selection_value(value.binding) : nothing
    parameters = hasproperty(value, :parameters) ? _selection_value(value.parameters) : (;)
    hooks = hasproperty(value, :hooks) ? _selection_value(value.hooks) : (;)
    options = hasproperty(value, :options) ? _selection_value(value.options) : (;)
    equivalent_earth = if hasproperty(value, :equivalent_earth) &&
                          value.equivalent_earth !== nothing
        sequence = value.equivalent_earth
        (order = nameof(typeof(sequence)), rule = _selection_record(sequence.rule))
    else
        nothing
    end
    return (;
        identifier, assumptions, parameters, hooks, binding, options, equivalent_earth)
end

function formulation_record(formulation::LineCableModels.Engine.LineParametersFormulation)
    methods = formulation.methods
    return (
        schema_version = 2,
        backend = :coaxial,
        internal_impedance = _selection_record(methods.internal_impedance),
        insulation_impedance = _selection_record(methods.insulation_impedance),
        earth_impedance = _selection_record(methods.earth_impedance),
        insulation_admittance = _selection_record(methods.insulation_admittance),
        semicon_admittance = _selection_record(methods.semicon_admittance),
        earth_admittance = _selection_record(methods.earth_admittance),
        earth_properties = _selection_record(methods.earth_properties),
        pipe_impedance = _selection_record(methods.pipe_impedance),
        temperature_dependence = _selection_record(methods.temperature_dependence),
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
    append!(paths,
        (
            "src/engine/internalimpedance/interface.jl",
            "src/engine/insulationimpedance/interface.jl",
            "src/engine/insulationadmittance/interface.jl",
            "src/engine/semiconadmittance/interface.jl",
            "src/engine/earthimpedance/interface.jl",
            "src/engine/earthimpedance/homogeneous.jl",
            "src/engine/earthadmittance/interface.jl",
            "src/engine/earthadmittance/homogeneous.jl",
            "src/engine/pipeimpedance/interface.jl",
            "src/earth/equivalenthomogeneous/interface.jl",
            "src/earth/frequencydependent/interface.jl"
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
        _ehem_path(methods.earth_impedance.equivalent_earth),
        _ehem_path(methods.earth_admittance.equivalent_earth)
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
