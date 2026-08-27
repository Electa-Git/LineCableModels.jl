const CASE_ROOT = joinpath(GAUNTLET_ROOT, "cases")
const CASE_INDEX_PATH = joinpath(CASE_ROOT, "index.toml")
const _CASE_IDENTIFIER = r"^[a-z][a-z0-9_]*$"

struct CaseParameter{T, Tags <: Tuple}
    id::Symbol
    nominal::T
    tags::Tags

    function CaseParameter(id::Symbol, nominal::T, tags::Tags) where {T, Tags <: Tuple}
        occursin(_CASE_IDENTIFIER, string(id)) || throw(ArgumentError(
            "case parameter identifiers must be lowercase; got $(repr(id))",
        ))
        isempty(tags) && throw(ArgumentError(
            "case parameter :$id must declare at least one tag",
        ))
        all(tag -> tag isa Symbol && occursin(_CASE_IDENTIFIER, string(tag)), tags) ||
            throw(ArgumentError(
                "case parameter tags must be lowercase symbols; got $(repr(tags))",
            ))
        length(unique(tags)) == length(tags) || throw(ArgumentError(
            "case parameter :$id contains duplicate tags",
        ))
        return new{T, Tags}(id, nominal, tags)
    end
end

case_parameter(id::Symbol, nominal; tags) =
    CaseParameter(id, nominal, Tuple(tags))

struct CaseDefinition{P <: NamedTuple, F}
    id::Symbol
    parameters::P
    build::F
    port_order::Vector{String}

    function CaseDefinition(
            id::Symbol,
            parameters::P,
            build::F,
            port_order::Vector{String}
    ) where {P <: NamedTuple, F}
        occursin(_CASE_IDENTIFIER, string(id)) || throw(ArgumentError(
            "case identifiers must be lowercase; got $(repr(id))",
        ))
        isempty(parameters) && throw(ArgumentError(
            "case :$id must declare at least one parameter",
        ))
        parameter_ids = Tuple(parameter.id for parameter in values(parameters))
        parameter_ids == keys(parameters) || throw(ArgumentError(
            "case :$id parameter field names must match their identifiers",
        ))
        length(unique(port_order)) == length(port_order) || throw(ArgumentError(
            "case :$id port order contains duplicate entries",
        ))
        isempty(port_order) && throw(ArgumentError(
            "case :$id must declare at least one terminal",
        ))
        return new{P, F}(id, parameters, build, port_order)
    end
end

function case_definition(
        build::F,
        id::Symbol,
        parameters::NamedTuple,
        port_order::AbstractVector{<:AbstractString}
) where {F}
    return CaseDefinition(id, parameters, build, String.(port_order))
end

abstract type AbstractCaseVariation end

struct NoVariation <: AbstractCaseVariation end

struct ExactOverrides{O <: NamedTuple} <: AbstractCaseVariation
    values::O
end
ExactOverrides(; values...) = ExactOverrides((; values...))

struct ParameterGrids{G <: NamedTuple} <: AbstractCaseVariation
    values::G

    function ParameterGrids(values::G) where {G <: NamedTuple}
        isempty(values) && throw(ArgumentError(
            "parameter grids must contain at least one source",
        ))
        all(source -> source isa Union{AbstractGrid, Gridspace}, Base.values(values)) ||
            throw(ArgumentError(
                "parameter-grid values must be Grid or Gridspace sources",
            ))
        return new{G}(values)
    end
end
ParameterGrids(; values...) = ParameterGrids((; values...))

struct RelativeStandardUncertainty{Tags <: Tuple} <: AbstractCaseVariation
    percent::Float64
    tags::Tags

    function RelativeStandardUncertainty(percent::Real; tags)
        isfinite(percent) && percent > 0 || throw(ArgumentError(
            "relative standard uncertainty must be positive and finite",
        ))
        selected_tags = Tuple(tags)
        isempty(selected_tags) && throw(ArgumentError(
            "relative standard uncertainty must select at least one tag",
        ))
        all(tag -> tag isa Symbol, selected_tags) || throw(ArgumentError(
            "relative standard uncertainty tags must be symbols",
        ))
        return new{typeof(selected_tags)}(Float64(percent), selected_tags)
    end
end

struct CompositeVariation{V <: Tuple} <: AbstractCaseVariation
    variations::V

    function CompositeVariation(variations::V) where {V <: Tuple}
        isempty(variations) && throw(ArgumentError(
            "a composite variation cannot be empty",
        ))
        all(variation -> variation isa AbstractCaseVariation, variations) ||
            throw(ArgumentError(
                "composite variations must contain case variation policies",
            ))
        return new{V}(variations)
    end
end

compose_variations(variations::AbstractCaseVariation...) =
    CompositeVariation(variations)

struct LoadedCase{D <: CaseDefinition, N, P, S <: NamedTuple, V}
    id::Symbol
    definition::D
    nominal_problem::N
    problem::P
    sources::S
    selected_parameters::Vector{Symbol}
    variation::V
    port_order::Vector{String}
    expected_size::NTuple{3, Int}
    source_file::String
    source_sha256::String
end

function _source_record(source)
    source isa Gridspace && return (
        kind = :gridspace,
        cardinality = length(source),
        type = string(typeof(source))
    )
    source isa AbstractGrid || return (kind = :exact, value = source)
    points = map(collect(source)) do point
        if point isa LineCableModels.UncertainValue
            return (
                nominal = LineCableModels.nominal(point),
                standard_uncertainty = LineCableModels.standard_uncertainty(point)
            )
        end
        return (value = point,)
    end
    return (
        kind = source isa LineCableModels.AbstractUncertainGrid ?
               :uncertain_grid : :deterministic_grid,
        cardinality = length(source),
        points
    )
end

variation_record(::NoVariation) = (kind = :none,)
variation_record(variation::ExactOverrides) = (
    kind = :exact_overrides,
    parameters = collect(keys(variation.values))
)
variation_record(variation::ParameterGrids) = (
    kind = :parameter_grids,
    parameters = collect(keys(variation.values))
)
variation_record(variation::RelativeStandardUncertainty) = (
    kind = :relative_standard_uncertainty,
    percent = variation.percent,
    tags = collect(variation.tags),
    interpretation = :standard_uncertainty
)
variation_record(variation::CompositeVariation) = (
    kind = :composite,
    variations = variation_record.(collect(variation.variations))
)

function parameter_manifest(model::LoadedCase)
    selected = Set(model.selected_parameters)
    return NamedTuple[
        (
            id = parameter.id,
            nominal = parameter.nominal,
            tags = collect(parameter.tags),
            selected = parameter.id in selected,
            source = _source_record(getproperty(model.sources, parameter.id))
        ) for parameter in values(model.definition.parameters)
    ]
end

function correlation_record(model::LoadedCase)
    return (
        rule = :parameter_identity,
        uncertain_primitives = copy(model.selected_parameters),
        note = "each selected parameter ID is one primitive; repeated builder uses are shared"
    )
end

function _case_index_document(path::AbstractString)
    isfile(path) || throw(ArgumentError("case index is missing: $path"))
    document = TOML.parsefile(path)
    haskey(document, "cases") || throw(ArgumentError(
        "case index must contain a [cases] table",
    ))
    entries = document["cases"]
    entries isa AbstractDict || throw(ArgumentError(
        "case index [cases] entry must be a table",
    ))
    return entries
end

function _case_source_path(relative::AbstractString; case_root::AbstractString = CASE_ROOT)
    isabspath(relative) && throw(ArgumentError(
        "case index paths must be relative; got $relative",
    ))
    root = realpath(case_root)
    candidate = normpath(joinpath(root, relative))
    separator = Base.Filesystem.path_separator
    startswith(candidate, root * separator) || throw(ArgumentError(
        "case index path escapes the cases directory: $relative",
    ))
    isfile(candidate) || throw(ArgumentError("indexed case file is missing: $candidate"))
    real = realpath(candidate)
    startswith(real, root * separator) || throw(ArgumentError(
        "indexed case file resolves outside the cases directory: $relative",
    ))
    return real
end

function case_index(
        ; index_path::AbstractString = CASE_INDEX_PATH,
        case_root::AbstractString = dirname(index_path)
)
    entries = _case_index_document(index_path)
    index = Dict{Symbol, String}()
    paths = Set{String}()
    for (raw_id, raw_path) in entries
        id = Symbol(raw_id)
        occursin(_CASE_IDENTIFIER, raw_id) || throw(ArgumentError(
            "case index identifier must be lowercase; got $(repr(raw_id))",
        ))
        raw_path isa AbstractString || throw(ArgumentError(
            "case index path for :$id must be a string",
        ))
        path = _case_source_path(raw_path; case_root)
        path in paths && throw(ArgumentError(
            "case index maps more than one identifier to $path",
        ))
        index[id] = path
        push!(paths, path)
    end
    isempty(index) && throw(ArgumentError("case index cannot be empty"))
    return index
end

_matches(parameter::CaseParameter, tags::Tuple) =
    all(tag -> tag in parameter.tags, tags)

function _apply_variation(
        ::NoVariation,
        parameter::CaseParameter,
        value,
        matched::Vector{Int},
        offset::Int
)
    return value
end

function _apply_variation(
        variation::ExactOverrides,
        parameter::CaseParameter,
        value,
        matched::Vector{Int},
        offset::Int
)
    haskey(variation.values, parameter.id) || return value
    matched[offset] += 1
    replacement = getproperty(variation.values, parameter.id)
    replacement isa Union{AbstractGrid, Gridspace} && throw(ArgumentError(
        "exact override :$(parameter.id) cannot be a Grid source; use ParameterGrids",
    ))
    return replacement
end

function _apply_variation(
        variation::ParameterGrids,
        parameter::CaseParameter,
        value,
        matched::Vector{Int},
        offset::Int
)
    haskey(variation.values, parameter.id) || return value
    matched[offset] += 1
    return getproperty(variation.values, parameter.id)
end

function _apply_variation(
        variation::RelativeStandardUncertainty,
        parameter::CaseParameter,
        value,
        matched::Vector{Int},
        offset::Int
)
    _matches(parameter, variation.tags) || return value
    value isa Real || throw(ArgumentError(
        "relative uncertainty selected non-real case parameter :$(parameter.id)",
    ))
    iszero(value) && throw(ArgumentError(
        "relative uncertainty selected zero-valued case parameter :$(parameter.id)",
    ))
    matched[offset] += 1
    return Grid(value, variation.percent)
end

function _apply_variation(
        variation::CompositeVariation,
        parameter::CaseParameter,
        value,
        matched::Vector{Int},
        offset::Int
)
    current = value
    position = offset
    for child in variation.variations
        current = _apply_variation(child, parameter, current, matched, position)
        position += _variation_count(child)
    end
    return current
end

_variation_count(::AbstractCaseVariation) = 1
_variation_count(variation::CompositeVariation) =
    sum(_variation_count, variation.variations; init = 0)

_variation_labels(::NoVariation) = ["no variation"]
_variation_labels(::ExactOverrides) = ["exact overrides"]
_variation_labels(::ParameterGrids) = ["parameter grids"]
function _variation_labels(variation::RelativeStandardUncertainty)
    ["relative standard uncertainty for tags $(variation.tags)"]
end
function _variation_labels(variation::CompositeVariation)
    collect(Iterators.flatten(_variation_labels(child) for child in variation.variations))
end

function _validate_variation(
        variation::AbstractCaseVariation,
        definition::CaseDefinition,
        matched::Vector{Int}
)
    labels = _variation_labels(variation)
    length(labels) == length(matched) || throw(AssertionError(
        "case variation accounting is inconsistent",
    ))
    for (label, count) in zip(labels, matched)
        label == "no variation" && continue
        count > 0 || throw(ArgumentError("$label matched no case parameters"))
    end
    overrides = _override_values(variation)
    unknown = filter(key -> !haskey(definition.parameters, key), keys(overrides))
    isempty(unknown) || throw(ArgumentError(
        "unknown case parameter overrides: $(sort!(collect(unknown)))",
    ))
    return nothing
end

_override_values(::AbstractCaseVariation) = (;)
_override_values(variation::ExactOverrides) = variation.values
_override_values(variation::ParameterGrids) = variation.values
function _override_values(variation::CompositeVariation)
    return merge((_override_values(child) for child in variation.variations)...)
end

function _case_sources(definition::CaseDefinition, variation::AbstractCaseVariation)
    count = _variation_count(variation)
    matched = zeros(Int, count)
    source_values = map(definition.parameters) do parameter
        _apply_variation(variation, parameter, parameter.nominal, matched, 1)
    end
    _validate_variation(variation, definition, matched)
    return source_values
end

function _materialize_case(
        definition::CaseDefinition,
        sources::NamedTuple,
        nominal_problem
)
    any(source -> source isa Union{AbstractGrid, Gridspace}, values(sources)) ||
        return Base.invokelatest(definition.build, sources)
    names = keys(sources)
    grids = map(values(sources)) do source
        source isa Union{AbstractGrid, Gridspace} ? source : Grid((source,))
    end
    materializer = function (args...)
        Base.invokelatest(definition.build, NamedTuple{names}(args))
    end
    return Gridspace{LineCableModels.Engine.LineParametersProblem}(
        materializer,
        grids
    )
end

function _validate_loaded_problem(
        definition::CaseDefinition,
        nominal_problem::LineCableModels.Engine.LineParametersProblem
)
    nominal_problem.system.system_id == string(definition.id) || throw(ArgumentError(
        "case :$(definition.id) built system $(repr(nominal_problem.system.system_id))",
    ))
    assignments = collect(Iterators.flatten(
        position.conn for position in nominal_problem.system.cables
    ))
    isempty(assignments) && throw(ArgumentError(
        "case :$(definition.id) requires at least one explicit terminal",
    ))
    any(iszero, assignments) && throw(ArgumentError(
        "case :$(definition.id) may not eliminate a declared terminal",
    ))
    sort(assignments) == collect(1:length(assignments)) || throw(ArgumentError(
        "case :$(definition.id) phase assignments must be contiguous from one",
    ))
    length(unique(assignments)) == length(assignments) || throw(ArgumentError(
        "case :$(definition.id) may not bundle declared terminals",
    ))
    length(definition.port_order) == length(assignments) || throw(DimensionMismatch(
        "case :$(definition.id) port order does not match its terminals",
    ))
    return (length(assignments), length(assignments), length(nominal_problem.frequencies))
end

function load_case(
        id::Symbol;
        variation::AbstractCaseVariation = NoVariation(),
        index_path::AbstractString = CASE_INDEX_PATH,
        case_root::AbstractString = dirname(index_path)
)
    index = case_index(; index_path, case_root)
    haskey(index, id) || throw(ArgumentError("unknown Gauntlet case :$id"))
    source_file = index[id]
    definition = Base.include(@__MODULE__, source_file)
    definition isa CaseDefinition || throw(ArgumentError(
        "case file $source_file must evaluate to one CaseDefinition",
    ))
    definition.id === id || throw(ArgumentError(
        "case file $source_file declares :$(definition.id), expected :$id",
    ))
    nominal_values = map(parameter -> parameter.nominal, definition.parameters)
    nominal_problem = Base.invokelatest(definition.build, nominal_values)
    nominal_problem isa LineCableModels.Engine.LineParametersProblem ||
        throw(ArgumentError(
            "case :$id must build an Engine.LineParametersProblem",
        ))
    expected_size = _validate_loaded_problem(definition, nominal_problem)
    sources = _case_sources(definition, variation)
    problem = _materialize_case(definition, sources, nominal_problem)
    loaded_size = problem isa LineCableModels.Engine.LineParametersProblem ?
                  _validate_loaded_problem(definition, problem) : expected_size
    selected = Symbol[
        parameter.id for (parameter, source) in zip(values(definition.parameters), values(sources))
        if !isequal(source, parameter.nominal)
    ]
    return LoadedCase(
        id,
        definition,
        nominal_problem,
        problem,
        sources,
        selected,
        variation,
        copy(definition.port_order),
        loaded_size,
        source_file,
        bytes2hex(sha256(read(source_file)))
    )
end
