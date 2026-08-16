"""
$(TYPEDEF)

Typed blueprint interface for staged construction of a concrete `Target`.
"""
abstract type AbstractSpec{Target} end

target_type(::Type{<:AbstractSpec{Target}}) where {Target} = Target
target_type(spec::AbstractSpec) = target_type(typeof(spec))

struct ConstantAxis{T}
    value::T
end

Base.iterate(axis::ConstantAxis) = (axis.value, nothing)
Base.iterate(::ConstantAxis, ::Nothing) = nothing
Base.length(::ConstantAxis) = 1
Base.getindex(axis::ConstantAxis, index::Integer) =
    index == 1 ? axis.value : throw(BoundsError(axis, index))

struct GridBinding{K}
    key::K
    index::Int
    cardinality::Int
end

"""
$(TYPEDEF)

Represent one resolved deterministic choice while retaining any uncertainty
descriptors for later realization.

$(TYPEDFIELDS)
"""
struct Configuration{Target,F,V<:Tuple,N<:Tuple,B<:Tuple}
    "Callable that constructs `Target` from the resolved axis values."
    target::F

    "Resolved values supplied to `target`."
    values::V

    "Parameter names corresponding to `values`."
    names::N

    "Selections retained for coupled axes."
    bindings::B
end

target_type(::Type{<:Configuration{Target}}) where {Target} = Target
target_type(configuration::Configuration) = target_type(typeof(configuration))

"""
$(TYPEDEF)

A lazy space of complete `Target` configurations. `combine` is local to this
node and may be `:product` or `:zip`.

$(TYPEDFIELDS)
"""
struct Gridspace{Target,F,A<:Tuple,N<:Tuple,C} <: AbstractSpec{Target}
    "Callable that constructs `Target` from one selection of the direct axes."
    target::F

    "Direct parameter or object-valued axes."
    axes::A

    "Parameter names corresponding to `axes`."
    names::N

    "Local composition rule, represented by `Val{:product}` or `Val{:zip}`."
    combine::C
end

"""
$(TYPEDSIGNATURES)

Construct a lazy space of complete `Target` configurations from `axes`.

`combine=:product` forms the Cartesian product of the direct axes.
`combine=:zip` pairs axes of equal cardinality while treating singleton axes
as constants. A nested [`Gridspace`](@ref) applies its own composition rule
before entering its parent as one object-valued axis.
"""
function Gridspace{Target}(
    target::F,
    axes::A,
    names::N=();
    combine::Symbol=:product,
) where {Target,F,A<:Tuple,N<:Tuple}
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip; got :$combine"))
    isempty(names) || length(names) == length(axes) ||
        throw(DimensionMismatch("Gridspace names and axes must have equal lengths"))
    normalized_axes = map(_gridspace_axis, axes)
    normalized_names = isempty(names) ?
        ntuple(index -> Symbol(:arg, index), length(axes)) : names
    return Gridspace{
        Target,
        F,
        typeof(normalized_axes),
        typeof(normalized_names),
        Val{combine},
    }(
        target,
        normalized_axes,
        normalized_names,
        Val(combine),
    )
end

Gridspace{Target}(axes::Tuple; combine::Symbol=:product) where {Target} =
    Gridspace{Target}(Target, axes; combine)

Grid(space::Gridspace; key=nothing) = key === nothing ? space :
    throw(ArgumentError("Gridspace coupling is defined by its child Grids"))

_gridspace_axis(value::ConstantAxis) = value
_gridspace_axis(value::Union{AbstractGrid,AbstractSpec}) = value
_gridspace_axis(value) = ConstantAxis(value)

struct AxisSelection{V,K}
    value::V
    key::K
    index::Int
    cardinality::Int
end

"""A resolved axis value retaining its coupling identity for realization."""
struct ResolvedGridValue{V,K}
    value::V
    key::K
end

_axis_cases(axis::ConstantAxis) = axis

function _axis_cases(grid::AbstractGrid)
    return (
        AxisSelection(value, grid.key, index, length(grid))
        for (index, value) in enumerate(grid)
    )
end

_axis_cases(spec::AbstractSpec) = configurations(spec)

_axis_value(selection::AxisSelection) =
    ResolvedGridValue(selection.value, selection.key)
_axis_value(configuration::Configuration) = configuration
_axis_value(value) = value

_axis_bindings(selection::AxisSelection) =
    (GridBinding(selection.key, selection.index, selection.cardinality),)
_axis_bindings(configuration::Configuration) = configuration.bindings
_axis_bindings(::Any) = ()

_same_grid_key(left, right) = left == right

function _compatible_bindings(items::Tuple)
    bindings = tuple((_axis_bindings(item) for item in items)...)
    flat = tuple(Iterators.flatten(bindings)...)
    for i in eachindex(flat), j in (i + 1):length(flat)
        if _same_grid_key(flat[i].key, flat[j].key)
            flat[i].cardinality == flat[j].cardinality || throw(DimensionMismatch(
                "coupled Grids have incompatible cardinalities $(flat[i].cardinality) and $(flat[j].cardinality)",
            ))
            flat[i].index != flat[j].index && return false
        end
    end
    return true
end

function _merged_bindings(items::Tuple)
    groups = tuple((_axis_bindings(item) for item in items)...)
    all_bindings = tuple(Iterators.flatten(groups)...)
    return tuple((
        binding for (index, binding) in pairs(all_bindings)
        if all(
            previous -> !_same_grid_key(previous.key, binding.key),
            all_bindings[1:(index - 1)],
        )
    )...)
end

function _product_combinations(axes::Tuple)
    isempty(axes) && return ((),)
    iterators = map(_axis_cases, axes)
    return Iterators.product(iterators...)
end

function _nth(iterator, index::Int)
    item = iterate(Iterators.drop(iterator, index - 1))
    item === nothing && throw(BoundsError(iterator, index))
    return item[1]
end

function _zip_combinations(axes::Tuple, names::Tuple)
    isempty(axes) && return ((),)
    iterators = map(_axis_cases, axes)
    counts = map(length, axes)
    target_count = maximum(counts)
    for index in eachindex(counts)
        counts[index] in (1, target_count) || throw(DimensionMismatch(
            "zip axis $(names[index]) has cardinality $(counts[index]); expected 1 or $target_count",
        ))
    end
    return (
        map(
            (iterator, count) -> _nth(iterator, count == 1 ? 1 : row),
            iterators,
            counts,
        )
        for row in 1:target_count
    )
end

_combinations(space::Gridspace{<:Any,<:Any,<:Any,<:Any,Val{:product}}) =
    _product_combinations(space.axes)
_combinations(space::Gridspace{<:Any,<:Any,<:Any,<:Any,Val{:zip}}) =
    _zip_combinations(space.axes, space.names)

"""
$(TYPEDSIGNATURES)

Return a lazy iterator over the resolved configurations admitted by a
[`Gridspace`](@ref).
"""
function configurations(space::Gridspace{Target}) where {Target}
    compatible = Iterators.filter(_compatible_bindings, _combinations(space))
    return (
        Configuration{Target}(
            space.target,
            map(_axis_value, items),
            space.names,
            _merged_bindings(items),
        )
        for items in compatible
    )
end

# Partial parameter constructor used above while retaining fully concrete field
# parameters in the actual value.
function Configuration{Target}(target::F, values::V, names::N, bindings::B) where {
    Target,F,V<:Tuple,N<:Tuple,B<:Tuple
}
    return Configuration{Target,F,V,N,B}(target, values, names, bindings)
end

"""
$(TYPEDSIGNATURES)

Return the [`Gridspace`](@ref) backing a staged specification. Concrete
specifications implement this method.
"""
function gridspace(spec::AbstractSpec)
    throw(MethodError(gridspace, (spec,)))
end

configurations(spec::AbstractSpec) = configurations(gridspace(spec))

_direct_value(value::UncertainValue) = throw(ArgumentError(
    "direct materialization of uncertain configurations requires Measurements.jl; load it before using FullParametric",
))
_direct_value(configuration::Configuration) = materialize(configuration)
_direct_value(value) = value

function _resolved_direct(value::ResolvedGridValue, cache::Dict)
    return get!(cache, value.key) do
        _direct_value(value.value)
    end
end
_resolved_direct(configuration::Configuration, cache::Dict) =
    _materialize(configuration, cache)
_resolved_direct(value, ::Dict) = _direct_value(value)

function _materialize(configuration::Configuration, cache::Dict)
    values = map(value -> _resolved_direct(value, cache), configuration.values)
    return configuration.target(values...)
end

"""
$(TYPEDSIGNATURES)

Materialize a resolved configuration through its target constructor.
"""
function materialize(configuration::Configuration)
    return _materialize(configuration, Dict{Any,Any}())
end

function _random_value(
    rng::Random.AbstractRNG,
    value::ResolvedGridValue,
    distribution,
    cache::Dict,
)
    return get!(cache, value.key) do
        value.value isa UncertainValue ?
            rand(rng, value.value; distribution) : value.value
    end
end
_random_value(
    rng::Random.AbstractRNG,
    configuration::Configuration,
    distribution,
    cache::Dict,
) =
    _random_materialize(rng, configuration, distribution, cache)
_random_value(::Random.AbstractRNG, value, distribution, ::Dict) = value

function _random_materialize(rng, configuration, distribution, cache)
    values = map(
        value -> _random_value(rng, value, distribution, cache),
        configuration.values,
    )
    return configuration.target(values...)
end

function Base.rand(
    rng::Random.AbstractRNG,
    configuration::Configuration;
    distribution=:normal,
)
    return _random_materialize(
        rng,
        configuration,
        distribution,
        Dict{Any,Any}(),
    )
end

Base.rand(configuration::Configuration; kwargs...) =
    rand(Random.default_rng(), configuration; kwargs...)

function Base.rand(rng::Random.AbstractRNG, spec::AbstractSpec; distribution=:normal)
    iterator = configurations(spec)
    first_item = iterate(iterator)
    first_item === nothing && throw(ArgumentError("cannot sample an empty Gridspace"))
    configuration, state = first_item
    iterate(iterator, state) === nothing || throw(ArgumentError(
        "rand(Gridspace) requires exactly one outer configuration; enumerate configurations and sample one explicitly",
    ))
    return rand(rng, configuration; distribution)
end

Base.rand(spec::AbstractSpec; kwargs...) =
    rand(Random.default_rng(), spec; kwargs...)

function Base.iterate(spec::AbstractSpec)
    iterator = configurations(spec)
    item = iterate(iterator)
    item === nothing && return nothing
    configuration, state = item
    return materialize(configuration), state
end

function Base.iterate(spec::AbstractSpec, state)
    iterator = configurations(spec)
    item = iterate(iterator, state)
    item === nothing && return nothing
    configuration, next_state = item
    return materialize(configuration), next_state
end

Base.IteratorSize(::Type{<:AbstractSpec}) = Base.HasLength()
Base.IteratorEltype(::Type{<:AbstractSpec}) = Base.HasEltype()
Base.eltype(::Type{<:AbstractSpec{Target}}) where {Target} = Target
Base.length(spec::AbstractSpec) = count(_ -> true, configurations(spec))
Base.size(spec::AbstractSpec) = (length(spec),)
Base.getindex(spec::AbstractSpec, index::Integer) =
    first(Iterators.drop(spec, index - 1))

"""
$(TYPEDSIGNATURES)

Return whether a value, configuration, or specification contains an
uncertainty descriptor.
"""
has_uncertainty(value::UncertainValue) = true
has_uncertainty(value::ResolvedGridValue) = has_uncertainty(value.value)
has_uncertainty(configuration::Configuration) =
    any(has_uncertainty, configuration.values)
has_uncertainty(::Any) = false
has_uncertainty(spec::AbstractSpec) =
    any(has_uncertainty, configurations(spec))

_manifest_value(value::UncertainValue) = (
    nominal=value.nominal,
    sigma=value.sigma,
    style=value.style isa RelativeUncertainty ? :relative : :absolute,
)
_manifest_value(value::ResolvedGridValue) = _manifest_value(value.value)
_manifest_value(configuration::Configuration) =
    configuration_manifest(configuration)
_manifest_value(value) = value

"""
$(TYPEDSIGNATURES)

Return the named resolved parameterization of a configuration in a stable,
serializable form.
"""
function configuration_manifest(configuration::Configuration)
    values = map(_manifest_value, configuration.values)
    return NamedTuple{configuration.names}(values)
end
