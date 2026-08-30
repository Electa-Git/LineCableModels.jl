"A selected but unresolved argument tuple from a [`Gridspace`](@ref)."
struct Gridpoint{F, A <: Tuple}
    build::F
    args::A
end

"""
$(TYPEDEF)

Represent a typed finite space assembled from explicit [`Grid`](@ref) or
nested `Gridspace` sources. `combine` is encoded in the type and may be
`:product` or `:zip`.

$(TYPEDFIELDS)
"""
struct Gridspace{Target, F, G <: Tuple, C}
    "Callable that constructs `Target` from one selected argument tuple."
    build::F

    "Explicit Grid or nested Gridspace sources."
    grids::G
end

"""
$(TYPEDSIGNATURES)

Construct a lazy space from explicit finite sources. Product composition uses
Julia's product order, with the first source varying fastest. Zip composition
pairs equal-length sources and broadcasts singleton sources.

# Errors

- Throws `ArgumentError` when `combine` is unsupported or an input is not a
  `Grid` or nested `Gridspace`.
- Throws `DimensionMismatch` when zip source lengths are incompatible.
"""
function Gridspace{Target}(
        build::F,
        grids::G;
        combine::Symbol = :product
) where {Target, F, G <: Tuple}
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip; got :$combine"))
    all(grid -> grid isa Union{AbstractGrid, Gridspace}, grids) || throw(
        ArgumentError("Gridspace sources must be Grid or Gridspace values"),
    )
    space = Gridspace{Target, F, G, Val{combine}}(build, grids)
    combine === :zip && _zip_length(space)
    return space
end

function Gridspace{Target}(grids::Tuple; combine::Symbol = :product) where {Target}
    Gridspace{Target}(Target, grids; combine)
end

Grid(space::Gridspace) = space

const _FiniteSource = Union{AbstractGrid, Gridspace}

_construction_axis(::_FiniteSource) = true

function _finite_construction(
        ::Type{Target}, caller, values::Tuple; combine::Symbol
) where {Target}
    sources = map(values) do value
        _construction_axis(value) ? value : Grid((value,))
    end
    return Gridspace{Target}(caller, sources; combine)
end

points(grid::AbstractGrid) = grid

_product_combinations(::Tuple{}) = ((),)
function _product_combinations(grids::Tuple)
    Iterators.product(map(points, grids)...)
end

_zip_length(::Gridspace{<:Any, <:Any, Tuple{}, Val{:zip}}) = 1
function _zip_length(space::Gridspace{<:Any, <:Any, <:Tuple, Val{:zip}})
    counts = map(length, space.grids)
    target_count = maximum(counts)
    all(count -> count == 1 || count == target_count, counts) || throw(
        DimensionMismatch(
        "zip sources must have equal cardinality or be singletons; got $(Tuple(counts))",
    ),
    )
    return target_count
end

function _zip_source(source, target_count::Int)
    length(source) == 1 &&
        return Iterators.repeated(only(points(source)), target_count)
    return points(source)
end

_zip_combinations(::Gridspace{<:Any, <:Any, Tuple{}, Val{:zip}}) = ((),)
function _zip_combinations(space::Gridspace{<:Any, <:Any, <:Tuple, Val{:zip}})
    target_count = _zip_length(space)
    iterators = map(source -> _zip_source(source, target_count), space.grids)
    return Iterators.zip(iterators...)
end

function _combinations(space::Gridspace{<:Any, <:Any, <:Any, Val{:product}})
    _product_combinations(space.grids)
end
function _combinations(space::Gridspace{<:Any, <:Any, <:Any, Val{:zip}})
    _zip_combinations(space)
end

"Return the lazy unresolved points of a Gridspace."
function points(space::Gridspace)
    (Gridpoint(space.build, args) for args in _combinations(space))
end

"Return a value unchanged during deterministic point materialisation."
materialize(value) = value

function materialize(value::UncertainValue)
    throw(ArgumentError(
        "materialising uncertainty descriptors requires Measurements.jl",
    ))
end

"Recursively materialise a selected Gridspace point."
function materialize(point::Gridpoint)
    point.build(map(materialize, point.args)...)
end

"Return every uncertainty descriptor in recursive realization order."
uncertainties(::Any) = ()
uncertainties(value::UncertainValue) = (value,)

_point_uncertainties(::Tuple{}) = ()
function _point_uncertainties(arguments::Tuple)
    return (
        uncertainties(first(arguments))...,
        _point_uncertainties(Base.tail(arguments))...
    )
end

uncertainties(point::Gridpoint) = _point_uncertainties(point.args)

"Return a deterministic value unchanged during stochastic realisation."
realize(::Random.AbstractRNG, value, _) = value

function realize(rng::Random.AbstractRNG, value::UncertainValue, distribution)
    rand(rng, value; distribution)
end

"Draw the arguments of a selected Gridspace point without invoking its builder."
function realize_arguments(rng::Random.AbstractRNG, point::Gridpoint, distribution)
    return map(value -> realize(rng, value, distribution), point.args)
end

"Build a selected Gridspace point from an already realised argument tuple."
build(point::Gridpoint, arguments::Tuple) = point.build(arguments...)

_realize_supplied(value, values::Tuple) = (value, values)
function _realize_supplied(::UncertainValue, values::Tuple)
    return first(values), Base.tail(values)
end

_realize_supplied_arguments(::Tuple{}, values::Tuple) = ((), values)
function _realize_supplied_arguments(arguments::Tuple, values::Tuple)
    value, remaining = _realize_supplied(first(arguments), values)
    tail, remaining = _realize_supplied_arguments(
        Base.tail(arguments),
        remaining
    )
    return ((value, tail...), remaining)
end

function _realize_supplied(point::Gridpoint, values::Tuple)
    arguments, remaining = _realize_supplied_arguments(point.args, values)
    return build(point, arguments), remaining
end

"Build a selected Gridspace point from its previously realised arguments."
realize(point::Gridpoint, arguments::Tuple) = build(point, arguments)

"Recursively realise a selected Gridspace point from physical scalar values."
function realize(point::Gridpoint, values::Tuple{Vararg{Real}})
    expected = length(uncertainties(point))
    length(values) == expected || throw(DimensionMismatch(
        "supplied realization values must match the $expected uncertainty coordinates",
    ))
    realized, remaining = _realize_supplied(point, values)
    isempty(remaining) || throw(DimensionMismatch(
        "supplied realization values were not consumed exactly",
    ))
    return realized
end

"Recursively realise a selected Gridspace point using the caller's RNG."
function realize(rng::Random.AbstractRNG, point::Gridpoint, distribution)
    return build(point, realize_arguments(rng, point, distribution))
end

function Base.iterate(space::Gridspace{Target}, state...) where {Target}
    item = iterate(points(space), state...)
    item === nothing && return nothing
    point, next_state = item
    return materialize(point)::Target, next_state
end

#! explicit-imports: off
# Base's iterator trait protocol exposes these values without public bindings.
Base.IteratorSize(::Type{<:Gridspace}) = Base.HasShape{1}()
Base.IteratorEltype(::Type{<:Gridspace}) = Base.HasEltype()
#! explicit-imports: on
Base.eltype(::Type{<:Gridspace{Target}}) where {Target} = Target
function Base.length(space::Gridspace{<:Any, <:Any, <:Any, Val{:product}})
    prod(length, space.grids; init = 1)
end
Base.length(space::Gridspace{<:Any, <:Any, <:Any, Val{:zip}}) = _zip_length(space)
Base.size(space::Gridspace) = (length(space),)

"Draw and realise one product point by sampling each source once."
function Base.rand(
        rng::Random.AbstractRNG,
        space::Gridspace{<:Any, <:Any, <:Any, Val{:product}};
        distribution = :normal
)
    isempty(space) && throw(ArgumentError("cannot sample an empty Gridspace"))
    arguments = map(
        source -> rand(rng, source; distribution),
        space.grids
    )
    return space.build(arguments...)
end

"Draw one zipped point and realise its aligned sources without collecting."
function Base.rand(
        rng::Random.AbstractRNG,
        space::Gridspace{<:Any, <:Any, <:Any, Val{:zip}};
        distribution = :normal
)
    isempty(space) && throw(ArgumentError("cannot sample an empty Gridspace"))
    offset = rand(rng, 0:(length(space) - 1))
    arguments = map(space.grids) do source
        source_offset = length(source) == 1 ? 0 : offset
        selected = first(Iterators.drop(points(source), source_offset))
        realize(rng, selected, distribution)
    end
    return space.build(arguments...)
end

function Base.rand(space::Gridspace; distribution = :normal)
    return rand(Random.default_rng(), space; distribution)
end

"Return whether a value structurally contains an uncertainty descriptor."
has_uncertainty(::UncertainValue) = true
has_uncertainty(grid::AbstractUncertainGrid) = true
has_uncertainty(grid::DeterministicGrid) = any(has_uncertainty, grid.vals)
has_uncertainty(point::Gridpoint) = any(has_uncertainty, point.args)
has_uncertainty(space::Gridspace) = any(has_uncertainty, space.grids)
has_uncertainty(::Any) = false
