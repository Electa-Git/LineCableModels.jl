"Concrete callable representation of a type constructor used by a Gridspace."
struct _TypeConstructor{Target} end

(::_TypeConstructor{Target})(arguments...) where {Target} = Target(arguments...)

_gridspace_callable(build) = build
_gridspace_callable(::Type{Target}) where {Target} = _TypeConstructor{Target}()

"""
$(TYPEDEF)

Represent a typed finite space assembled from explicit [`Grid`](@ref), nested
`Gridspace`, or admitted completed result-space sources. `combine` is encoded
in the type and may be `:product` or `:zip`. `Target` identifies the semantic
result family. A nonempty deterministic space records the concrete type
returned by its callable when Julia can prove it without evaluating a point.
Otherwise the space declares its iterator element type unknown until values
are materialised; it never advertises a `UnionAll` result type.

$(TYPEDFIELDS)
"""
struct Gridspace{Target, F, G <: Tuple, C, R}
    "Callable that constructs `Target` from one selected argument tuple."
    build::F

    "Explicit Grid or nested Gridspace sources."
    grids::G
end

"Sentinel result type for a Gridspace whose materialised element type is unknown."
struct _UnknownGridspaceEltype end

const _GridspaceSource = Union{AbstractGrid, Gridspace, AbstractResultSpace}

function _gridspace_eltype(::Type{Tuple}, ::typeof(tuple), grids::Tuple)
    (any(has_uncertainty, grids) || any(grid -> iszero(length(grid)), grids)) &&
        return _UnknownGridspaceEltype
    argument_types = map(grid -> eltype(typeof(grid)), grids)
    all(isconcretetype, argument_types) || return _UnknownGridspaceEltype
    return Tuple{argument_types...}
end

function _gridspace_eltype(::Type{Target}, build, grids::Tuple) where {Target}
    (any(has_uncertainty, grids) || any(grid -> iszero(length(grid)), grids)) &&
        return _UnknownGridspaceEltype

    argument_types = map(grid -> eltype(typeof(grid)), grids)
    inferred = Base.code_typed(
        build,
        Tuple{argument_types...};
        optimize = false
    )
    result_type = isempty(inferred) ? Any :
                  foldl(typejoin, (last(entry) for entry in inferred))
    isconcretetype(result_type) || return _UnknownGridspaceEltype
    result_type <: Target || throw(ArgumentError(
        "Gridspace callable result $result_type is not a subtype of target $Target",
    ))
    return result_type
end

"""
$(TYPEDSIGNATURES)

Construct a lazy space from explicit finite sources. Product composition uses
Julia's product order, with the first source varying fastest. Zip composition
pairs equal-length sources and broadcasts singleton sources. A completed
result space is admitted only through a result-owner transport method of the
form `Gridspace{Target}(source)`.

# Errors

- Throws `ArgumentError` when `combine` is unsupported or an input is not a
  `Grid` or nested `Gridspace`.
- Throws `ArgumentError` when an inferred concrete callable result does not
  belong to `Target`.
- Throws `DimensionMismatch` when zip source lengths are incompatible.
"""
function Gridspace{Target}(
        build,
        grids::G;
        combine::Symbol = :product
) where {Target, G <: Tuple}
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip; got :$combine"))
    all(grid -> grid isa _GridspaceSource, grids) || throw(
        ArgumentError(
        "Gridspace sources must be Grid, Gridspace, or completed result-space values",
    ),
    )
    callable = _gridspace_callable(build)
    result_type = _gridspace_eltype(Target, callable, grids)
    space = Gridspace{
        Target, typeof(callable), G, Val{combine}, result_type
    }(callable, grids)
    combine === :zip && _zip_length(space)
    return space
end

function Gridspace{Target}(grids::Tuple; combine::Symbol = :product) where {Target}
    Gridspace{Target}(Target, grids; combine)
end

Grid(space::Gridspace) = space

const _FiniteSource = Union{AbstractGrid, Gridspace}
const _FiniteCollection = Union{Tuple, AbstractVector}

_construction_axis(::_FiniteSource) = true
_construction_axis(values::_FiniteCollection) = any(_construction_axis, values)

_vector(values...) = collect(values)

_finite_source(source::_FiniteSource; combine::Symbol) = source
function _finite_source(values::Tuple; combine::Symbol)
    sources = map(values) do value
        _construction_axis(value) ? _finite_source(value; combine) : Grid((value,))
    end
    return Gridspace{Tuple}(tuple, sources; combine)
end
function _finite_source(values::AbstractVector; combine::Symbol)
    sources = map(values) do value
        _construction_axis(value) ? _finite_source(value; combine) : Grid((value,))
    end
    return Gridspace{Vector}(_vector, Tuple(sources); combine)
end

function _finite_construction(
        ::Type{Target}, caller, values::Tuple; combine::Symbol
) where {Target}
    sources = map(values) do value
        _construction_axis(value) ? _finite_source(value; combine) : Grid((value,))
    end
    return Gridspace{Target}(caller, sources; combine)
end

points(grid::AbstractGrid) = grid
points(result::AbstractResultSpace) = result

_product_combinations(::Tuple{}) = ((),)
function _product_combinations(grids::Tuple)
    Iterators.product(map(points, grids)...)
end

_zip_length(::Gridspace{<:Any, <:Any, Tuple{}, Val{:zip}, <:Any}) = 1
function _zip_length(space::Gridspace{<:Any, <:Any, <:Tuple, Val{:zip}, <:Any})
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

_zip_combinations(::Gridspace{<:Any, <:Any, Tuple{}, Val{:zip}, <:Any}) = ((),)
function _zip_combinations(
        space::Gridspace{<:Any, <:Any, <:Tuple, Val{:zip}, <:Any}
)
    target_count = _zip_length(space)
    iterators = map(source -> _zip_source(source, target_count), space.grids)
    return Iterators.zip(iterators...)
end

function _combinations(
        space::Gridspace{<:Any, <:Any, <:Any, Val{:product}, <:Any}
)
    _product_combinations(space.grids)
end
function _combinations(space::Gridspace{<:Any, <:Any, <:Any, Val{:zip}, <:Any})
    _zip_combinations(space)
end

"Return the lazy unresolved points of a Gridspace."
function points(space::Gridspace{Target}) where {Target}
    (Gridpoint{Target}(space.build, args) for args in _combinations(space))
end

"Return a value unchanged during deterministic point materialisation."
materialize(value) = value

function materialize(value::UncertainValue)
    throw(ArgumentError(
        "direct materialisation of an uncertainty-bearing Gridspace requires " *
        "Measurements.jl; load it with `using Measurements` before " *
        "iteration, or use stochastic realisation through MonteCarlo",
    ))
end

"Recursively materialise a selected Gridspace point."
function materialize(point::Gridpoint)
    point.build(map(materialize, point.args)...)
end

function materialize(
        point::Gridpoint{Target}
) where {Target <: AbstractProblemDefinition}
    return validate(point.build(map(materialize, point.args)...))::Target
end

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
realize(point::Gridpoint, arguments::Tuple) = point.build(arguments...)

function realize(
        point::Gridpoint{Target},
        arguments::Tuple
) where {Target <: AbstractProblemDefinition}
    return validate(point.build(arguments...))::Target
end

"Recursively realise a selected Gridspace point using the caller's RNG."
function realize(rng::Random.AbstractRNG, point::Gridpoint, distribution)
    return realize(point, realize_arguments(rng, point, distribution))
end

function _materialized_result(
        ::Type{_UnknownGridspaceEltype}, ::Type{Target}, point
) where {Target}
    return materialize(point)::Target
end
function _materialized_result(::Type{Result}, ::Type, point) where {Result}
    return materialize(point)::Result
end

function Base.iterate(
        space::Gridspace{Target, <:Any, <:Any, <:Any, Result}, state...
) where {Target, Result}
    item = iterate(points(space), state...)
    item === nothing && return nothing
    point, next_state = item
    return _materialized_result(Result, Target, point), next_state
end

_gridspace_iterator_eltype(::Type{_UnknownGridspaceEltype}) = Base.IteratorEltype(Any)
_gridspace_iterator_eltype(::Type) = Base.IteratorEltype(Vector{Int})

_gridspace_eltype_trait(::Type{_UnknownGridspaceEltype}) = Any
_gridspace_eltype_trait(::Type{Result}) where {Result} = Result

#! explicit-imports: off
# Base's iterator trait protocol exposes these values without public bindings.
Base.IteratorSize(::Type{<:Gridspace}) = Base.HasShape{1}()
function Base.IteratorEltype(
        ::Type{<:Gridspace{<:Any, <:Any, <:Any, <:Any, Result}}
) where {Result}
    return _gridspace_iterator_eltype(Result)
end
#! explicit-imports: on
function Base.eltype(
        ::Type{<:Gridspace{<:Any, <:Any, <:Any, <:Any, Result}}
) where {Result}
    _gridspace_eltype_trait(Result)
end
function Base.length(
        space::Gridspace{<:Any, <:Any, <:Any, Val{:product}, <:Any}
)
    prod(length, space.grids; init = 1)
end
Base.length(space::Gridspace{<:Any, <:Any, <:Any, Val{:zip}, <:Any}) = _zip_length(space)
Base.size(space::Gridspace) = (length(space),)

"Draw and realise one product point by sampling each source once."
function Base.rand(
        rng::Random.AbstractRNG,
        space::Gridspace{<:Any, <:Any, <:Any, Val{:product}, Result};
        distribution = :normal
) where {Result}
    isempty(space) && throw(ArgumentError("cannot sample an empty Gridspace"))
    arguments = map(
        source -> rand(rng, source; distribution),
        space.grids
    )
    value = space.build(arguments...)
    Result === _UnknownGridspaceEltype && return value
    return value::Result
end

"Draw one zipped point and realise its aligned sources without collecting."
function Base.rand(
        rng::Random.AbstractRNG,
        space::Gridspace{<:Any, <:Any, <:Any, Val{:zip}, Result};
        distribution = :normal
) where {Result}
    isempty(space) && throw(ArgumentError("cannot sample an empty Gridspace"))
    offset = rand(rng, 0:(length(space) - 1))
    arguments = map(space.grids) do source
        source_offset = length(source) == 1 ? 0 : offset
        selected = first(Iterators.drop(points(source), source_offset))
        realize(rng, selected, distribution)
    end
    value = space.build(arguments...)
    Result === _UnknownGridspaceEltype && return value
    return value::Result
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
