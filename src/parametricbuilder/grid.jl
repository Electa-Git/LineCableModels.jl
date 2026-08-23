import Base: eltype, extrema, getindex, iterate, length, rand, size
import Random

"""
$(TYPEDEF)

Supertype for the deterministic and uncertainty-bearing parameter axes created
by [`Grid`](@ref).
"""
abstract type AbstractGrid end

"""
$(TYPEDEF)

Supertype for parameter axes whose configurations retain an uncertainty
descriptor.
"""
abstract type AbstractUncertainGrid <: AbstractGrid end

"""
$(TYPEDEF)

A dependency-free uncertainty descriptor. `sigma` is always an absolute
standard deviation.

$(TYPEDFIELDS)
"""
struct UncertainValue{T, E}
    "Nominal parameter value."
    nominal::T

    "Absolute standard uncertainty in the same unit as `nominal`."
    sigma::E

    function UncertainValue(nominal::T, sigma::E) where {T, E}
        if nominal isa Real && sigma isa Real
            isfinite(nominal) || throw(ArgumentError(
                "uncertain nominal values must be finite; got $nominal",
            ))
            isfinite(sigma) || throw(ArgumentError(
                "uncertainty must be finite; got $sigma",
            ))
            sigma >= zero(sigma) || throw(ArgumentError(
                "uncertainty must be nonnegative; got $sigma",
            ))
        end
        return new{T, E}(nominal, sigma)
    end
end

"""
$(TYPEDSIGNATURES)

Return the nominal value stored in an [`UncertainValue`](@ref).
"""
nominal(value::UncertainValue) = value.nominal

"""
$(TYPEDSIGNATURES)

Return the absolute standard uncertainty stored in an
[`UncertainValue`](@ref), in the same unit as its nominal value.
"""
standard_uncertainty(value::UncertainValue) = value.sigma

# A reused Grid instance carries the same automatic key. Separately-created
# equal Grids deliberately do not couple. Named keys provide explicit coupling.
struct AutomaticGridKey
    token::Base.RefValue{Nothing}
end
Base.:(==)(left::AutomaticGridKey, right::AutomaticGridKey) = left.token === right.token
Base.isequal(left::AutomaticGridKey, right::AutomaticGridKey) = left == right
Base.hash(key::AutomaticGridKey, seed::UInt) = hash(objectid(key.token), seed)

struct NamedGridKey{K}
    value::K
end

_grid_key(::Nothing) = AutomaticGridKey(Ref(nothing))
_grid_key(key) = NamedGridKey(key)

"""
$(TYPEDEF)

Represent one axis of deterministic parameter values.

$(TYPEDFIELDS)
"""
struct DeterministicGrid{V <: Tuple, K} <: AbstractGrid
    "Values admitted by the axis."
    vals::V

    "Identity used to couple a selection with another axis."
    key::K
end

"""
$(TYPEDEF)

Represent nominal parameter values and relative standard uncertainties.

$(TYPEDFIELDS)
"""
struct RelativeGrid{V <: Tuple, P <: Tuple, K} <: AbstractUncertainGrid
    "Nominal values admitted by the axis."
    vals::V

    "Relative standard uncertainties expressed as percentages."
    rel_err::P

    "Identity used to couple a selection with another axis."
    key::K

    function RelativeGrid(vals::V, rel_err::P, key::K) where {V <: Tuple, P <: Tuple, K}
        _validate_uncertainty(vals, rel_err, "relative")
        return new{V, P, K}(vals, rel_err, key)
    end
end

"""
$(TYPEDEF)

Represent nominal parameter values and absolute standard uncertainties.

$(TYPEDFIELDS)
"""
struct AbsoluteGrid{V <: Tuple, P <: Tuple, K} <: AbstractUncertainGrid
    "Nominal values admitted by the axis."
    vals::V

    "Absolute standard uncertainties in the same unit as `vals`."
    abs_err::P

    "Identity used to couple a selection with another axis."
    key::K

    function AbsoluteGrid(vals::V, abs_err::P, key::K) where {V <: Tuple, P <: Tuple, K}
        _validate_uncertainty(vals, abs_err, "absolute")
        return new{V, P, K}(vals, abs_err, key)
    end
end

"""
$(TYPEDEF)

Mark values as absolute standard uncertainties for
`Grid(values, AbsoluteError(...))`.

$(TYPEDFIELDS)
"""
struct AbsoluteError{T <: Tuple}
    "Absolute standard uncertainties."
    vals::T

    function AbsoluteError(vals::T) where {T <: Tuple}
        _validate_errors(vals, "absolute")
        return new{T}(vals)
    end
end

_grid_values(value::Tuple) = value
_grid_values(value::AbstractArray) = Tuple(value)
_grid_values(value) = (value,)

function _validate_errors(errors::Tuple, kind::AbstractString)
    isempty(errors) && throw(ArgumentError("$kind uncertainty cannot be empty"))
    for error in errors
        error isa Real ||
            throw(ArgumentError("$kind errors must be real; got $(typeof(error))"))
        isfinite(error) || throw(ArgumentError("$kind errors must be finite; got $error"))
        error >= zero(error) ||
            throw(ArgumentError("$kind errors must be nonnegative; got $error"))
    end
    return nothing
end

function _validate_uncertainty(values::Tuple, errors::Tuple, kind::AbstractString)
    isempty(values) &&
        throw(ArgumentError("$kind uncertainty cannot have an empty nominal axis"))
    _validate_errors(errors, kind)
    for value in values
        value isa Real ||
            throw(ArgumentError("$kind uncertainty requires real nominal values; got $(typeof(value))"))
        isfinite(value) ||
            throw(ArgumentError("$kind nominal values must be finite; got $value"))
    end
    return nothing
end

AbsoluteError(value) = AbsoluteError(_grid_values(value))

"""
$(TYPEDSIGNATURES)

Create the sole representation of parameter variation. Collections are axes
only when passed explicitly to `Grid`; constructors otherwise treat them as
ordinary atomic values. Reusing one Grid instance couples its selections.
`key` couples separately-created Grids deliberately.
"""
function Grid(grid::AbstractGrid; key = nothing)
    key === nothing ? grid :
    throw(ArgumentError("cannot replace the coupling key of an existing Grid"))
end
Grid(value; key = nothing) = DeterministicGrid(_grid_values(value), _grid_key(key))
function Grid(value, error::AbsoluteError; key = nothing)
    AbsoluteGrid(_grid_values(value), error.vals, _grid_key(key))
end
function Grid(value, relative_error; key = nothing)
    RelativeGrid(_grid_values(value), _grid_values(relative_error), _grid_key(key))
end

iterate(grid::DeterministicGrid, state...) = iterate(grid.vals, state...)
length(grid::DeterministicGrid) = length(grid.vals)
size(grid::DeterministicGrid) = (length(grid),)
getindex(grid::DeterministicGrid, index::Integer) = grid.vals[index]
eltype(::Type{<:DeterministicGrid{V}}) where {V} = eltype(V)
Base.IteratorSize(::Type{<:DeterministicGrid}) = Base.HasShape{1}()

function iterate(grid::RelativeGrid, state...)
    item = iterate(Iterators.product(grid.vals, grid.rel_err), state...)
    item === nothing && return nothing
    (value, percent), next_state = item
    sigma = abs(value) * percent / 100
    return UncertainValue(value, sigma), next_state
end

function iterate(grid::AbsoluteGrid, state...)
    item = iterate(Iterators.product(grid.vals, grid.abs_err), state...)
    item === nothing && return nothing
    (value, error), next_state = item
    return UncertainValue(value, error), next_state
end

length(grid::RelativeGrid) = length(grid.vals) * length(grid.rel_err)
length(grid::AbsoluteGrid) = length(grid.vals) * length(grid.abs_err)
size(grid::AbstractUncertainGrid) = (length(grid),)
Base.IteratorSize(::Type{<:AbstractUncertainGrid}) = Base.HasShape{1}()

function getindex(grid::AbstractUncertainGrid, index::Integer)
    1 <= index <= length(grid) || throw(BoundsError(grid, index))
    return first(Iterators.drop(grid, index - 1))
end

extrema(grid::DeterministicGrid) = extrema(grid.vals)

function extrema(grid::RelativeGrid)
    bounds = map(Iterators.product(grid.vals, grid.rel_err)) do (value, percent)
        delta = abs(value) * percent / 100
        (value - delta, value + delta)
    end
    return minimum(first, bounds), maximum(last, bounds)
end

function extrema(grid::AbsoluteGrid)
    bounds = map(Iterators.product(grid.vals, grid.abs_err)) do (value, error)
        (value - error, value + error)
    end
    return minimum(first, bounds), maximum(last, bounds)
end

function _sample_uncertainty(
        rng::Random.AbstractRNG,
        value::UncertainValue{<:Real, <:Real},
        distribution::Symbol
)
    distribution === :normal && return value.nominal + value.sigma * randn(rng)
    distribution === :uniform &&
        return value.nominal + sqrt(3) * value.sigma * (2 * rand(rng) - 1)
    throw(ArgumentError("unsupported distribution :$distribution; expected :normal or :uniform"))
end

function _sample_uncertainty(rng::Random.AbstractRNG, value::UncertainValue{<:Real}, sampler::Function)
    sampler(rng, value.nominal, value.sigma)
end

function _sample_uncertainty(::Random.AbstractRNG, ::UncertainValue, distribution)
    throw(ArgumentError("unsupported distribution $(typeof(distribution)); load its package extension or pass a sampler function"))
end

function rand(rng::Random.AbstractRNG, value::UncertainValue{<:Real}; distribution = :normal)
    iszero(value.sigma) && return float(value.nominal)
    return _sample_uncertainty(rng, value, distribution)
end

rand(value::UncertainValue; kwargs...) = rand(Random.default_rng(), value; kwargs...)

function rand(rng::Random.AbstractRNG, grid::AbstractGrid; distribution = :normal)
    length(grid) == 1 ||
        throw(ArgumentError("rand(Grid) requires one configuration; select a configuration before sampling"))
    value = first(grid)
    return value isa UncertainValue ? rand(rng, value; distribution) : value
end

rand(grid::AbstractGrid; kwargs...) = rand(Random.default_rng(), grid; kwargs...)
