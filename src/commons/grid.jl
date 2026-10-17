"""
$(TYPEDEF)

Supertype for the deterministic and uncertainty-bearing finite sources created
by [`Grid`](@ref).
"""
abstract type AbstractGrid end

"""
$(TYPEDEF)

Supertype for finite sources whose points retain an uncertainty descriptor.
"""
abstract type AbstractUncertainGrid <: AbstractGrid end

"""
$(TYPEDEF)

Store a nominal value and its absolute standard uncertainty.

$(TYPEDFIELDS)
"""
struct UncertainValue{T, E}
    "Nominal parameter value."
    nominal::T

    "Absolute standard uncertainty in the same unit as `nominal`."
    sigma::E

    function UncertainValue(nominal::T, sigma::E) where {T, E}
        return validate(new{T, E}(nominal, sigma))
    end
end

# A real uncertain value has a finite nominal value and a finite, nonnegative uncertainty.
function validate(value::UncertainValue{<:Real, <:Real})
    isfinite(value.nominal) || throw(ArgumentError(
        "uncertain nominal values must be finite; got $(value.nominal)",
    ))
    isfinite(value.sigma) || throw(ArgumentError(
        "uncertainty must be finite; got $(value.sigma)",
    ))
    value.sigma >= zero(value.sigma) || throw(ArgumentError(
        "uncertainty must be nonnegative; got $(value.sigma)",
    ))
    return value
end
validate(value::UncertainValue) = value

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
uncertainty(value::UncertainValue) = value.sigma

"""
$(TYPEDEF)

Represent one finite source of deterministic values.

$(TYPEDFIELDS)
"""
struct DeterministicGrid{V <: Tuple} <: AbstractGrid
    "Values admitted by the source."
    vals::V
end

"""
$(TYPEDEF)

Represent the Cartesian product of nominal values and relative standard
uncertainties.

$(TYPEDFIELDS)
"""
struct RelativeGrid{V <: Tuple{Vararg{Real}}, P <: Tuple{Vararg{Real}}} <: AbstractUncertainGrid
    "Nominal values admitted by the source."
    vals::V

    "Relative standard uncertainties expressed as percentages."
    rel_err::P

    function RelativeGrid(vals::V,
            rel_err::P) where {V <: Tuple{Vararg{Real}}, P <: Tuple{Vararg{Real}}}
        return validate(new{V, P}(vals, rel_err))
    end
end
function RelativeGrid(vals::Tuple, rel_err::Tuple)
    throw(ArgumentError("relative uncertainty requires real nominal values and errors; " *
                        "got $(typeof(vals)) and $(typeof(rel_err))"))
end

"""
$(TYPEDEF)

Represent the Cartesian product of nominal values and absolute standard
uncertainties.

$(TYPEDFIELDS)
"""
struct AbsoluteGrid{V <: Tuple{Vararg{Real}}, P <: Tuple{Vararg{Real}}} <: AbstractUncertainGrid
    "Nominal values admitted by the source."
    vals::V

    "Absolute standard uncertainties in the same unit as `vals`."
    abs_err::P

    function AbsoluteGrid(vals::V,
            abs_err::P) where {V <: Tuple{Vararg{Real}}, P <: Tuple{Vararg{Real}}}
        return validate(new{V, P}(vals, abs_err))
    end
end
function AbsoluteGrid(vals::Tuple, abs_err::Tuple)
    throw(ArgumentError("absolute uncertainty requires real nominal values and errors; " *
                        "got $(typeof(vals)) and $(typeof(abs_err))"))
end

"""
$(TYPEDEF)

Mark values as absolute standard uncertainties for
`Grid(values, AbsoluteError(...))`.

$(TYPEDFIELDS)
"""
struct AbsoluteError{T <: Tuple{Vararg{Real}}}
    "Absolute standard uncertainties."
    vals::T

    function AbsoluteError(vals::T) where {T <: Tuple{Vararg{Real}}}
        validate(vals, AbsoluteError)
        return new{T}(vals)
    end
end
function AbsoluteError(vals::Tuple)
    throw(ArgumentError("absolute errors must be real; got $(typeof(vals))"))
end

"""
$(TYPEDSIGNATURES)

Collect the values of one deterministic source. A tuple keeps its values, an array
contributes its elements, and any other value is the one value of the source.
"""
DeterministicGrid(values::AbstractArray) = DeterministicGrid(Tuple(values))
DeterministicGrid(value) = DeterministicGrid((value,))

# Standard uncertainties of the record being built. Relative uncertainties are percentages.
function validate(errors::Tuple,
        ::Type{R}) where {R <: Union{AbsoluteError, AbsoluteGrid, RelativeGrid}}
    kind = R <: RelativeGrid ? "relative" : "absolute"
    isempty(errors) && throw(ArgumentError("$kind uncertainty cannot be empty"))
    for error in errors
        isfinite(error) || throw(ArgumentError("$kind errors must be finite; got $error"))
        error >= zero(error) ||
            throw(ArgumentError("$kind errors must be nonnegative; got $error"))
    end
    return errors
end

function validate(grid::RelativeGrid)
    isempty(grid.vals) &&
        throw(ArgumentError("relative uncertainty cannot have an empty nominal source"))
    validate(grid.rel_err, typeof(grid))
    for value in grid.vals
        isfinite(value) ||
            throw(ArgumentError("relative nominal values must be finite; got $value"))
    end
    return grid
end

function validate(grid::AbsoluteGrid)
    isempty(grid.vals) &&
        throw(ArgumentError("absolute uncertainty cannot have an empty nominal source"))
    validate(grid.abs_err, typeof(grid))
    for value in grid.vals
        isfinite(value) ||
            throw(ArgumentError("absolute nominal values must be finite; got $value"))
    end
    return grid
end

AbsoluteError(value) = AbsoluteError(DeterministicGrid(value).vals)

"""
$(TYPEDSIGNATURES)

Create an explicit finite source. Collections vary only when passed to
`Grid`. Ordinary builder arguments remain atomic domain values.
"""
Grid(grid::AbstractGrid) = grid
Grid(value) = DeterministicGrid(value)
Grid(value, error::AbsoluteError) = AbsoluteGrid(DeterministicGrid(value).vals, error.vals)
function Grid(value, relative_error)
    RelativeGrid(
        DeterministicGrid(value).vals,
        DeterministicGrid(relative_error).vals
    )
end

Base.iterate(grid::DeterministicGrid, state...) = iterate(grid.vals, state...)
Base.length(grid::DeterministicGrid) = length(grid.vals)
Base.size(grid::DeterministicGrid) = (length(grid),)
Base.getindex(grid::DeterministicGrid, index::Integer) = grid.vals[index]
Base.eltype(::Type{<:DeterministicGrid{V}}) where {V} = eltype(V)
#! explicit-imports: off
# Base's iterator trait protocol exposes this value without a public binding.
Base.IteratorSize(::Type{<:DeterministicGrid}) = Base.HasShape{1}()
#! explicit-imports: on

function Base.iterate(grid::RelativeGrid, state...)
    item = iterate(Iterators.product(grid.vals, grid.rel_err), state...)
    item === nothing && return nothing
    (value, percent), next_state = item
    return UncertainValue(value, abs(value) * percent / 100), next_state
end

function Base.iterate(grid::AbsoluteGrid, state...)
    item = iterate(Iterators.product(grid.vals, grid.abs_err), state...)
    item === nothing && return nothing
    (value, error), next_state = item
    return UncertainValue(value, error), next_state
end

Base.length(grid::RelativeGrid) = length(grid.vals) * length(grid.rel_err)
Base.length(grid::AbsoluteGrid) = length(grid.vals) * length(grid.abs_err)
Base.size(grid::AbstractUncertainGrid) = (length(grid),)
#! explicit-imports: off
# Base's iterator trait protocol exposes this value without a public binding.
Base.IteratorSize(::Type{<:AbstractUncertainGrid}) = Base.HasShape{1}()
#! explicit-imports: on

function Base.getindex(grid::RelativeGrid, index::Integer)
    1 <= index <= length(grid) || throw(BoundsError(grid, index))
    value_index = mod1(index, length(grid.vals))
    error_index = cld(index, length(grid.vals))
    value = grid.vals[value_index]
    return UncertainValue(value, abs(value) * grid.rel_err[error_index] / 100)
end

function Base.getindex(grid::AbsoluteGrid, index::Integer)
    1 <= index <= length(grid) || throw(BoundsError(grid, index))
    value_index = mod1(index, length(grid.vals))
    error_index = cld(index, length(grid.vals))
    return UncertainValue(grid.vals[value_index], grid.abs_err[error_index])
end

Base.extrema(grid::DeterministicGrid) = extrema(grid.vals)

function Base.extrema(grid::RelativeGrid)
    bounds = map(Iterators.product(grid.vals, grid.rel_err)) do (value, percent)
        delta = abs(value) * percent / 100
        (value - delta, value + delta)
    end
    return minimum(first, bounds), maximum(last, bounds)
end

function Base.extrema(grid::AbsoluteGrid)
    bounds = map(Iterators.product(grid.vals, grid.abs_err)) do (value, error)
        (value - error, value + error)
    end
    return minimum(first, bounds), maximum(last, bounds)
end

"""
$(TYPEDSIGNATURES)

Draw one scalar realization from an uncertainty descriptor using a caller-owned
random number generator and the selected sampling family.
"""
function sample_uncertainty(
        rng::Random.AbstractRNG,
        value::UncertainValue{<:Real, <:Real},
        distribution::Symbol
)
    distribution === :normal && return value.nominal + value.sigma * randn(rng)
    distribution === :uniform &&
        return value.nominal + sqrt(3) * value.sigma * (2 * rand(rng) - 1)
    throw(ArgumentError(
        "unsupported distribution :$distribution; expected :normal or :uniform",
    ))
end

function sample_uncertainty(
        rng::Random.AbstractRNG,
        value::UncertainValue{<:Real},
        sampler::Function
)
    sampler(rng, value.nominal, value.sigma)
end

function sample_uncertainty(::Random.AbstractRNG, ::UncertainValue, distribution)
    throw(ArgumentError(
        "unsupported distribution $(typeof(distribution)); load its package extension or pass a sampler function",
    ))
end

function Base.rand(
        rng::Random.AbstractRNG,
        value::UncertainValue{<:Real};
        distribution = :normal
)
    iszero(value.sigma) && return float(value.nominal)
    return sample_uncertainty(rng, value, distribution)
end

function Base.rand(value::UncertainValue{<:Real}; distribution = :normal)
    return rand(Random.default_rng(), value; distribution)
end

function Base.rand(
        rng::Random.AbstractRNG,
        grid::AbstractGrid;
        distribution = :normal
)
    isempty(grid) && throw(ArgumentError("cannot sample an empty Grid"))
    return realize(rng, grid[rand(rng, 1:length(grid))], distribution)
end

function Base.rand(grid::AbstractGrid; distribution = :normal)
    return rand(Random.default_rng(), grid; distribution)
end
