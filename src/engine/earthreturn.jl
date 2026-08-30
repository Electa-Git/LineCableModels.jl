@inline function earth!(
        destination::AbstractMatrix{Complex{T}},
        pairs::AbstractVector{<:EarthPair{T}},
        input::NamedTuple,
        earth,
        frequency::Int,
        formula,
        Γ,
        segments,
        thickness = nothing
) where {T <: Real}
    s = input.jω[frequency]
    @inbounds for index in eachindex(pairs)
        pair = pairs[index]
        rho = @view earth.rho[:, index]
        epsilon = @view earth.epsilon[:, index]
        mu = @view earth.mu[:, index]
        functor = formula(rho, epsilon, mu, s, Γ, segments, thickness)
        kind = pair.row == pair.column ? Val(:self) : Val(:mutual)
        value = functor(kind, pair)
        destination[pair.row, pair.column] = value
        destination[pair.column, pair.row] = value
    end
    return destination
end

@inline function earth!(
        destination::AbstractMatrix{Complex{T}},
        pairs::AbstractVector{<:EarthPair{T}},
        input::NamedTuple,
        earth,
        frequency::Int,
        formula,
        Γ,
        segments,
        thickness::AbstractVector
) where {T <: Real}
    s = input.jω[frequency]
    rho = @view earth.rho[:, firstindex(pairs)]
    epsilon = @view earth.epsilon[:, firstindex(pairs)]
    mu = @view earth.mu[:, firstindex(pairs)]
    functor = formula(rho, epsilon, mu, s, Γ, segments, thickness)
    @inbounds for pair in pairs
        kind = pair.row == pair.column ? Val(:self) : Val(:mutual)
        value = functor(kind, pair)
        destination[pair.row, pair.column] = value
        destination[pair.column, pair.row] = value
    end
    return destination
end

@inline _gamma(::Nothing, frequency::Int) = nothing
@inline _gamma(values::AbstractVector, frequency::Int) = values[frequency]
