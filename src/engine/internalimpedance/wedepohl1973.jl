defaults(::Type{Formula{:wedepohl1973}}) = (;)
domain(::Formula{:wedepohl1973}) = (wall = 0.0 .. 0.125, order = 0 .. 0)

function (::Formula{:wedepohl1973})(
        r_in::T, r_ex::T, ρ::T, μ_r::T, s::Complex{T}) where {T <: Real}
    μ = μ₀ * μ_r
    σ = inv(ρ)
    m = sqrt(s * μ * σ)
    if iszero(r_in)
        outer = (ρ * m / (2π * r_ex)) * coth(0.777 * m * r_ex) +
                0.356 * ρ / (π * r_ex^2)
        return (inner = zero(outer), outer, transfer = zero(outer))
    end
    Δ = r_ex - r_in
    inner = (ρ * m / (2π * r_in)) * coth(m * Δ) -
            ρ / (2π * r_in * (r_in + r_ex))
    outer = (ρ * m / (2π * r_ex)) * coth(m * Δ) +
            ρ / (2π * r_ex * (r_in + r_ex))
    transfer = (ρ * m / (π * (r_in + r_ex))) * csch(m * Δ)
    return (; inner, outer, transfer)
end
