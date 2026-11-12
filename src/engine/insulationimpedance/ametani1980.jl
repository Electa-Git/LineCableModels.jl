defaults(::Type{Formula{:ametani1980}}) = (;)

function (::Formula{:ametani1980})(
        r_in::T, r_ex::T, μ_r::T, s::Complex{T}) where {T <: Real}
    factor = s * μ₀ * μ_r / 2π
    if iszero(r_in) || r_in == r_ex
        return zero(factor)
    end
    return factor * log(r_ex / r_in)
end

function (::Formula{:ametani1980})(
        n::Int, r_in::T, r_ex::T, μ_r::T, s::Complex{T}) where {T <: Real}
    z = s * μ₀ * μ_r / (2π * n)
    ρ = (r_in / r_ex)^(2n)
    inner = z * (1 + ρ) / (1 - ρ)
    return (; inner, outer = inner, transfer = 2z * (r_in / r_ex)^n / (1 - ρ))
end
