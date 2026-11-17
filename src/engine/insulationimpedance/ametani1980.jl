defaults(::Type{Formula{:ametani1980}}) = (;)

function (::Formula{:ametani1980})(
        r_in::T, r_ex::T, μ_r::T, s::Complex{T}) where {T <: Real}
    factor = s * μ₀ * μ_r / 2π
    iszero(r_in) && return zero(factor)
    return factor * log(r_ex / r_in)
end
