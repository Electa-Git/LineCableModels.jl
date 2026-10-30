defaults(::Type{Formula{:ametani1980}}) = (;)

function (::Formula{:ametani1980})(
        r_in::T, r_ex::T, μ_r::T, s::Complex{T}) where {T <: Real}
    factor = s * vacuum_permeability * μ_r / 2π
    if iszero(r_in) || r_in == r_ex
        return zero(factor)
    end
    return factor * log(r_ex / r_in)
end
