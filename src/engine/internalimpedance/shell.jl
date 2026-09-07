"Return m*coth(m*t) and m*csch(m*t) without large-exponential overflow."
function _shell_factors(m::Complex{T},thickness::T) where {T<:Real}
    thickness>0 && isfinite(thickness) ||
        throw(DomainError(thickness,"shell thickness must be positive and finite"))
    iszero(m) && return (coth=complex(inv(thickness)),csch=complex(inv(thickness)))
    x=m*thickness
    denominator=-expm1(-2x)
    return (coth=Complex{T}(m*(1+2exp(-2x)/denominator)),
        csch=Complex{T}(2m*exp(-x)/denominator))
end
