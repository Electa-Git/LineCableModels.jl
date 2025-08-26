# Make Material fieldwise-coercible with specialized dispatch
import ..Utils: coerce_to_T

coerce_to_T(m::Material, ::Type{T}) where {T<:REALSCALAR} = Material(
    coerce_to_T(m.rho, T),
    coerce_to_T(m.eps_r, T),
    coerce_to_T(m.mu_r, T),
    coerce_to_T(m.T0, T),
    coerce_to_T(m.alpha, T),
)