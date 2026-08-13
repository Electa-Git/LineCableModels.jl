# Identity: no allocation if already at T
@inline coerce_to_T(m::Material{T}, ::Type{T}) where {T} = m

# Cross-T rebuild: use the TYPED constructor to avoid surprise promotion
@inline coerce_to_T(m::Material{S}, ::Type{T}) where {S,T} = Material{T}(
    coerce_to_T(m.rho, T),
    coerce_to_T(m.eps_r, T),
    coerce_to_T(m.mu_r, T),
    coerce_to_T(m.T0, T),
    coerce_to_T(m.alpha, T),
)