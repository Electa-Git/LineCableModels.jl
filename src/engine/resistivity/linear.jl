defaults(::Type{Formula{:linear}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the resistivity of `material` at `temperature` by the linear law:

```math
\\rho(\\theta) = \\rho_0 \\left[1 + \\alpha (\\theta - T_0)\\right],
```

where ``\\rho_0`` is the resistivity `rho` of the material at its reference temperature
``T_0``, the field `T0`. ``\\alpha`` is the coefficient `alpha`, and ``\\theta`` is the
temperature. The law returns an infinite resistivity unchanged.

# Arguments

- `material`: material of the conductor.
- `temperature`: temperature of the conductor \\[°C\\].

# Returns

- Resistivity \\[Ω·m\\].

# Errors

- Throws `DomainError` when `temperature` lies 150 K or more from `T0`, the range of the
  law, or when the factor in brackets is not positive and finite.
"""
function resistivity(
        ::Formula{:linear}, material::Material{T}, temperature::T) where {T <: Real}
    difference = temperature - material.T0
    abs(difference) < 150 || throw(DomainError(temperature,
        "the linear resistivity law holds within 150 K of Material.T0 = $(material.T0) °C"))
    factor = 1 + material.alpha * difference
    isfinite(factor) && factor > 0 || throw(DomainError(factor,
        "the factor 1 + α (θ − T0) of the linear law must be positive and finite"))
    return isinf(material.rho) ? material.rho : material.rho * factor
end
