defaults(::Type{Formula{:linear}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the resistivity of `material` at `temperature` by the linear law:

```math
\\rho(\\theta) = \\rho_0 \\left[1 + \\alpha (\\theta - T_0)\\right],
```

Here ``\\rho_0`` is material resistivity `rho` \\[Ω·m\\] at the reference temperature
``T_0``, the field `T0` \\[°C\\]. The coefficient ``\\alpha`` is `alpha` \\[1/K\\], and
``\\theta`` is conductor temperature \\[°C\\]. The law returns infinite resistivity
unchanged.

# Arguments

- `Formula{:linear}`: linear temperature law.
- `material`: material of the conductor.
- `temperature`: temperature of the conductor \\[°C\\].

# Returns

- Resistivity \\[Ω·m\\].

# Errors

- Throws `DomainError` when the temperature difference from `T0` is nonfinite or at least
  150 K.
- Throws `DomainError` when the factor in brackets is nonpositive or nonfinite.
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
