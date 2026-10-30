defaults(::Type{Formula{:linear}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the resistivity of `material` at the temperature `θ` by the formula `Formula{:linear}`,
whose equation is in the section "Formulas" of [`Resistivity`](@ref).

# Arguments

- `Formula{:linear}`: linear temperature law.
- `material`: material of the conductor.
- `θ`: temperature of the conductor \\[°C\\].

# Returns

- Resistivity \\[Ω·m\\].

# Errors

- Throws `DomainError` when the temperature difference from `T0` is nonfinite or at least
  150 K.
- Throws `DomainError` when the factor in brackets is nonpositive or nonfinite.
"""
function (::Formula{:linear})(material::Material{T}, θ::T) where {T <: Real}
    difference = θ - material.T0
    abs(difference) < 150 || throw(DomainError(θ,
        "the linear resistivity law holds within 150 K of Material.T0 = $(material.T0) °C"))
    factor = 1 + material.alpha * difference
    isfinite(factor) && factor > 0 || throw(DomainError(factor,
        "the factor 1 + α (θ − T0) of the linear law must be positive and finite"))
    return isinf(material.rho) ? material.rho : material.rho * factor
end
