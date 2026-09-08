
"""
$(TYPEDSIGNATURES)

**Identification.** Bottommost-layer equivalent homogeneous earth.

**Expression.** Represent the earth by the bottommost soil layer's resistivity, relative
permittivity, and relative permeability.

The property vectors use index 1 for air and indices 2 through N for soil.
This package default selects index N for every conductor-pair layout.

**Reference.** Package layer-selection policy with no literature approximation.
"""
description(::Formula{:default}) = "Bottommost earth layer"

function equivalent_material(
        ::Val{:default}, ::Union{Val{:self}, Val{:mutual}}, ::Val{S}, ::Val{T},
        rho, eps_r, mu_r, model, pair, frequency,
        values, options, workspace
) where {S, T}
    S >= 1 && T >= 1 || throw(ArgumentError("physical layer indices must be positive"))
    return EarthMaterial(rho[end], eps_r[end], mu_r[end])
end

computation_options(::FormulaMethod{:default, typeof(equivalent_material)}) = (;)

:default
