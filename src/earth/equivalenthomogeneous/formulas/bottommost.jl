
"""
$(TYPEDSIGNATURES)

**Identification.** Homogeneous earth represented by the deepest soil layer.

**Expression.** Select the deepest layer's resistivity \\[Ω·m\\], relative
permittivity, and relative permeability (both dimensionless). Property
vectors use index 1 for air and indices 2 through N for soil; this default
selects index N for every conductor-pair layout.

**Scope.** Martins-Britto et al. found that deep-layer conductivity predominated
in the magnetic ground-return impedance of the multilayer soil cases they
studied. This supports using deep-layer resistivity as a default, subject to
the soil structure and frequency range. Large conductivity contrasts and
high-frequency effects limit that approximation. Selecting one layer does
not implement their equivalent-conductivity formula, and their result does
not establish the accuracy of selecting its permittivity or permeability.

Select a different rule to change the equivalent material; see [`Formula`](@ref).

**Reference.** A. G. Martins-Britto, F. V. Lopes, and S. R. M. J. Rondineau,
“Multilayer Earth Structure Approximation by a Homogeneous Conductivity Soil
for Ground Return Impedance Calculations,” *IEEE Transactions on Power
Delivery*, 35(2), 881–891, 2020.
[DOI: 10.1109/TPWRD.2019.2930406](https://doi.org/10.1109/TPWRD.2019.2930406).
"""
description(::Type{<:Formula{:bottommost}}; compact::Bool=false) = compact ? "Bottommost" : "Bottommost earth layer"

function equivalent_material(
        ::Formula{:bottommost}, ::Union{Val{:self}, Val{:mutual}}, ::Val{S}, ::Val{T},
        rho, eps_r, mu_r, model, pair, frequency,
        values, options, workspace
) where {S, T}
    S >= 1 && T >= 1 || throw(ArgumentError("physical layer indices must be positive"))
    return EarthMaterial(rho[end], eps_r[end], mu_r[end])
end

computation_options(::FormulaMethod{<:Formula{:bottommost}, typeof(equivalent_material)}) = (;)

:bottommost
