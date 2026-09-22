function assumptions(::Val{:gary1976})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Complex-depth logarithmic approximation for overhead
conductors.

**Availability.** Registered scientific identity; the coaxial implementation is
not yet implemented. No numerical fallback is provided.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\ln\\frac{S_{ij}}{d_{ij}},\\qquad
S_{ij}=\\sqrt{(H+2h_e)^2+y_{ij}^2},\\qquad
h_e=(j\\omega\\mu_0\\sigma_1)^{-1/2}.
```

**Reference.** C. Gary, “Approche complète de la propagation multifilaire en
haute fréquence par utilisation des matrices complexes,” *EDF Bulletin de la
Direction des Études et Recherches*, série B, 1976; formula as reproduced in
Ametani et al., IET, 2021.
"""
function description(::Type{<:Formula{:gary1976}}; compact::Bool = false)
    compact ? "Gary" : "Gary complex-depth approximation (1976)"
end

function earth_impedance(
        ::Formula{:gary1976}, kind::Union{Val{:self}, Val{:mutual}}, ::Val{1}, ::Val{1},
        functor, pair, workspace
)
    throw(ArgumentError("earth_impedance :gary1976 ($kind), source layer 1, target layer 1: not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:gary1976}, typeof(earth_impedance)}) =
    FormulationOptions()

:gary1976
