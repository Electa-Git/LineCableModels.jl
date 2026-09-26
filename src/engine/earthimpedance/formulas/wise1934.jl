function assumptions(::Val{:wise1934})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Homogeneous-earth wideband overhead integral retaining
earth displacement current and magnetic permeability.

**Availability.** Registered scientific identity; the coaxial implementation is
not yet implemented. No numerical fallback is provided.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+2\\int_0^\\infty
\\frac{\\mu_1e^{-\\lambda H}}
{\\lambda\\mu_1+a_1\\mu_0}\\cos(y_{ij}\\lambda)d\\lambda\\right],
\\quad a_1=\\sqrt{\\lambda^2+\\gamma_1^2-\\gamma_0^2}.
```

**Reference.** W. H. Wise, “Propagation of High-Frequency Currents in Ground
Return Circuits,” *Proceedings of the Institute of Radio Engineers*, 22,
522–527, 1934.
"""
function description(::Type{<:Formula{:wise1934}}; compact::Bool = false)
    compact ? "Wise" : "Wise homogeneous-earth overhead impedance (1934) — not yet implemented"
end

function earth_impedance(
        ::Formula{:wise1934}, kind::Union{Val{:self}, Val{:mutual}}, ::Val{1}, ::Val{1},
        functor, pair, workspace
)
    throw(ArgumentError("earth_impedance :wise1934 ($kind), source layer 1, target layer 1: not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:wise1934}, typeof(earth_impedance)}) =
    FormulationOptions()

:wise1934
