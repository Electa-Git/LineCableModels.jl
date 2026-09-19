function assumptions(::Val{:carson1926})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Classical homogeneous, conductive-earth overhead
impedance. Displacement currents and longitudinal propagation are neglected.

**Availability.** Registered scientific identity; the coaxial implementation is
not yet implemented. No numerical fallback is provided.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+2\\int_0^\\infty
\\frac{e^{-H\\lambda}\\cos(y_{ij}\\lambda)}
{\\lambda+\\sqrt{\\lambda^2+\\gamma_g^2}}d\\lambda\\right],
\\qquad \\gamma_g^2=j\\omega\\mu_0\\sigma_g.
```

**Reference.** J. R. Carson, “Wave Propagation in Overhead Wires with Ground
Return,” *Bell System Technical Journal*, 5, 539–554, 1926.
"""
function description(::Type{<:Formula{:carson1926}}; compact::Bool = false)
    compact ? "Carson" : "Carson homogeneous-earth overhead impedance (1926) — not yet implemented"
end

function earth_impedance(
        ::Formula{:carson1926}, kind::Union{Val{:self}, Val{:mutual}}, ::Val{1}, ::Val{1},
        functor, pair, workspace
)
    throw(ArgumentError("earth_impedance :carson1926 ($kind), source layer 1, target layer 1: not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:carson1926}, typeof(earth_impedance)}) =
    FormulationOptions()

:carson1926
