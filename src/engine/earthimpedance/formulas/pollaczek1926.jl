function assumptions(::Val{:pollaczek1926})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Classical homogeneous-earth underground integral.

**Availability.** Registered scientific identity; the coaxial implementation is
not yet implemented. No numerical fallback is provided.

**Expression.** The underground term is

```math
Z_{e,ij}^{11}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
K_0(\\gamma_1d_{ij})-K_0(\\gamma_1D_{ij})+2\\int_0^\\infty
\\frac{e^{-H\\sqrt{\\lambda^2+\\gamma_1^2}}}
{\\lambda+\\sqrt{\\lambda^2+\\gamma_1^2}}
\\cos(y_{ij}\\lambda)d\\lambda\\right],
```

**Reference.** F. Pollaczek, “Über das Feld einer unendlich langen
wechselstromdurchflossenen Einfachleitung,” *Elektrische Nachrichtentechnik*,
3, 339–360, 1926.
"""
function description(::Type{<:Formula{:pollaczek1926}}; compact::Bool = false)
    compact ? "Pollaczek" : "Pollaczek homogeneous-earth underground impedance (1926)"
end

function earth_impedance(
        ::Formula{:pollaczek1926}, kind::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    throw(ArgumentError("earth_impedance :pollaczek1926 ($kind), source layer 2, target layer 2: not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:pollaczek1926}, typeof(earth_impedance)}) =
    FormulationOptions()

:pollaczek1926
