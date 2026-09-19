function assumptions(::Val{:pollaczek1926})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Classical homogeneous-earth underground potential coefficient.

**Availability.** Registered scientific identity; the coaxial implementation is
not yet implemented. No numerical fallback is provided.

**Expression.**

```math
P_{e,ij}^{11}=\\frac{j\\omega}{2\\pi(\\sigma_1+j\\omega\\varepsilon_1)}
[K_0(\\gamma_0d_{ij})-K_0(\\gamma_0D_{ij})].
```

**Reference.** F. Pollaczek, “Über das Feld einer unendlich langen
wechselstromdurchflossenen Einfachleitung,” *Elektrische Nachrichtentechnik*,
3, 339–360, 1926; potential-coefficient transcription follows Ametani et al.,
IET, 2021.
"""
function description(::Type{<:Formula{:pollaczek1926}}; compact::Bool = false)
    compact ? "Pollaczek" : "Pollaczek underground potential coefficients (1926) — not yet implemented"
end

function earth_potential_coefficient(
        ::Formula{:pollaczek1926}, kind::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    throw(ArgumentError("earth_potential_coefficient :pollaczek1926 ($kind), source layer 2, target layer 2: not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:pollaczek1926}, typeof(earth_potential_coefficient)}) =
    FormulationOptions()

:pollaczek1926
