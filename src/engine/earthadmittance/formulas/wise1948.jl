function assumptions(::Val{:wise1948})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Wideband homogeneous-earth overhead potential
coefficient.

**Availability.** Registered scientific identity; the coaxial implementation is
not yet implemented. No numerical fallback is provided.

**Expression.**

```math
P_{e,ij}=\\frac{P_{0,ij}+M_{ij}+jN_{ij}}{2\\pi\\varepsilon_0},
```

```math
M_{ij}+jN_{ij}=2\\int_0^\\infty
\\frac{e^{-H\\lambda}\\cos(y_{ij}\\lambda)}
{(\\gamma_1^2/\\gamma_0^2)\\lambda+
\\sqrt{\\lambda^2+\\gamma_1^2-\\gamma_0^2}}d\\lambda,
\\quad P_{0,ij}=\\ln(D_{ij}/d_{ij}).
```

**Reference.** W. H. Wise, “Potential Coefficients for Ground Return
Circuits,” *Bell System Technical Journal*, 27, 365–371, 1948.
"""
function description(::Type{<:Formula{:wise1948}}; compact::Bool = false)
    compact ? "Wise" : "Wise homogeneous-earth overhead potential coefficient (1948) — not yet implemented"
end

function earth_potential_coefficient(
        ::Formula{:wise1948}, kind::Union{Val{:self}, Val{:mutual}}, ::Val{1}, ::Val{1},
        functor, pair, workspace
)
    throw(ArgumentError("earth_potential_coefficient :wise1948 ($kind), source layer 1, target layer 1: not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:wise1948}, typeof(earth_potential_coefficient)}) =
    FormulationOptions()

:wise1948
