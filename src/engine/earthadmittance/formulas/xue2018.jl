function assumptions(::Val{:xue2018})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Generalized underground potential coefficient. The
registered expression is referenced to infinite earth depth; the infinite-depth expression preserves the former default.

**Availability.** Registered scientific identity; the coaxial implementation is
not yet implemented. No numerical fallback is provided.

**Expression.**

```math
P_{e,ij}^{\\infty}=\\frac{j\\omega}{2\\pi(\\sigma_1+j\\omega\\varepsilon_1)}
\\left[K_0(\\gamma_1d_{ij})-K_0(\\gamma_1D_{ij})+2S_{12}^c+
2\\gamma_1^2S_{13}^c\\right],
```

```math
S_{12}^c=\\int_0^\\infty\\frac{e^{-Hu_1}\\lambda^2\\cos(y\\lambda)}
{(\\lambda^2+\\gamma_1^2)[u_0+(\\gamma_0^2/\\gamma_1^2)u_1]}d\\lambda,
\\quad
S_{13}^c=\\int_0^\\infty\\frac{e^{-Hu_1}\\cos(y\\lambda)}
{(\\lambda^2+\\gamma_1^2)(u_0+u_1)}d\\lambda.
```

**Reference.** Haoyan Xue, *General Formulation and Accurate Evaluation of
Earth-Return Parameters for Overhead / Underground Cables*, doctoral thesis,
Polytechnique Montréal, 2018.
[Primary source](https://publications.polymtl.ca/3190/1/2018_HaoyanXue.pdf).
"""
function description(::Type{<:Formula{:xue2018}}; compact::Bool = false)
    compact ? "Xue" : "Xue homogeneous-earth underground potential coefficient (2018) — not yet implemented"
end

function earth_potential_coefficient(
        ::Formula{:xue2018}, kind::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    throw(ArgumentError("earth_potential_coefficient :xue2018 ($kind), source layer 2, target layer 2: not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:xue2018}, typeof(earth_potential_coefficient)}) =
    FormulationOptions()

:xue2018
