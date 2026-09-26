function assumptions(::Val{:wedepohl1973})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Low-frequency underground expansion for conductive,
nonmagnetic earth.

**Availability.** Registered scientific identity; the coaxial implementation is
not yet implemented. No numerical fallback is provided.

**Expression.**

```math
Z_{e,ii}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[-\\ln
\\left(\\frac{e_c\\gamma_1r_i}{2}\\right)+\\frac12-
\\frac43\\gamma_1h_i\\right],
```

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[-\\ln
\\left(\\frac{e_c\\gamma_1d_{ij}}{2}\\right)+\\frac12-
\\frac23\\gamma_1H\\right],\\qquad e_c=1.7811.
```

Here ``d_{ij}=\\sqrt{x_{ij}^2+(h_i-h_j)^2}`` is the distance between
cable axes and ``H=h_i+h_j`` is the sum of burial depths, all in meters.

**Reference.** L. M. Wedepohl and D. J. Wilcox, “Transient Analysis of
Underground Power-Transmission Systems: System-Model and Wave-Propagation
Characteristics,” *Proceedings of the IEE*, 120, 253–260, 1973, Eqs. (7)–(8).
DOI: 10.1049/piee.1973.0056.
"""
function description(::Type{<:Formula{:wedepohl1973}}; compact::Bool = false)
    compact ? "Wedepohl" : "Wedepohl-Wilcox low-frequency underground approximation (1973)"
end

function earth_impedance(
        ::Formula{:wedepohl1973}, kind::Val{:self}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    throw(ArgumentError("earth_impedance :wedepohl1973 ($kind), source layer 2, target layer 2: not yet implemented for the coaxial backend"))
end

function earth_impedance(
        ::Formula{:wedepohl1973}, kind::Val{:mutual}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    throw(ArgumentError("earth_impedance :wedepohl1973 ($kind), source layer 2, target layer 2: not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:wedepohl1973}, typeof(earth_impedance)}) =
    FormulationOptions()

:wedepohl1973
