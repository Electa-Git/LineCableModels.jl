function assumptions(::Val{:saad1996})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Closed-form underground approximation combining the
direct cylindrical term and an interface correction.

**Availability.** Registered scientific identity; the coaxial implementation is
not yet implemented. No numerical fallback is provided.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
K_0(\\gamma_1R_{ab})+
\\frac{2e^{-H\\gamma_1}}{4+\\gamma_1^2R_{ab}^2}\\right].
```

**Reference.** O. Saad, G. Gaba, and M. Giroux, “A Closed-Form Approximation
for Ground Return Impedance of Underground Cables,” *IEEE Transactions on
Power Delivery*, 11(3), 1536–1545, 1996.
"""
function description(::Type{<:Formula{:saad1996}}; compact::Bool = false)
    compact ? "Saad" : "Saad underground closed form (1996) — not yet implemented"
end

function earth_impedance(
        ::Formula{:saad1996}, kind::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    throw(ArgumentError("earth_impedance :saad1996 ($kind), source layer 2, target layer 2: not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:saad1996}, typeof(earth_impedance)}) =
    FormulationOptions()

:saad1996
