function assumptions(::Val{:lucca1994})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Homogeneous-earth mixed-pair model with a corrected
complex-depth approximation for mixed overhead-underground coupling.

**Availability.** Registered scientific identity; the coaxial implementation is
not yet implemented. No numerical fallback is provided.

**Expression.** Its distinctive mixed term is

```math
Z_{e,ij}^{01}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{S}{D}-\\frac23\\left(\\frac{h_e}{S^2}\\right)^3
H(H^2-3y_{ij}^2)\\right],
```

```math
h_e=(j\\omega\\mu_0\\sigma_g)^{-1/2},\\quad
H=h_a+h_g+2h_e,\\quad S=\\sqrt{H^2+y_{ij}^2},\\quad
D=\\sqrt{(h_a+h_g)^2+y_{ij}^2}.
```

**Reference.** G. Lucca, “Mutual Impedance Between an Overhead and a Buried
Line with Earth Return,” *9th International Conference on Electromagnetic
Compatibility*, 1994. DOI: 10.1049/cp:19940679.
"""
function description(::Type{<:Formula{:lucca1994}}; compact::Bool = false)
    compact ? "Lucca" : "Lucca mixed-pair homogeneous-earth impedance (1994)"
end

function earth_impedance(
        ::Formula{:lucca1994}, kind::Val{:mutual}, ::Val{1}, ::Val{2},
        functor, pair, workspace
)
    throw(ArgumentError("earth_impedance :lucca1994 ($kind), source layer 1, target layer 2: not yet implemented for the coaxial backend"))
end

function earth_impedance(
        ::Formula{:lucca1994}, kind::Val{:mutual}, ::Val{2}, ::Val{1},
        functor, pair, workspace
)
    throw(ArgumentError("earth_impedance :lucca1994 ($kind), source layer 2, target layer 1: not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:lucca1994}, typeof(earth_impedance)}) =
    FormulationOptions()

:lucca1994
