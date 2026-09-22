function assumptions(::Val{:ideal})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

Identify electrostatic image potential coefficients above an ideal conducting
earth plane. For aerial conductors in layer 1, the Maxwell coefficients are

```math
P_{ii}=\\frac{1}{2\\pi\\varepsilon_0}\\ln\\frac{2h_i}{r_i},\\qquad
P_{ij}=\\frac{1}{2\\pi\\varepsilon_0}\\ln\\frac{D_{ij}}{d_{ij}},
```

where ``h_i`` is height, ``r_i`` is the exterior radius, and ``d_{ij}`` and
``D_{ij}`` are the distances to the real conductor and its image, respectively,
all in meters. The coefficients have units of m/F. The external potential
coefficient is zero whenever either conductor is buried, including buried
self and mutual interactions. Insulation potential coefficients remain separate.
For lossless materials, ``Y=j\\omega P^{-1}`` is purely imaginary; capacitance
itself is real.

This scientific identity is registered for backend dispatch. Evaluation by the
owned coaxial backend is not yet implemented and has no numerical fallback.
PSCAD maps this identity to its native potential model. Its direct-integration
setting can produce nonzero aerial conductance even with lossless insulation;
the adapter preserves that native deviation rather than changing these equations.

Reference: PSCAD 5.1 help, *Deriving System Y and Z Matrices*, Eq. (8-25), and
*Mutual Impedance with Earth Return*, Eq. (8-35).
"""
function description(::Type{<:Formula{:ideal}}; compact::Bool = false)
    compact ? "Ideal earth" : "Ideal-earth electrostatic image potential coefficients"
end

function earth_potential_coefficient(
        ::Formula{:ideal}, kind::Union{Val{:self}, Val{:mutual}}, ::Val{1}, ::Val{1},
        functor, pair, workspace)
    throw(ArgumentError("earth_potential_coefficient :ideal ($kind), source layer 1, target layer 1: not yet implemented for the coaxial backend"))
end

function earth_potential_coefficient(
        ::Formula{:ideal}, kind::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace)
    throw(ArgumentError("earth_potential_coefficient :ideal ($kind), source layer 2, target layer 2: not yet implemented for the coaxial backend"))
end

function earth_potential_coefficient(
        ::Formula{:ideal}, kind::Val{:mutual}, ::Val{1}, ::Val{2},
        functor, pair, workspace)
    throw(ArgumentError("earth_potential_coefficient :ideal ($kind), source layer 1, target layer 2: not yet implemented for the coaxial backend"))
end

function earth_potential_coefficient(
        ::Formula{:ideal}, kind::Val{:mutual}, ::Val{2}, ::Val{1},
        functor, pair, workspace)
    throw(ArgumentError("earth_potential_coefficient :ideal ($kind), source layer 2, target layer 1: not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:ideal}, typeof(earth_potential_coefficient)}) =
    FormulationOptions()

:ideal
