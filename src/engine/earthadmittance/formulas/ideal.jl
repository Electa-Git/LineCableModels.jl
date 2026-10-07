function description(::Type{<:Formula{:ideal}}; compact::Bool = false)
    compact ? "Ideal earth" : "Ideal-earth electrostatic image potential coefficients"
end

"""
$(TYPEDSIGNATURES)

Evaluate Maxwell's electrostatic image potential coefficients above an
ideal conducting earth plane, in m/F.

# Assumptions

Aerial conductors in layer 1 above an ideal conducting plane, with air
permittivity ``\\varepsilon_0`` [F/m]. The external potential coefficient is
zero whenever either conductor is buried, including buried self and mutual
interactions. Insulation potential coefficients remain separate.

# Expression

```math
P_{ii}=\\frac{1}{2\\pi\\varepsilon_0}\\ln\\frac{2h_i}{r_i},\\qquad
P_{ij}=\\frac{1}{2\\pi\\varepsilon_0}\\ln\\frac{D_{ij}}{d_{ij}},
```

```math
d_{ij}=\\sqrt{x_{ij}^2+(h_i-h_j)^2},\\qquad
D_{ij}=\\sqrt{x_{ij}^2+(h_i+h_j)^2}.
```

Here ``h_i`` is aerial height, ``r_i`` is the exterior radius, and ``x_{ij}``
is horizontal separation. ``d_{ij}`` and ``D_{ij}`` are distances to the real
conductor and its image, respectively, all in meters.
For lossless materials, the complete potential matrix gives purely imaginary
admittance ``Y=j\\omega P^{-1}``. Capacitance itself is real.

The coaxial backend evaluates these equations directly. PSCAD maps this
selection to its native potential model. Its direct-integration setting can
produce nonzero aerial conductance even with lossless insulation. The adapter
preserves that native deviation rather than changing these equations.

# Reference

PSCAD 5.1 help, *Deriving System Y and Z Matrices*, Eq. (8-25), and
*Mutual Impedance with Earth Return*, Eq. (8-35).
"""
function earth_potential_coefficient(
        ::Formula{:ideal}, ::Val{:self}, ::Val{1}, ::Val{1},
        functor, pair, workspace)
    ε0 = vacuum_permittivity(typeof(real(functor.state.jω)))
    return log(2 * pair.heights[1] / pair.radius) / (2 * (one(ε0) * π) * ε0)
end

function earth_potential_coefficient(
        ::Formula{:ideal}, ::Val{:mutual}, ::Val{1}, ::Val{1},
        functor, pair, workspace)
    ε0 = vacuum_permittivity(typeof(real(functor.state.jω)))
    hi, hj = pair.heights
    D = hypot(pair.separation, hi + hj)
    d = hypot(pair.separation, hi - hj)
    return log(D / d) / (2 * (one(ε0) * π) * ε0)
end

function earth_potential_coefficient(
        ::Formula{:ideal}, ::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace)
    return zero(functor.state.jω)
end

function earth_potential_coefficient(
        ::Formula{:ideal}, ::Val{:mutual}, ::Val{1}, ::Val{2},
        functor, pair, workspace)
    return zero(functor.state.jω)
end

function earth_potential_coefficient(
        ::Formula{:ideal}, ::Val{:mutual}, ::Val{2}, ::Val{1},
        functor, pair, workspace)
    return zero(functor.state.jω)
end

function formulation_options(::Expression{
        <:Formula{:ideal}, typeof(earth_potential_coefficient)})
    FormulationOptions()
end

:ideal
