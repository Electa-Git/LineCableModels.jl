function routes(identifier::Val{:Wise1948})
    return (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Wise1948})
    (
        air = _full,
        earth = _full,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Wise1948}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Wideband homogeneous-earth overhead potential
coefficient.

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
function description(::Formula{:Wise1948})
    "Wise wideband homogeneous-earth overhead potential coefficient (1948)"
end

function propagation_constant(
        ::Val{:Wise1948}, jω, permeability, permittivity
)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Wise1948})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(Val(:Wise1948), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate Wise's wideband overhead earth potential coefficient:

```math
P_{e,ij}=\frac{P_{0,ij}+M_{ij}+jN_{ij}}{2\pi\varepsilon_0},
```

```math
M_{ij}+jN_{ij}=2\int_0^\infty
\frac{e^{-(h_i+h_j)\lambda}\cos(y_{ij}\lambda)}
{(\gamma_1^2/\gamma_0^2)\lambda+
\sqrt{\lambda^2+\gamma_1^2-\gamma_0^2}}d\lambda,
\qquad P_{0,ij}=\ln(D_{ij}/d_{ij}).
```
"""
function earth_potential_coefficient(
        ::Val{:Wise1948}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    ratio = gamma_1_squared / gamma_0_squared
    integral = _quadrature(state) do lambda
        radial = sqrt(lambda^2 + gamma_1_squared - gamma_0_squared)
        exp(-geometry.H * lambda) * cos(geometry.y_ij * lambda) /
        (ratio * lambda + radial)
    end
    return (log(geometry.D_ij / geometry.d_ij) + 2 * integral) /
           (2π * state.epsilon[1])
end

:Wise1948
