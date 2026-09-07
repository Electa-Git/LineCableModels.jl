function routes(identifier::Val{:Wise1934})
    return (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

assumptions(::Val{:Wise1934}) = (
    air = _lossless,
    earth = _full,
    permeability = vacuum_permeability
)

propagation(::Val{:Wise1934}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Straight infinitely long thin wires; the external formula uses direct/image distances ``\\rho'`` and ``\\rho''`` and no insulation layer. |
| Calculated quantities | P.u.l. mutual series impedance of two parallel overhead ground-return wires, with earth polarization current retained |
| Earth structure | Plane homogeneous earth half-space below air. |
| Model and approximation | The earth-polarization integral is not approximated after the source fixes ``\\gamma=jk`` and unit permeability. Wise then states that Carson's asymptotic and convergent series are reused by replacing Carson's ``r`` by ``r\\sqrt{1+j(\\varepsilon-1)/(2c\\lambda\\sigma)}``; those evaluator series are not reproduced as separate physical formulations here. |
| Main source | W. H. Wise (1934), modifying Carson's homogeneous-earth formula to retain earth dielectric polarization |
| Citation key(s) | `:Wise1934` |
| Evidence status | Original publication page images checked |

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+2\\int_0^\\infty
\\frac{\\mu_1e^{-\\lambda H}}
{\\lambda\\mu_1+a_1\\mu_0}\\cos(y_{ij}\\lambda)d\\lambda\\right],
\\quad a_1=\\sqrt{\\lambda^2+\\gamma_1^2-\\gamma_0^2}.
```

**Reference.** W. H. Wise, “Propagation of High-Frequency Currents in Ground
Return Circuits,” *Proceedings of the Institute of Radio Engineers*, 22,
522–527, 1934.
"""
description(::Formula{:Wise1934}) = "Wise homogeneous-earth overhead impedance (1934)"

function propagation_constant(::Val{:Wise1934}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Wise1934})(rho, epsilon, mu, jω, Γ, segments = nothing)
    isinf(first(rho)) || throw(ArgumentError(":Wise1934 requires lossless air"))
    return _homogeneous_functor(Val(:Wise1934), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate Wise's homogeneous-earth overhead impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[\ln\frac{D_{ij}}{d_{ij}}+
2\int_0^\infty F_{ij}^{W}(\lambda)\cos(y_{ij}\lambda)\,d\lambda\right],
```

```math
F_{ij}^{W}=\frac{\mu_1e^{-\lambda(h_i+h_j)}}
{\lambda\mu_1+a_1\mu_0},\qquad
a_1=\sqrt{\lambda^2+\gamma_1^2-\gamma_0^2}.
```
"""
function earth_impedance(
        ::Val{:Wise1934}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    contrast = state.gamma_medium_squared[2] - state.gamma_medium_squared[1]
    self=pair.row==pair.column
    lateral=self ? zero(geometry.H) : geometry.y_ij
    ideal=self ? log(geometry.H/geometry.y_ij) : log(geometry.D_ij/geometry.d_ij)
    integral = _height_quadrature(state,geometry.H) do lambda
        a_1 = spectral_root(lambda^2 + contrast,state.jω)
        state.mu[2] * exp(-lambda * geometry.H) * cos(lambda * lateral) /
        (lambda * state.mu[2] + a_1 * state.mu[1])
    end
    return state.jω * state.mu[1] / (2*(one(geometry.H)*π)) *
           (ideal + 2 * integral)
end

:Wise1934
