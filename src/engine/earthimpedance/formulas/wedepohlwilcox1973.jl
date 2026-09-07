function routes(identifier::Val{:WedepohlWilcox1973})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:WedepohlWilcox1973})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:WedepohlWilcox1973}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Low-frequency underground expansion for conductive,
nonmagnetic earth.

**Expression.**

```math
Z_{e,ii}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[-\\ln
\\left(\\frac{e_c\\gamma_1r_i}{2}\\right)+\\frac12-
\\frac43\\gamma_1h_i\\right],
```

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[-\\ln
\\left(\\frac{e_c\\gamma_1y_{ij}}{2}\\right)+\\frac12-
\\frac23\\gamma_1H\\right],\\qquad e_c=1.7811.
```

**Reference.** L. M. Wedepohl and D. J. Wilcox, “Transient Analysis of
Underground Power-Transmission Systems: System-Model and Wave-Propagation
Characteristics,” *Proceedings of the IEE*, 120, 253–260, 1973.
"""
function description(::Formula{:WedepohlWilcox1973})
    "Wedepohl-Wilcox low-frequency underground approximation (1973)"
end

function propagation_constant(
        ::Val{:WedepohlWilcox1973}, jω, permeability, permittivity
)
    (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:WedepohlWilcox1973})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:WedepohlWilcox1973), formula,
        rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Wedepohl-Wilcox low-frequency underground terms:

```math
Z_{e,ii}=\frac{j\omega\mu_0}{2\pi}
\left[-\ln\left(\frac{e_c\gamma_1r_i}{2}\right)+\frac12
-\frac43\gamma_1h_i\right],
```

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}
\left[-\ln\left(\frac{e_c\gamma_1y_{ij}}{2}\right)+\frac12
-\frac23\gamma_1(h_i+h_j)\right],\qquad e_c=1.7811.
```
"""
function earth_impedance(
        ::Val{:WedepohlWilcox1973}, ::Val{:self}, functor, pair
)
    validate(pair, FormulaMethod(Val(:WedepohlWilcox1973), earth_impedance, Val(:self)), functor)
    state = functor.state
    geometry = _geometry(pair)
    e_c = oftype(geometry.h_i, 1.7811)
    bracket = -log(e_c * state.gamma[2] * geometry.y_ij / 2) +
              one(geometry.h_i) / 2 -
              (4 * one(geometry.h_i) / 3) * state.gamma[2] * geometry.h_i
    return state.jω * state.mu[1] / (2π) * bracket
end

function earth_impedance(
        ::Val{:WedepohlWilcox1973}, ::Val{:mutual}, functor, pair
)
    validate(pair, FormulaMethod(Val(:WedepohlWilcox1973), earth_impedance, Val(:mutual)), functor)
    state = functor.state
    geometry = _geometry(pair)
    e_c = oftype(geometry.H, 1.7811)
    bracket = -log(e_c * state.gamma[2] * geometry.y_ij / 2) +
              one(geometry.H) / 2 -
              (2 * one(geometry.H) / 3) * state.gamma[2] * geometry.H
    return state.jω * state.mu[1] / (2π) * bracket
end

function validate(
        pair::EarthPair, route::FormulaMethod{:WedepohlWilcox1973, typeof(earth_impedance)}, formula
)
    validate(pair)
    (pair.layers[1] > 1 && pair.layers[2] > 1) || throw(ArgumentError(
        ":WedepohlWilcox1973 earth impedance requires underground conductors; pair ($(pair.row), $(pair.column)) has layers $(pair.layers)"))
    route.arguments != (Val(:self),) && iszero(pair.separation) && throw(DomainError(
        pair.separation, ":WedepohlWilcox1973 mutual closed form requires nonzero horizontal cable separation"))
    return pair
end

:WedepohlWilcox1973
