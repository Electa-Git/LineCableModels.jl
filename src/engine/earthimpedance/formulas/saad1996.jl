function routes(identifier::Val{:Saad1996})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Saad1996})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Saad1996}) = Val(:zero)

function Formula(::Val{:Saad1996}; approximation::Symbol = :closed_form, kwargs...)
    approximation in (:closed_form, :small_argument) || throw(ArgumentError(
        ":Saad1996 approximation must be :closed_form or :small_argument"
    ))
    identifier = Val(:Saad1996)
    defaults = approximation === :closed_form ? routes(identifier) : (
        self = FormulaMethod(identifier, earth_impedance, Val(:small_argument)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:small_argument)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
    overrides = (; kwargs...)
    isempty(setdiff(keys(overrides), keys(defaults))) || throw(ArgumentError(
        "unknown routes for earth-impedance formula :Saad1996"
    ))
    return Formula(identifier, merge(defaults, overrides), assumptions(identifier))
end
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Cables are treated as filamentary for the earth-return development; the self term substitutes the cable radius ``R``. Mutual geometry uses depths ``h_i,h_j`` and horizontal separation ``x``. |
| Calculated quantities | Per-unit-length self and mutual earth-return impedances of buried horizontal cables |
| Earth structure | Homogeneous semi-infinite earth below air. |
| Model and approximation | Starting from the Pollaczek/Wedepohl expressions (1)–(4), the authors deform the integration path, approximate ``\\sqrt{\\delta^2+1}/(\\delta+\\sqrt{\\delta^2+1})`` by ``(1+e^{-2\\delta})/2`` in (15), and then use ``\\sqrt{\\delta^2+1}\\simeq1`` in the rapidly decaying part, equations (20)–(21). The paper reports about 3% maximum relative error for the first kernel approximation and restricts the contour proof to ``x/\\ell<1``. A further small-argument reduction is kept in a separate record. |
| Main source | O. Saad, G. Gaba, and M. Giroux (1996) |
| Citation key(s) | `:Saad1996` |
| Evidence status | Original PDF page images checked |

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
K_0(\\gamma_1d_{ij})+
\\frac{2e^{-H\\gamma_1}}{4+\\gamma_1^2x_{ij}^2}\\right].
```

**Reference.** O. Saad, G. Gaba, and M. Giroux, “A Closed-Form Approximation
for Ground Return Impedance of Underground Cables,” *IEEE Transactions on
Power Delivery*, 11(3), 1536–1545, 1996.

## Small-argument reduction

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filamentary mutual earth-return model with depths ``h_i,h_j`` and separation ``x``; self term uses radius ``R`` and depth ``h``. |
| Calculated quantities | Further reduced per-unit-length mutual and self earth-return impedances of buried horizontal cables |
| Earth structure | Homogeneous semi-infinite earth. |
| Model and approximation | Further approximation of source (26)–(27). It substitutes the first small-argument term ``K_0(md)=-\\{\\ln(md/2)+\\gamma\\}``, expands ``e^{-m\\ell}\\simeq1-m\\ell``, and treats ``m^2x^2\\ll4``; analogous substitutions give the self formula. Retained order is first order in the displayed exponent arguments, with higher Bessel and exponential terms discarded. |
| Main source | O. Saad, G. Gaba, and M. Giroux (1996) |
| Citation key(s) | `:Saad1996` |
| Evidence status | Original PDF page images checked |

Select `Formula(:Saad1996; approximation=:small_argument)` for equations
(31)–(32), not the default Bessel/exponential expression. The two routes retain
one publication identifier. The small-argument conditions in the source apply.

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}
\\left[-\\ln\\frac{md_{ij}}{2}-\\gamma_{\\mathrm E}
+\\frac12-\\frac{m(h_i+h_j)}{2}\\right].
```

**Reference.** [Saad1996](@cite), equations (31)–(32).
"""
description(::Formula{:Saad1996}) = "Saad underground closed form (1996)"

function propagation_constant(::Val{:Saad1996}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Saad1996})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Saad1996), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate the Saad et al. underground approximation:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}
\left[K_0(\gamma_1d_{ij})+
\frac{2e^{-(h_i+h_j)\gamma_1}}{4+\gamma_1^2x_{ij}^2}\right].
```

For a self term the direct distance is the cable radius. For a mutual
term it is the distance between the two axes. The rational correction
uses horizontal separation, including when the depths differ.

# Reference

O. Saad, G. Gaba, and M. Giroux, "A closed-form approximation for ground
return impedance of underground cables," *IEEE Transactions on Power
Delivery*, vol. 11, no. 3, pp. 1536-1545, 1996.
"""
function earth_impedance(
        ::Val{:Saad1996}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    gamma = state.gamma[2]
    radius = geometry.y_ij
    correction = 2exp(-geometry.H * gamma) / (4 + gamma^2 * radius^2)
    direct = _complex_result(
        state.jω, special_besselk(0, gamma * geometry.d_ij)
    )
    πT = one(radius) * π
    return state.jω * state.mu[1] / (2πT) * (direct + correction)
end

function earth_impedance(
        ::Val{:Saad1996}, ::Val{:small_argument}, functor, pair
)
    _require(pair, Val(:underground))
    state, geometry = functor.state, _geometry(pair)
    m = state.gamma[2]
    πT = one(geometry.d_ij) * π
    euler = one(geometry.d_ij) * Base.MathConstants.eulergamma
    return state.jω * state.mu[1] / (2πT) *
           (-log(m * geometry.d_ij / 2) - euler +
            one(geometry.d_ij) / 2 - m * geometry.H / 2)
end


:Saad1996
