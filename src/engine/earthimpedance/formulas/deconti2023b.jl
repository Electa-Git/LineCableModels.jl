function routes(identifier::Val{:DeConti2023b})
    return (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

assumptions(::Val{:DeConti2023b}) = (
    air = _lossless, earth = _full, permeability = vacuum_permeability
)
propagation(::Val{:DeConti2023b}) = Val(:zero)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Circular insulated cable positions; self replaces ``r`` by total cable radius including jacket and sets equal depths. |
| Calculated quantities | Self and mutual Bessel-free small-argument approximation of the 2023 De Conti–Duarte–Alipio closed form |
| Earth structure | Homogeneous earth below air. |
| Model and approximation | The only additional operation on parent equation (1) is ``K_0(z)\\approx-\\ln(z/2)-\\gamma_E`` for ``0<z\\ll1``. All other terms are retained unchanged. |
| Main source | A. De Conti, N. Duarte, R. Alipio, and O. E. Leal (2023) |
| Citation key(s) | `:DeConti2023b` |
| Evidence status | Original publication page images checked |

**Expression.** Equation (8); the direct-distance Bessel argument must be small.

**Reference.** [DeConti2023b](@cite).
"""
description(::Formula{:DeConti2023b}) = "De Conti underground impedance (2023)"

function propagation_constant(::Val{:DeConti2023b}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:DeConti2023b})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(Val(:DeConti2023b), formula, rho, epsilon, mu, jω, Γ, segments)
end

function earth_impedance(::Val{:DeConti2023b}, ::Val{:mutual}, functor, pair)
    _require(pair, Val(:underground))
    state, geometry = functor.state, _geometry(pair)
    gamma_air, gamma = state.gamma
    image = (gamma - gamma_air) / (gamma + gamma_air) *
            exp(-geometry.H * gamma) * 2 / (4 + (gamma * geometry.y_ij)^2)
    direct = -log(gamma * geometry.d_ij / 2) - one(real(state.jω)) * Base.MathConstants.eulergamma
    πT = one(real(state.jω)) * π
    return _complex_result(state.jω, state.jω * state.mu[1] / (2πT) * (direct + image))
end

:DeConti2023b
