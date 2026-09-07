function routes(identifier::Val{:DeConti2023a})
    return (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

assumptions(::Val{:DeConti2023a}) = (
    air = _lossless, earth = _full, permeability = vacuum_permeability
)
propagation(::Val{:DeConti2023a}) = Val(:zero)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Insulated circular cables at depths ``h_m,h_n``; mutual horizontal separation ``r``; self sets ``h_m=h_n`` and ``r=r_o``, the external radius including insulation. |
| Calculated quantities | Self and mutual closed-form approximation of Xue's homogeneous-earth underground-cable ground-return impedance |
| Earth structure | Homogeneous earth below homogeneous air. |
| Model and approximation | The source begins with ``Z_{g(m,n)}=j\\omega\\mu_0[\\Lambda+\\Theta_1]/(2\\pi)`` and approximates the square-root ratio in ``d\\Theta_1/dH`` by a constant plus an exponentially decaying term. It then replaces ``e^{-H\\sqrt{\\lambda^2+\\gamma_1^2}}`` by ``e^{-H\\gamma_1}`` only in that residual and integrates using the Bessel identity (12), yielding (13). No series order or formal remainder bound is provided. |
| Main source | A. De Conti, N. Duarte, and R. Alipio (2023) |
| Citation key(s) | `:DeConti2023a` |
| Evidence status | Original publication page images checked |

**Expression.** Equation (13); finite air propagation and full earth admittivity.

**Reference.** [DeConti2023a](@cite).
"""
description(::Formula{:DeConti2023a}) = "De Conti underground impedance (2023)"

function propagation_constant(::Val{:DeConti2023a}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:DeConti2023a})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(Val(:DeConti2023a), formula, rho, epsilon, mu, jω, Γ, segments)
end

function earth_impedance(::Val{:DeConti2023a}, ::Val{:mutual}, functor, pair)
    _require(pair, Val(:underground))
    state, geometry = functor.state, _geometry(pair)
    gamma_air, gamma = state.gamma
    image = (gamma - gamma_air) / (gamma + gamma_air) *
            exp(-geometry.H * gamma) * 2 / (4 + (gamma * geometry.y_ij)^2)
    direct = special_besselk(0, gamma * geometry.d_ij)
    πT = one(real(state.jω)) * π
    return _complex_result(state.jω, state.jω * state.mu[1] / (2πT) * (direct + image))
end

:DeConti2023a
