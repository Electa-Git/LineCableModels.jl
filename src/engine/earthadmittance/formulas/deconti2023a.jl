function routes(identifier::Val{:DeConti2023a})
    return (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

assumptions(::Val{:DeConti2023a}) = (
    air = _vacuum, earth = _full, permeability = vacuum_permeability
)
propagation(::Val{:DeConti2023a}) = Val(:zero)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Insulated circular cables; self uses ``h_m=h_n`` and outer insulated radius ``r=r_o``. |
| Calculated quantities | Self and mutual entries of ground-return potential-coefficient matrix ``P_g`` and source-prescribed assembled admittance ``Y_g=j\\omega P_g^{-1}`` |
| Earth structure | Homogeneous earth below homogeneous air. |
| Model and approximation | The source applies the asymptote ``u_0/(u_0+\\gamma_0^2\\gamma_1^{-2}u_1)\\approx\\gamma_1^2/(\\gamma_1^2+\\gamma_0^2)`` for ``\\|\\gamma_1\\|\\gg\\|\\gamma_0\\|`` to the parent integral ``\\Theta_2``. The remaining Fourier–Bessel integral is evaluated as ``K_0(\\gamma_1D)``, yielding (14) and (15). |
| Main source | A. De Conti, N. Duarte, and R. Alipio (2023) |
| Citation key(s) | `:DeConti2023a` |
| Evidence status | Original publication page images checked |

**Expression.** Equation (15); the earth bulk constant must dominate the air bulk constant. Assemble the potential coefficients before matrix inversion.

**Reference.** [DeConti2023a](@cite).
"""
description(::Formula{:DeConti2023a}) = "De Conti underground potential coefficient (2023)"

function propagation_constant(::Val{:DeConti2023a}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:DeConti2023a})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(Val(:DeConti2023a), formula, rho, epsilon, mu, jω, Γ, segments)
end

function earth_potential_coefficient(::Val{:DeConti2023a}, ::Val{:mutual}, functor, pair)
    _require(pair, Val(:underground))
    state, geometry = functor.state, _geometry(pair)
    gamma_air_squared, gamma_squared = state.gamma_medium_squared
    alpha = (gamma_squared - gamma_air_squared) / (gamma_squared + gamma_air_squared)
    gamma = state.gamma[2]
    bracket = special_besselk(0, gamma * geometry.d_ij) +
              alpha * special_besselk(0, gamma * geometry.D_ij)
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    πT = one(real(state.jω)) * π
    return _complex_result(state.jω, state.jω / (2πT * kappa) * bracket)
end

:DeConti2023a
