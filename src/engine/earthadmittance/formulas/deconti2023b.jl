function routes(identifier::Val{:DeConti2023b})
    return (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

assumptions(::Val{:DeConti2023b}) = (
    air = _vacuum, earth = _full, permeability = vacuum_permeability
)
propagation(::Val{:DeConti2023b}) = Val(:zero)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Circular insulated cables; self uses total cable radius and equal depths. |
| Calculated quantities | Self and mutual Bessel-free small-argument entries of ``P_g``; assembled conversion ``Y_g=j\\omega P_g^{-1}`` |
| Earth structure | Homogeneous earth below air. |
| Model and approximation | Both Bessel functions in parent ``P_{g(m,n)}=j\\omega[K_0(\\gamma_1d)+\\alpha K_0(\\gamma_1D)]/[2\\pi(\\sigma_1+j\\omega\\varepsilon_1)]`` receive the leading small-argument expansion (7); logarithms are then collected to give (9). |
| Main source | A. De Conti, N. Duarte, R. Alipio, and O. E. Leal (2023) |
| Citation key(s) | `:DeConti2023b` |
| Evidence status | Original publication page images checked |

**Expression.** Equation (9); both direct and image Bessel arguments must be small. Assemble the potential coefficients before matrix inversion.

**Reference.** [DeConti2023b](@cite).
"""
description(::Formula{:DeConti2023b}) = "De Conti underground potential coefficient (2023)"

function propagation_constant(::Val{:DeConti2023b}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:DeConti2023b})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(Val(:DeConti2023b), formula, rho, epsilon, mu, jω, Γ, segments)
end

function earth_potential_coefficient(::Val{:DeConti2023b}, ::Val{:mutual}, functor, pair)
    _require(pair, Val(:underground))
    state, geometry = functor.state, _geometry(pair)
    gamma_air_squared, gamma_squared = state.gamma_medium_squared
    alpha = (gamma_squared - gamma_air_squared) / (gamma_squared + gamma_air_squared)
    gamma = state.gamma[2]
    euler = one(real(state.jω)) * Base.MathConstants.eulergamma
    bracket = log(geometry.D_ij / geometry.d_ij) -
              (alpha + 1) * (euler + log(gamma * geometry.D_ij / 2))
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    πT = one(real(state.jω)) * π
    return _complex_result(state.jω, state.jω / (2πT * kappa) * bracket)
end

:DeConti2023b
