function routes(identifier::Val{:MartinsBritto2024})
    (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        overhead = FormulaMethod(
            identifier, earth_potential_coefficient, Val(:same_medium)
        ),
        underground = FormulaMethod(
            identifier, earth_potential_coefficient, Val(:same_medium)
        ),
        mixed = FormulaMethod(identifier, earth_potential_coefficient, Val(:mixed)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:MartinsBritto2024})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:MartinsBritto2024}) = Val(:explicit)
"""
$(TYPEDSIGNATURES)

**Identification.** Wideband, pair-complete potential-coefficient formulation
with explicit longitudinal propagation.

**Expression.** For conductors in medium ``m``,

```math
P_{e,ij}^{mm}=\\frac{j\\omega}{2\\pi\\kappa_m}\\left[
K_0(a_md_{ij})-K_0(a_mD_{ij})+2\\int_0^\\infty
I_{ij}^{mm}(\\lambda)\\cos(y_{ij}\\lambda)d\\lambda\\right],
\\quad \\kappa_m=\\sigma_m+j\\omega\\varepsilon_m.
```

For a mixed pair,

```math
P_{e,ij}^{01}=\\frac{j\\omega}{\\pi(\\sigma_0+j\\omega\\varepsilon_0)}
\\int_0^\\infty I_{ij}^{01,MPC}(\\lambda)\\cos(y_{ij}\\lambda)d\\lambda,
```

```math
I_{ij}^{01,MPC}=\\gamma_0^2\\mu_1
\\frac{a_0\\mu_0+a_1\\mu_1}
{(a_0\\gamma_1^2\\mu_0+a_1\\gamma_0^2\\mu_1)
(a_0\\mu_1+a_1\\mu_0)}e^{-a_0|h_i|-a_1|h_j|}.
```

**Reference.** A. G. Martins-Britto, T. A. Papadopoulos, and A. I.
Chrysochos, “Transient Electromagnetic Interference Between Overhead and
Underground Conductors,” *IEEE Transactions on Electromagnetic Compatibility*,
66(3), 983–992, 2024.
"""
function description(::Formula{:MartinsBritto2024})
    "Martins-Britto, Papadopoulos, and Chrysochos wideband homogeneous-earth potential coefficient (2024)"
end

function propagation_constant(
        ::Val{:MartinsBritto2024}, jω, permeability, permittivity
)
    squared = oftype(jω, (-jω^2) * permeability * permittivity)
    return (Γ = sqrt(squared), squared)
end

function (formula::Formula{:MartinsBritto2024})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:MartinsBritto2024), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

"""
$(TYPEDSIGNATURES)

Select the wideband potential-coefficient leaf from conductor placement.

The recipe uses the same general homogeneous kernel for overhead and
underground pairs and the published transmission kernel for mixed pairs. The
three placement routes are individually replaceable without changing the
public `:MartinsBritto2024` identity.
"""
function earth_potential_coefficient(
        ::Val{:MartinsBritto2024}, ::Val{:mutual}, functor, pair
)
    placement = _placement(pair)
    typeof(placement) === Val{:overhead} &&
        return functor.routes.overhead(functor, pair)
    typeof(placement) === Val{:underground} &&
        return functor.routes.underground(functor, pair)
    return functor.routes.mixed(functor, pair)
end

raw"""
Evaluate the general same-medium wideband potential coefficient:

```math
P_{e,ij}^{mm}=\frac{j\omega}{2\pi\kappa_m}
\left[K_0(a_md_{ij})-K_0(a_mD_{ij})+
2\int_0^\infty I_{ij}^{mm}(\lambda)
\cos(y_{ij}\lambda)d\lambda\right],
```

where ``\kappa_m=\sigma_m+j\omega\varepsilon_m`` and
``a_q=\sqrt{\lambda^2+\gamma_q^2+k_x^2}``.
"""
function earth_potential_coefficient(
        ::Val{:MartinsBritto2024}, ::Val{:same_medium}, functor, pair
)
    state = functor.state
    placement = _placement(pair)
    source = typeof(placement) === Val{:overhead} ? 1 : 2
    other = source == 1 ? 2 : 1
    geometry = _geometry(pair)
    gamma_source_squared = state.gamma_medium_squared[source]
    gamma_other_squared = state.gamma_medium_squared[other]
    source_squared = gamma_source_squared + state.gamma_squared
    other_squared = gamma_other_squared + state.gamma_squared
    source_radial = sqrt(source_squared)
    direct = special_besselk(0, source_radial * geometry.d_ij) -
             special_besselk(0, source_radial * geometry.D_ij)
    integral = _quadrature(state) do lambda
        a_source = sqrt(lambda^2 + source_squared)
        a_other = sqrt(lambda^2 + other_squared)
        common = a_source * state.mu[other] + a_other * state.mu[source]
        decay = exp(-a_source * geometry.H)
        magnetic = state.mu[other] * decay / common
        coupling_numerator = state.mu[other] * state.mu[source] *
                             a_source *
                             (gamma_source_squared - gamma_other_squared) *
                             decay
        coupling_denominator = common * (
            a_source * gamma_other_squared * state.mu[source] +
            a_other * gamma_source_squared * state.mu[other]
        )
        (magnetic + coupling_numerator / coupling_denominator) *
        cos(geometry.y_ij * lambda)
    end
    kappa = state.sigma[source] + state.jω * state.epsilon[source]
    return state.jω / (2π * kappa) * (direct + 2 * integral)
end

raw"""
Evaluate the Martins-Britto, Papadopoulos, and Chrysochos wideband mixed-media
potential coefficient:

```math
P_{e,ij}^{01}=\frac{j\omega}{\pi(\sigma_0+j\omega\varepsilon_0)}
\int_0^\infty I_{ij}^{01,MPC}(\lambda)\cos(y_{ij}\lambda)d\lambda,
```

```math
I_{ij}^{01,MPC}=\gamma_0^2\mu_1
\frac{a_0\mu_0+a_1\mu_1}
{(a_0\gamma_1^2\mu_0+a_1\gamma_0^2\mu_1)
(a_0\mu_1+a_1\mu_0)}
e^{-a_0|h_i|-a_1|h_j|}.
```

The implementation preserves the product—not the typographical sum—of the
two parenthesized denominator factors in the source equation. The
longitudinal ``k_x`` is exposed through `Γ`.
"""
function earth_potential_coefficient(
        ::Val{:MartinsBritto2024}, ::Val{:mixed}, functor, pair
)
    state = functor.state
    air = pair.layers[1] == 1 ? 1 : 2
    earth = air == 1 ? 2 : 1
    h_air = abs(pair.heights[air])
    h_earth = abs(pair.heights[earth])
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    radial_squared = (
        gamma_0_squared + state.gamma_squared,
        gamma_1_squared + state.gamma_squared
    )
    integral = _quadrature(state) do lambda
        a_0 = sqrt(lambda^2 + radial_squared[1])
        a_1 = sqrt(lambda^2 + radial_squared[2])
        numerator = gamma_0_squared * state.mu[2] *
                    (a_0 * state.mu[1] + a_1 * state.mu[2])
        denominator = (
            a_0 * gamma_1_squared * state.mu[1] +
            a_1 * gamma_0_squared * state.mu[2]
        ) * (a_0 * state.mu[2] + a_1 * state.mu[1])
        numerator / denominator * exp(-a_0 * h_air - a_1 * h_earth) *
        cos(pair.separation * lambda)
    end
    kappa_0 = state.sigma[1] + state.jω * state.epsilon[1]
    return state.jω / (π * kappa_0) * integral
end

:MartinsBritto2024
