function routes(identifier::Val{:MartinsBritto2024})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        overhead = FormulaMethod(
            identifier, earth_impedance, Val(:same_medium)
        ),
        underground = FormulaMethod(
            identifier, earth_impedance, Val(:same_medium)
        ),
        mixed = FormulaMethod(identifier, earth_impedance, Val(:mixed)),
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

**Identification.** Wideband, pair-complete homogeneous-earth formulation
with explicit longitudinal propagation.

**Expression.** For two conductors in medium ``m`` with the other medium
``n``,

```math
Z_{e,ij}^{mm}=\\frac{j\\omega\\mu_m}{2\\pi}\\left[
K_0(a_md_{ij})-K_0(a_mD_{ij})+2\\int_0^\\infty
\\mu_n\\frac{e^{-a_mH}}{a_m\\mu_n+a_n\\mu_m}
\\cos(y_{ij}\\lambda)d\\lambda\\right],
```

and for a mixed pair,

```math
Z_{e,ij}^{01}=\\frac{j\\omega\\mu_0}{\\pi}\\int_0^\\infty
\\mu_1\\frac{e^{-a_0|h_i|-a_1|h_j|}}
{a_0\\mu_1+a_1\\mu_0}\\cos(y_{ij}\\lambda)d\\lambda,
\\quad a_m=\\sqrt{\\lambda^2+\\gamma_m^2+k_x^2}.
```

**Reference.** A. G. Martins-Britto, T. A. Papadopoulos, and A. I.
Chrysochos, “Transient Electromagnetic Interference Between Overhead and
Underground Conductors,” *IEEE Transactions on Electromagnetic Compatibility*,
66(3), 983–992, 2024.
"""
function description(::Formula{:MartinsBritto2024})
    "Martins-Britto, Papadopoulos, and Chrysochos wideband homogeneous-earth impedance (2024)"
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

Select the wideband impedance leaf from the physical conductor placement.

The recipe uses the same general homogeneous kernel for overhead and
underground pairs and the published transmission kernel for mixed pairs. The
three placement routes are individually replaceable without changing the
public `:MartinsBritto2024` identity.
"""
function earth_impedance(
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
Evaluate the general same-medium wideband impedance:

```math
Z_{e,ij}^{mm}=\frac{j\omega\mu_m}{2\pi}
\left[K_0(a_md_{ij})-K_0(a_mD_{ij})+
2\int_0^\infty\mu_n
\frac{e^{-a_m(h_i+h_j)}}{a_m\mu_n+a_n\mu_m}
\cos(y_{ij}\lambda)d\lambda\right],
```

where ``m`` is the conductor medium, ``n`` is the other medium, and
``a_q=\sqrt{\lambda^2+\gamma_q^2+k_x^2}``.
"""
function earth_impedance(
        ::Val{:MartinsBritto2024}, ::Val{:same_medium}, functor, pair
)
    state = functor.state
    placement = _placement(pair)
    source = typeof(placement) === Val{:overhead} ? 1 : 2
    other = source == 1 ? 2 : 1
    geometry = _geometry(pair)
    source_squared = state.gamma_medium_squared[source] + state.gamma_squared
    other_squared = state.gamma_medium_squared[other] + state.gamma_squared
    source_radial = sqrt(source_squared)
    direct = special_besselk(0, source_radial * geometry.d_ij) -
             special_besselk(0, source_radial * geometry.D_ij)
    integral = _quadrature(state) do lambda
        a_source = sqrt(lambda^2 + source_squared)
        a_other = sqrt(lambda^2 + other_squared)
        state.mu[other] * exp(-a_source * geometry.H) /
        (a_source * state.mu[other] + a_other * state.mu[source]) *
        cos(geometry.y_ij * lambda)
    end
    return state.jω * state.mu[source] / (2π) * (direct + 2 * integral)
end

raw"""
Evaluate the Martins-Britto, Papadopoulos, and Chrysochos wideband
overhead-underground mutual impedance:

```math
Z_{e,ij}^{01}=\frac{j\omega\mu_0}{\pi}\int_0^\infty
\mu_1\frac{e^{-a_0|h_i|-a_1|h_j|}}
{a_0\mu_1+a_1\mu_0}\cos(y_{ij}\lambda)d\lambda,
```

where ``a_m=\sqrt{\lambda^2+\gamma_m^2+k_x^2}``. The default ``k_x``
is exposed through `Γ`.
"""
function earth_impedance(
        ::Val{:MartinsBritto2024}, ::Val{:mixed}, functor, pair
)
    state = functor.state
    air = pair.layers[1] == 1 ? 1 : 2
    earth = air == 1 ? 2 : 1
    h_air = abs(pair.heights[air])
    h_earth = abs(pair.heights[earth])
    radial_squared = (
        state.gamma_medium_squared[1] + state.gamma_squared,
        state.gamma_medium_squared[2] + state.gamma_squared
    )
    integral = _quadrature(state) do lambda
        a_0 = sqrt(lambda^2 + radial_squared[1])
        a_1 = sqrt(lambda^2 + radial_squared[2])
        state.mu[2] * exp(-a_0 * h_air - a_1 * h_earth) /
        (a_0 * state.mu[2] + a_1 * state.mu[1]) *
        cos(pair.separation * lambda)
    end
    return state.jω * state.mu[1] / π * integral
end

:MartinsBritto2024
