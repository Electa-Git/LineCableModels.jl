function routes(identifier::Val{:Pawlik2018})
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

function assumptions(::Val{:Pawlik2018})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Pawlik2018}) = Val(:explicit)
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinite parallel thin conductors placed in opposite homogeneous half-spaces. |
| Calculated quantities | Full-spectrum cross-boundary mutual per-unit-length impedance |
| Earth structure | Two homogeneous half-spaces with independent conductivity, permittivity, and permeability. |
| Model and approximation | Integral representation within the full TM/TE model retaining ``Γ``; thin-wire and infinite-line assumptions remain. |
| Main source | B. Pawlik, D. Woodhouse, and T. J. Summers (2018) |
| Citation key(s) | `:Pawlik2018` |
| Evidence status | Original publication page images checked; apparent printed LHS index defect in (88) retained |

**Numerical scope.** The source contribution is the mixed term. The recipe also supplies same-medium terms for assembling a complete matrix; those terms are not attributed as new mixed formulas to the 2018 paper.

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

**Propagation convention.** The engine input is the longitudinal spectral
wave number ``k_x``. For Pawlik's ``e^{-\\Gamma z}``, supply
``k_x=j\\Gamma``; then ``a_m^2=\\lambda^2+\\gamma_m^2+k_x^2``
matches equations (39)–(40). The inherited default numerical selection is
``k_x^2=\\omega^2\\mu_1\\varepsilon_1``; it is not a solved mode or a
source-imposed value. Supply an explicit value to compare a chosen mode;
``k_x=0`` gives the Dawalibi low-frequency restriction.

**Reference.** [Pawlik2018](@cite), equations (88)–(89).
[MartinsBritto2024](@cite), equation (5), provides an equivalent later
mixed-impedance expression. Its existing impedance selector remains an
alias, while its separately registered potential-matrix formula is retained.
[Pawlik2020](@cite), equations (14), (16)–(17), supplies an equivalent
same-medium self witness at the outer insulation radius. It retains the
horizontal surface sample in the cosine and image distance, and does not
require a new external series registration.
"""
function description(::Formula{:Pawlik2018})
    "Pawlik, Woodhouse, and Summers generalized mixed impedance (2018)"
end

function propagation_constant(
        ::Val{:Pawlik2018}, jω, permeability, permittivity
)
    squared = oftype(jω, (-jω^2) * permeability * permittivity)
    return (Γ = sqrt(squared), squared)
end

function (formula::Formula{:Pawlik2018})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Pawlik2018), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

"""
$(TYPEDSIGNATURES)

Select the wideband impedance leaf from the physical conductor placement.

The recipe uses the same general homogeneous kernel for overhead and
underground pairs and the published transmission kernel for mixed pairs. The
three placement routes are individually replaceable without changing the
public `:Pawlik2018` identity.
"""
function earth_impedance(
        ::Val{:Pawlik2018}, ::Val{:mutual}, functor, pair
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
        ::Val{:Pawlik2018}, ::Val{:same_medium}, functor, pair
)
    state = functor.state
    placement = _placement(pair)
    source = typeof(placement) === Val{:overhead} ? 1 : 2
    other = source == 1 ? 2 : 1
    geometry = _geometry(pair)
    source_squared = state.gamma_medium_squared[source] + state.gamma_squared
    other_squared = state.gamma_medium_squared[other] + state.gamma_squared
    source_radial = spectral_root(source_squared,state.jω)
    direct = iszero(source_radial) ? log(geometry.D_ij/geometry.d_ij) :
        special_besselk(0, source_radial * geometry.d_ij) -
        special_besselk(0, source_radial * geometry.D_ij)
    integral = _quadrature(state) do lambda
        a_source = spectral_root(lambda^2 + source_squared,state.jω)
        a_other = spectral_root(lambda^2 + other_squared,state.jω)
        state.mu[other] * exp(-a_source * geometry.H) /
        (a_source * state.mu[other] + a_other * state.mu[source]) *
        cos(geometry.y_ij * lambda)
    end
    return state.jω * state.mu[source] / (2π) * (direct + 2 * integral)
end

raw"""
Evaluate the Pawlik, Woodhouse, and Summers generalized
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
        ::Val{:Pawlik2018}, ::Val{:mixed}, functor, pair
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
        a_0 = spectral_root(lambda^2 + radial_squared[1],state.jω)
        a_1 = spectral_root(lambda^2 + radial_squared[2],state.jω)
        state.mu[2] * exp(-a_0 * h_air - a_1 * h_earth) /
        (a_0 * state.mu[2] + a_1 * state.mu[1]) *
        cos(pair.separation * lambda)
    end
    return state.jω * state.mu[1] / π * integral
end

:Pawlik2018
