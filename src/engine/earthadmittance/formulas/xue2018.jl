function routes(identifier::Val{:Xue2018b})
    (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        infinite = FormulaMethod(
            identifier, earth_potential_coefficient, Val(:infinite_depth)
        ),
        surface = FormulaMethod(
            identifier, earth_potential_coefficient, Val(:surface_reference)
        ),
        penetration = FormulaMethod(
            identifier, earth_potential_coefficient, Val(:penetration_depth)
        ),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Xue2018b})
    (
        air = _full,
        earth = _full,
        permeability = vacuum_permeability,
        quadrature = :adaptive
    )
end

propagation(::Val{:Xue2018b}) = Val(:zero)

function Formula(::Val{:Xue2018b};quadrature::Symbol=:adaptive,kwargs...)
    quadrature in (:adaptive,:double_exponential) ||
        throw(ArgumentError(":Xue2018b quadrature must be :adaptive or :double_exponential"))
    identifier=Val(:Xue2018b); defaults=routes(identifier); overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) ||
        throw(ArgumentError("unknown routes for earthadmittance formula :Xue2018b"))
    selected=merge(defaults,overrides)
    values=merge(assumptions(identifier),(;quadrature))
    return Formula{:Xue2018b,typeof(selected),typeof(values)}(selected,values)
end
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Circular insulated cables; self uses equal depths and outer radius; mutual uses centers. |
| Calculated quantities | Complete-field earth-return admittance matrix and qTEM self/mutual earth-return potential coefficients ``P_e`` with ``Y=j\\omega P^{-1}`` after total matrix assembly |
| Earth structure | Homogeneous earth beneath homogeneous air. |
| Model and approximation | Equation (34) is obtained from the complete-field admittance by setting ``k_v=0``; no integral is approximated further. The full-field formula retains modal ``k_v`` and the infinite-depth voltage reference. |
| Main source | H. Xue, A. Ametani, J. Mahseredjian, and I. Kocar (2018) |
| Citation key(s) | `:Xue2018b` |
| Evidence status | Original publication page images checked |

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Self/mutual matrix entries evaluated individually. Geometry enters parent integrand. |
| Calculated quantities | Numerical quadrature for earth potential/admittance cable integrals |
| Earth structure | Structure enters parent integrand. |
| Model and approximation | Finite node truncation only; no change to the selected admittance physics. |
| Main source | D. L. Pires, F. A. Moreira, M. G. Soares, and F. M. Vasconcellos (2026) |
| Citation key(s) | `:Pires2026` |
| Evidence status | Page-image verified |

**Numerical evaluation.** `quadrature=:double_exponential` selects the
Pires node sum on the existing physical kernels. The default
`quadrature=:adaptive` is unchanged. Positive spectral coordinates use
the appendix-map interpretation documented in the survey; refinement
and tail checks report nonconvergence without a silent fallback.

**Numerical scope.** The default evaluator supplies the quasi-TEM potential coefficient (34) referred to infinite depth. Surface and penetration-depth references remain selectable support terms.

**Expression.**

```math
P_{e,ij}^{\\infty}=\\frac{j\\omega}{2\\pi(\\sigma_1+j\\omega\\varepsilon_1)}
\\left[K_0(\\gamma_1d_{ij})-K_0(\\gamma_1D_{ij})+2S_{12}^c+
2\\gamma_1^2S_{13}^c\\right],
```

```math
S_{12}^c=\\int_0^\\infty\\frac{e^{-Hu_1}\\lambda^2\\cos(y\\lambda)}
{(\\lambda^2+\\gamma_1^2)[u_0+(\\gamma_0^2/\\gamma_1^2)u_1]}d\\lambda,
\\quad
S_{13}^c=\\int_0^\\infty\\frac{e^{-Hu_1}\\cos(y\\lambda)}
{(\\lambda^2+\\gamma_1^2)(u_0+u_1)}d\\lambda.
```

The surface reference is
``P_{e,ij}^{0}=P_{e,ij}^{\\infty}-P_{ec,ij}``; the penetration-depth reference
uses the same structure at the exact lossy-medium penetration depth.

**Reference.** H. Xue, *Electromagnetic Transients in Large HV Cable
Networks*, doctoral thesis, Delft University of Technology, 2018; equations
as consolidated in Ametani et al., IET, 2021.
"""
function description(::Formula{:Xue2018b})
    "Xue et al. generalized underground potential coefficient (2018)"
end

function propagation_constant(::Val{:Xue2018b}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Xue2018b})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Xue2018b), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Xue et al. underground potential coefficient referred to
infinite earth depth:

```math
P_{e,ij}^{\infty}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
\left[K_0(\gamma_1d_{ij})-K_0(\gamma_1D_{ij})+2S_{12}^c+
2\gamma_1^2S_{13}^c\right],
```

```math
S_{12}^c=\int_0^\infty\frac{e^{-Hu_1}\lambda^2\cos(y\lambda)}
{(\lambda^2+\gamma_1^2)[u_0+(\gamma_0^2/\gamma_1^2)u_1]}d\lambda,
\quad
S_{13}^c=\int_0^\infty\frac{e^{-Hu_1}\cos(y\lambda)}
{(\lambda^2+\gamma_1^2)(u_0+u_1)}d\lambda.
```

This infinite-depth reference is the registered default. The surface and
penetration-depth references are support leaves of the same `:Xue2018b`
entry, available through [`routes`](@ref) without creating more formula IDs.
"""
function integral_terms(
        ::Val{:Xue2018b}, state, geometry, height = geometry.H
)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    ratio = gamma_0_squared / gamma_1_squared
    S12 = _quadrature(state) do lambda
        u_0 = spectral_root(lambda^2 + gamma_0_squared,state.jω)
        u_1 = spectral_root(lambda^2 + gamma_1_squared,state.jω)
        exp(-height * u_1) * lambda^2 * cos(geometry.y_ij * lambda) /
        ((lambda^2 + gamma_1_squared) * (u_0 + ratio * u_1))
    end
    S13 = _quadrature(state) do lambda
        u_0 = spectral_root(lambda^2 + gamma_0_squared,state.jω)
        u_1 = spectral_root(lambda^2 + gamma_1_squared,state.jω)
        exp(-height * u_1) * cos(geometry.y_ij * lambda) /
        ((lambda^2 + gamma_1_squared) * (u_0 + u_1))
    end
    return S12, S13
end

function earth_potential_coefficient(
        identifier::Val{:Xue2018b}, ::Val{:mutual}, functor, pair
)
    return earth_potential_coefficient(
        identifier, Val(:infinite_depth), functor, pair
    )
end

function earth_potential_coefficient(
        identifier::Val{:Xue2018b}, ::Val{:infinite_depth}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    S12, S13 = integral_terms(identifier, state, geometry)
    gamma_1 = state.gamma[2]
    direct = special_besselk(0, gamma_1 * geometry.d_ij) -
             special_besselk(0, gamma_1 * geometry.D_ij)
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    return _complex_result(state.jω,state.jω / (2*(one(real(state.jω))*π) * kappa) *
           (direct + 2 * S12 + 2 * gamma_1^2 * S13))
end

raw"""
Evaluate the surface-referenced support coefficient from Xue et al.:

```math
P_{e,ij}^{0}=P_{e,ij}^{\infty}-P_{ec,ij},\qquad
P_{ec,ij}=\frac{j\omega}{\pi(\sigma_1+j\omega\varepsilon_1)}
(S_{14}^c+\gamma_1^2S_{15}^c),
```

where ``S_{14}^c`` and ``S_{15}^c`` use the infinite-depth denominators and
the attenuation ``e^{-\frac12(h_i+h_j)u_1}``.
"""
function earth_potential_coefficient(
        identifier::Val{:Xue2018b}, ::Val{:surface_reference}, functor, pair
)
    # S12 + gamma1^2*S13 = integral u0/[u1*(u0+n^2*u1)].
    # The surface subtraction is therefore exactly Magalhaes (2018).
    return earth_potential_coefficient(Val(:Magalhaes2018),Val(:mutual),functor,pair)
end

raw"""
Evaluate the penetration-depth-referenced support coefficient from Xue et al.:

```math
P_{e,ij}^{\delta}=P_{e,ij}^{\infty}-P_{\delta_e,ij},
```

```math
h_r=-\left\{\frac{\omega^2\varepsilon_1\mu_0}{2}
\left[\sqrt{1+(\sigma_1/(\omega\varepsilon_1))^2}-1\right]
\right\}^{-1/2}.
```

The correction retains the corpus distances ``d_{\delta,ij}``,
``D_{\delta,ij}`` and integrals ``S_{16}^c``, ``S_{17}^c`` verbatim.
"""
function earth_potential_coefficient(
        identifier::Val{:Xue2018b}, ::Val{:penetration_depth}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    S12, S13 = integral_terms(identifier, state, geometry)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    ratio = gamma_0_squared / gamma_1_squared
    # This is exactly the printed skin-depth radical, without cancellation
    # in sqrt(1+x^2)-1 for displacement-dominated earth.
    attenuation=real(state.gamma[2])
    attenuation>0 || throw(DomainError(attenuation,
        "a finite penetration-depth reference requires lossy earth"))
    h_r=-inv(attenuation)
    exponent_height = geometry.H / 2 - h_r
    S16 = _quadrature(state) do lambda
        u_0 = spectral_root(lambda^2 + gamma_0_squared,state.jω)
        u_1 = spectral_root(lambda^2 + gamma_1_squared,state.jω)
        exp(-exponent_height * u_1) * lambda^2 * cos(geometry.y_ij * lambda) /
        ((lambda^2 + gamma_1_squared) * (u_0 + ratio * u_1))
    end
    S17 = _quadrature(state) do lambda
        u_0 = spectral_root(lambda^2 + gamma_0_squared,state.jω)
        u_1 = spectral_root(lambda^2 + gamma_1_squared,state.jω)
        exp(-exponent_height * u_1) * cos(geometry.y_ij * lambda) /
        ((lambda^2 + gamma_1_squared) * (u_0 + u_1))
    end
    d_delta = hypot(geometry.y_ij, h_r - geometry.H / 2)
    D_delta = hypot(geometry.y_ij, h_r + geometry.H / 2)
    gamma_1 = state.gamma[2]
    direct = special_besselk(0, gamma_1 * geometry.d_ij) -
             special_besselk(0, gamma_1 * geometry.D_ij)
    reference = special_besselk(0, gamma_1 * d_delta) -
                special_besselk(0, gamma_1 * D_delta) +
                2 * S16 + 2 * gamma_1_squared * S17
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    infinite = direct + 2 * S12 + 2 * gamma_1_squared * S13
    return _complex_result(state.jω,state.jω /
        (2*(one(real(state.jω))*π) * kappa) * (infinite - reference))
end

:Xue2018b
