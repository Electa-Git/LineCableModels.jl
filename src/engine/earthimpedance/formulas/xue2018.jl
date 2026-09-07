function routes(identifier::Val{:Xue2018b})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
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
        throw(ArgumentError("unknown routes for earthimpedance formula :Xue2018b"))
    selected=merge(defaults,overrides)
    values=merge(assumptions(identifier),(;quadrature))
    return Formula{:Xue2018b,typeof(selected),typeof(values)}(selected,values)
end
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Circular cables with outer radius ``r_{on}``; self evaluates field at ``h_n=h_m`` and ``d_n-d_m=r_{on}``, mutual at cable centers. |
| Calculated quantities | Complete-field modal earth-return impedance matrix and its self/mutual quasi-TEM underground-cable specialization |
| Earth structure | Homogeneous earth half-space below homogeneous air. |
| Model and approximation | Equation (33) follows from the complete-field formula solely by setting each unknown modal ``k_v=0``. It is an exact spectral representation within that qTEM reduction, not a closed-form approximation of its integrals. |
| Main source | H. Xue, A. Ametani, J. Mahseredjian, and I. Kocar (2018) |
| Citation key(s) | `:Xue2018b` |
| Evidence status | Original publication page images checked for parent matrix, all displayed kernels and quasi-TEM formula |

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Cable self/mutual entries evaluated independently. Geometry enters parent integrands, not the transform. |
| Calculated quantities | Numerical quadrature for cable earth-return impedance integrals |
| Earth structure | The quadrature is structure-agnostic; applications use the stated cable media. |
| Model and approximation | Finite implementation truncates the formally infinite node sum; the physical parent kernel is unchanged. |
| Main source | D. L. Pires, F. A. Moreira, M. G. Soares, and F. M. Vasconcellos (2026) |
| Citation key(s) | `:Pires2026` |
| Evidence status | Page-image verified |

**Numerical evaluation.** `quadrature=:double_exponential` selects the
Pires node sum on the existing physical kernels. The default
`quadrature=:adaptive` is unchanged. Positive spectral coordinates use
the appendix-map interpretation documented in the survey; refinement
and tail checks report nonconvergence without a silent fallback.

**Numerical scope.** The evaluator supplies the quasi-TEM reduction (33), with the longitudinal modal wavenumber set to zero.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
K_0(\\gamma_1d_{ij})-K_0(\\gamma_1D_{ij})+2S_{11}^c+
2\\gamma_1^2S_{13}^c\\right],
```

```math
S_{11}^c=\\int_0^\\infty\\frac{e^{-Hu_1}\\lambda^2\\cos(y\\lambda)}
{(\\lambda^2+\\gamma_1^2)(u_0+u_1)}d\\lambda,\\qquad
S_{13}^c=\\int_0^\\infty\\frac{e^{-Hu_1}\\cos(y\\lambda)}
{(\\lambda^2+\\gamma_1^2)(u_0+u_1)}d\\lambda.
```

**Reference.** H. Xue, *Electromagnetic Transients in Large HV Cable
Networks*, doctoral thesis, Delft University of Technology, 2018; equations
as consolidated in Ametani et al., IET, 2021.
"""
function description(::Formula{:Xue2018b})
    "Xue et al. generalized homogeneous-earth underground impedance (2018)"
end

function propagation_constant(::Val{:Xue2018b}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Xue2018b})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Xue2018b), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate the Xue et al. generalized homogeneous-earth underground impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}
\left[K_0(\gamma_1d_{ij})-K_0(\gamma_1D_{ij})+2S_{11}^c+
2\gamma_1^2S_{13}^c\right],
```

```math
S_{11}^c=\int_0^\infty\frac{e^{-Hu_1}\lambda^2\cos(y\lambda)}
{(\lambda^2+\gamma_1^2)(u_0+u_1)}d\lambda,\qquad
S_{13}^c=\int_0^\infty\frac{e^{-Hu_1}\cos(y\lambda)}
{(\lambda^2+\gamma_1^2)(u_0+u_1)}d\lambda,
```

where ``u_m=\sqrt{\lambda^2+\gamma_m^2}``.
"""
function earth_impedance(
        ::Val{:Xue2018b}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    integral = _quadrature(state) do lambda
        u_0 = spectral_root(lambda^2 + gamma_0_squared,state.jω)
        u_1 = spectral_root(lambda^2 + gamma_1_squared,state.jω)
        # S11 + gamma_1_squared*S13 cancels the common u_1^2 factor.
        exp(-geometry.H * u_1) * cos(geometry.y_ij * lambda) / (u_0 + u_1)
    end
    gamma_1 = state.gamma[2]
    direct = special_besselk(0, gamma_1 * geometry.d_ij) -
             special_besselk(0, gamma_1 * geometry.D_ij)
    return _complex_result(state.jω,state.jω * state.mu[1] /
        (2*(one(real(state.jω))*π)) * (direct + 2 * integral))
end

:Xue2018b
