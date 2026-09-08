function assumptions(::Val{:default})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

Package default homogeneous-earth impedance for overhead or underground conductor pairs.
Resolve the placement before preparing frequency state. The identity and all
native routes remain `:default`; mixed placement has no default expression.

**Identification.** Homogeneous-earth wideband overhead integral retaining
earth displacement current and magnetic permeability.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+2\\int_0^\\infty
\\frac{\\mu_1e^{-\\lambda H}}
{\\lambda\\mu_1+a_1\\mu_0}\\cos(y_{ij}\\lambda)d\\lambda\\right],
\\quad a_1=\\sqrt{\\lambda^2+\\gamma_1^2-\\gamma_0^2}.
```

**Reference.** W. H. Wise, “Propagation of High-Frequency Currents in Ground
Return Circuits,” *Proceedings of the Institute of Radio Engineers*, 22,
522–527, 1934.

**Identification.** Generalized homogeneous-earth underground wideband
impedance.

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
description(::Formula{:default}) = "Package default homogeneous-earth impedance"

function Γ(::Val{:default}, jω, materials, layers)
    return zero(jω)
end

raw"""
Evaluate Wise's homogeneous-earth overhead impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[\ln\frac{D_{ij}}{d_{ij}}+
2\int_0^\infty F_{ij}^{W}(\lambda)\cos(y_{ij}\lambda)\,d\lambda\right],
```

```math
F_{ij}^{W}=\frac{\mu_1e^{-\lambda(h_i+h_j)}}
{\lambda\mu_1+a_1\mu_0},\qquad
a_1=\sqrt{\lambda^2+\gamma_1^2-\gamma_0^2}.
```
"""
function earth_impedance(
        ::Val{:default}, ::Union{Val{:self}, Val{:mutual}}, ::Val{1}, ::Val{1},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    contrast = state.gamma_medium_squared[2] - state.gamma_medium_squared[1]
    integral = integrate(functor.options.integration.method,
        SpectralIntegral(Val(:cosine),
            lambda -> begin
                a_1 = sqrt(lambda^2 + contrast)
                state.mu[2] /
                (lambda * state.mu[2] + a_1 * state.mu[1])
            end,
            (height = geometry.H, separation = geometry.y_ij),
            float(nominal(abs(state.gamma[2])));
            angle = min(pi/4, atan(float(nominal(geometry.H / (2geometry.y_ij)))))),
        functor.options.integration.options, workspace)
    return state.jω * state.mu[1] / (2π) *
           (log(geometry.D_ij / geometry.d_ij) + 2 * integral)
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
        ::Val{:default}, ::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    S11 = integrate(functor.options.integration.method,
        SpectralIntegral(Val(:cosine),
            lambda -> begin
                u_0 = sqrt(lambda^2 + gamma_0_squared)
                u_1 = sqrt(lambda^2 + gamma_1_squared)
                exp(-geometry.H * gamma_1_squared / (u_1 + lambda)) * lambda^2 /
                ((lambda^2 + gamma_1_squared) * (u_0 + u_1))
            end,
            (height = geometry.H, separation = geometry.y_ij),
            float(nominal(abs(state.gamma[2])));
            angle = min(pi/4, atan(float(nominal(geometry.H / (2geometry.y_ij)))))),
        functor.options.integration.options, workspace)
    S13 = integrate(functor.options.integration.method,
        SpectralIntegral(Val(:cosine),
            lambda -> begin
                u_0 = sqrt(lambda^2 + gamma_0_squared)
                u_1 = sqrt(lambda^2 + gamma_1_squared)
                exp(-geometry.H * gamma_1_squared / (u_1 + lambda)) /
                ((lambda^2 + gamma_1_squared) * (u_0 + u_1))
            end,
            (height = geometry.H, separation = geometry.y_ij),
            float(nominal(abs(state.gamma[2])));
            angle = min(pi/4, atan(float(nominal(geometry.H / (2geometry.y_ij)))))),
        functor.options.integration.options, workspace)
    gamma_1 = state.gamma[2]
    direct = special_besselk(0, gamma_1 * geometry.d_ij) -
             special_besselk(0, gamma_1 * geometry.D_ij)
    return state.jω * state.mu[1] / (2π) *
           (direct + 2 * S11 + 2 * gamma_1_squared * S13)
end

Formulation(::LineCableModelsCoaxial, selected::Formula{:default}) = selected

function hooks(::FormulaMethod{:default, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{1}, Val{1}}}
    return (configurable = (:Γ, :air, :earth, :permeability, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:default), Γ),
            air = FormulaMethod(Val(:full), propagation),
            earth = FormulaMethod(Val(:full), propagation),
            permeability = identity,
            contribution = nothing))
end

function computation_options(::FormulaMethod{:default, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{1}, Val{1}}}
    (integration = (method = :quad, options = (;)),)
end

function hooks(::FormulaMethod{:default, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{2}, Val{2}}}
    return (configurable = (:Γ, :air, :earth, :permeability, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:default), Γ),
            air = FormulaMethod(Val(:full), propagation),
            earth = FormulaMethod(Val(:full), propagation),
            permeability = vacuum_permeability,
            contribution = nothing))
end

function computation_options(::FormulaMethod{:default, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{2}, Val{2}}}
    (integration = (method = :quad, options = (;)),)
end

function validate(binding::FormulaMethod{:default, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:default})
    binding
end

:default
