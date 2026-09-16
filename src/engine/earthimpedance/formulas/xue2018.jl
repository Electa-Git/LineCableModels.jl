function assumptions(::Val{:xue2018})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

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

**Reference.** Haoyan Xue, *General Formulation and Accurate Evaluation of
Earth-Return Parameters for Overhead / Underground Cables*, doctoral thesis,
Polytechnique Montréal, 2018.
[Primary source](https://publications.polymtl.ca/3190/1/2018_HaoyanXue.pdf).
"""
description(::Type{<:Formula{:xue2018}}; compact::Bool=false) = compact ? "Xue" : "Xue homogeneous-earth underground impedance (2018)"


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
        ::Formula{:xue2018}, ::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    S11 = integrate(functor.options.data.integration.method,
        SpectralIntegral(Val(:cosine),
            lambda -> begin
                u_0 = sqrt(lambda^2 + gamma_0_squared)
                u_1 = sqrt(lambda^2 + gamma_1_squared)
                exp(-geometry.H * gamma_1_squared / (u_1 + lambda)) * lambda^2 /
                ((lambda^2 + gamma_1_squared) * (u_0 + u_1))
            end,
            (height = geometry.H, separation = geometry.y_ij),
            float(nominal(abs(state.gamma[2])));
            angle = min(pi/4, atan(float(nominal(geometry.H / (2geometry.y_ij))))),
            features = retained_earth_features(state, geometry, workspace)),
        functor.options.data.integration.options, workspace)
    S13 = integrate(functor.options.data.integration.method,
        SpectralIntegral(Val(:cosine),
            lambda -> begin
                u_0 = sqrt(lambda^2 + gamma_0_squared)
                u_1 = sqrt(lambda^2 + gamma_1_squared)
                exp(-geometry.H * gamma_1_squared / (u_1 + lambda)) /
                ((lambda^2 + gamma_1_squared) * (u_0 + u_1))
            end,
            (height = geometry.H, separation = geometry.y_ij),
            float(nominal(abs(state.gamma[2])));
            angle = min(pi/4, atan(float(nominal(geometry.H / (2geometry.y_ij))))),
            features = retained_earth_features(state, geometry, workspace)),
        functor.options.data.integration.options, workspace)
    gamma_1 = state.gamma[2]
    direct = special_besselk(0, gamma_1 * geometry.d_ij) -
             special_besselk(0, gamma_1 * geometry.D_ij)
    return state.jω * state.mu[1] / (2π) *
           (direct + 2 * S11 + 2 * gamma_1_squared * S13)
end



function formulation_options(::FormulaMethod{<:Formula{:xue2018}, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{2}, Val{2}}}
    return FormulationOptions((integration = (method = :quad, options = (;)),))
end

function validate(binding::FormulaMethod{<:Formula{:xue2018}, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:bottommost})
    binding
end

"""
$(TYPEDSIGNATURES)

Evaluate this formulation's medium state: absolute permeability \\[H/m\\]
and transverse propagation constant \\[1/m\\]. Air retains its prescribed
permeability; soil follows the selected source's magnetic approximation.
"""
function constitutive(::Formula{:xue2018}, ::Val{:air}, jω, μ, σ, ε)
    return (mu=μ, gamma=propagation(Val(:full), jω, μ, σ, ε))
end

function constitutive(::Formula{:xue2018}, ::Val{:earth}, jω, μ, σ, ε)
    permeability = vacuum_permeability(μ)
    return (mu=permeability, gamma=propagation(Val(:full), jω, permeability, σ, ε))
end

:xue2018
