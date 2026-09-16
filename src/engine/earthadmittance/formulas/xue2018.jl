function assumptions(::Val{:xue2018})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Generalized underground potential coefficient. The
registered expression is referenced to infinite earth depth; the infinite-depth expression preserves the former default.

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

**Reference.** Haoyan Xue, *General Formulation and Accurate Evaluation of
Earth-Return Parameters for Overhead / Underground Cables*, doctoral thesis,
Polytechnique Montréal, 2018.
[Primary source](https://publications.polymtl.ca/3190/1/2018_HaoyanXue.pdf).
"""
function description(::Type{<:Formula{:xue2018}}; compact::Bool=false)
    compact ? "Xue" : "Xue homogeneous-earth underground potential coefficient (2018)"
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

The underground branch uses the infinite-depth reference.
"""

function earth_potential_coefficient(
        ::Formula{:xue2018}, ::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    ratio = gamma_0_squared / gamma_1_squared
    S12 = integrate(functor.options.data.integration.method,
        SpectralIntegral(Val(:cosine),
            lambda -> begin
                u_0 = sqrt(lambda^2 + gamma_0_squared)
                u_1 = sqrt(lambda^2 + gamma_1_squared)
                exp(-geometry.H * gamma_1_squared / (u_1 + lambda)) * lambda^2 /
                ((lambda^2 + gamma_1_squared) * (u_0 + ratio * u_1))
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
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    return state.jω / (2π * kappa) *
           (direct + 2 * S12 + 2 * gamma_1^2 * S13)
end



function formulation_options(::FormulaMethod{<:Formula{:xue2018}, typeof(earth_potential_coefficient),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{2}, Val{2}}}
    return FormulationOptions((integration = (method = :quad, options = (;)),))
end

function validate(binding::FormulaMethod{<:Formula{:xue2018}, typeof(earth_potential_coefficient)},
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
