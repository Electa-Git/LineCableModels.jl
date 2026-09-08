function assumptions(::Val{:default})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

Package default homogeneous-earth potential coefficient for overhead or underground conductor pairs.
Resolve the placement before preparing frequency state. The identity and all
native routes remain `:default`; mixed placement has no default expression.

**Identification.** Wideband homogeneous-earth overhead potential
coefficient.

**Expression.**

```math
P_{e,ij}=\\frac{P_{0,ij}+M_{ij}+jN_{ij}}{2\\pi\\varepsilon_0},
```

```math
M_{ij}+jN_{ij}=2\\int_0^\\infty
\\frac{e^{-H\\lambda}\\cos(y_{ij}\\lambda)}
{(\\gamma_1^2/\\gamma_0^2)\\lambda+
\\sqrt{\\lambda^2+\\gamma_1^2-\\gamma_0^2}}d\\lambda,
\\quad P_{0,ij}=\\ln(D_{ij}/d_{ij}).
```

**Reference.** W. H. Wise, “Potential Coefficients for Ground Return
Circuits,” *Bell System Technical Journal*, 27, 365–371, 1948.

**Identification.** Generalized underground potential coefficient. The
registered expression is referenced to infinite earth depth; the infinite-depth expression is the retained package default.

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

**Reference.** H. Xue, *Electromagnetic Transients in Large HV Cable
Networks*, doctoral thesis, Delft University of Technology, 2018; equations
as consolidated in Ametani et al., IET, 2021.
"""
description(::Formula{:default}) = "Package default homogeneous-earth potential coefficient"

function Γ(::Val{:default}, jω, materials, layers)
    return zero(jω)
end

raw"""
Evaluate Wise's wideband overhead earth potential coefficient:

```math
P_{e,ij}=\frac{P_{0,ij}+M_{ij}+jN_{ij}}{2\pi\varepsilon_0},
```

```math
M_{ij}+jN_{ij}=2\int_0^\infty
\frac{e^{-(h_i+h_j)\lambda}\cos(y_{ij}\lambda)}
{(\gamma_1^2/\gamma_0^2)\lambda+
\sqrt{\lambda^2+\gamma_1^2-\gamma_0^2}}d\lambda,
\qquad P_{0,ij}=\ln(D_{ij}/d_{ij}).
```
"""
function earth_potential_coefficient(
        ::Val{:default}, ::Union{Val{:self}, Val{:mutual}}, ::Val{1}, ::Val{1},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    ratio = gamma_1_squared / gamma_0_squared
    integral = integrate(functor.options.integration.method,
        SpectralIntegral(Val(:cosine),
            lambda -> begin
                radial = sqrt(lambda^2 + gamma_1_squared - gamma_0_squared)
                origin = sqrt(gamma_1_squared - gamma_0_squared)
                # Exact pole extraction; rationalization avoids subtracting close kernels.
                -lambda^2 /
                ((origin + radial) * (ratio*lambda + radial) * (ratio*lambda + origin))
            end,
            (height = geometry.H, separation = geometry.y_ij),
            float(nominal(abs(sqrt(gamma_1_squared - gamma_0_squared))));
            pole = (residue = inv(ratio),
                location = sqrt(gamma_1_squared-gamma_0_squared)/ratio),
            angle = min(pi/4, atan(float(nominal(geometry.H / (2geometry.y_ij)))))),
        functor.options.integration.options, workspace)
    return (log(geometry.D_ij / geometry.d_ij) + 2 * integral) /
           (2π * state.epsilon[1])
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
        ::Val{:default}, ::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    ratio = gamma_0_squared / gamma_1_squared
    S12 = integrate(functor.options.integration.method,
        SpectralIntegral(Val(:cosine),
            lambda -> begin
                u_0 = sqrt(lambda^2 + gamma_0_squared)
                u_1 = sqrt(lambda^2 + gamma_1_squared)
                exp(-geometry.H * gamma_1_squared / (u_1 + lambda)) * lambda^2 /
                ((lambda^2 + gamma_1_squared) * (u_0 + ratio * u_1))
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
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    return state.jω / (2π * kappa) *
           (direct + 2 * S12 + 2 * gamma_1^2 * S13)
end

Formulation(::LineCableModelsCoaxial, selected::Formula{:default}) = selected

function hooks(::FormulaMethod{:default, typeof(earth_potential_coefficient),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{1}, Val{1}}}
    return (configurable = (:Γ, :air, :earth, :permeability, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:default), Γ),
            air = FormulaMethod(Val(:full), propagation),
            earth = FormulaMethod(Val(:full), propagation),
            permeability = vacuum_permeability,
            contribution = nothing))
end

function computation_options(::FormulaMethod{:default, typeof(earth_potential_coefficient),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{1}, Val{1}}}
    (integration = (method = :quad, options = (;)),)
end

function hooks(::FormulaMethod{:default, typeof(earth_potential_coefficient),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{2}, Val{2}}}
    return (configurable = (:Γ, :air, :earth, :permeability, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:default), Γ),
            air = FormulaMethod(Val(:full), propagation),
            earth = FormulaMethod(Val(:full), propagation),
            permeability = vacuum_permeability,
            contribution = nothing))
end

function computation_options(::FormulaMethod{:default, typeof(earth_potential_coefficient),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{2}, Val{2}}}
    (integration = (method = :quad, options = (;)),)
end

function validate(binding::FormulaMethod{:default, typeof(earth_potential_coefficient)},
        ::EquivalentHomogeneous.Formula{:default})
    binding
end

:default
