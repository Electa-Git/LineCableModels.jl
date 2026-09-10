function assumptions(::Val{:Wise1934})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

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
"""
description(::Formula{:Wise1934}) = "Wise1934 homogeneous-earth overhead impedance"

Γ(::Val{:Wise1934}, jω, materials, layers) = zero(jω)

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
        ::Val{:Wise1934}, ::Union{Val{:self}, Val{:mutual}}, ::Val{1}, ::Val{1},
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
            angle = min(pi/4, atan(float(nominal(geometry.H / (2geometry.y_ij))))),
            features = retained_earth_features(state, geometry, workspace)),
        functor.options.integration.options, workspace)
    return state.jω * state.mu[1] / (2π) *
           (log(geometry.D_ij / geometry.d_ij) + 2 * integral)
end

Formulation(::LineCableModelsCoaxial, selected::Formula{:Wise1934}) = selected

function hooks(::FormulaMethod{:Wise1934, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{1}, Val{1}}}
    return (configurable = (:Γ, :air, :earth, :permeability, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:Wise1934), Γ),
            air = FormulaMethod(Val(:full), propagation),
            earth = FormulaMethod(Val(:full), propagation),
            permeability = identity,
            contribution = nothing))
end

function computation_options(::FormulaMethod{:Wise1934, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{1}, Val{1}}}
    (integration = (method = :quad, options = (;)),)
end

function validate(binding::FormulaMethod{:Wise1934, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:default})
    binding
end

:Wise1934
