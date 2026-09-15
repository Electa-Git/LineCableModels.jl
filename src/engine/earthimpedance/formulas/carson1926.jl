function assumptions(::Val{:carson1926})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Classical homogeneous, conductive-earth overhead
impedance. Displacement currents and longitudinal propagation are neglected.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+2\\int_0^\\infty
\\frac{e^{-H\\lambda}\\cos(y_{ij}\\lambda)}
{\\lambda+\\sqrt{\\lambda^2+\\gamma_g^2}}d\\lambda\\right],
\\qquad \\gamma_g^2=j\\omega\\mu_0\\sigma_g.
```

**Reference.** J. R. Carson, “Wave Propagation in Overhead Wires with Ground
Return,” *Bell System Technical Journal*, 5, 539–554, 1926.
"""
description(::Type{<:Formula{:carson1926}}; compact::Bool=false) = compact ? "Carson" : "Carson homogeneous-earth overhead impedance (1926)"


raw"""
Evaluate Carson's homogeneous-earth overhead impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[
\ln\frac{D_{ij}}{d_{ij}}+
2\int_0^\infty
\frac{e^{-(h_i+h_j)\lambda}\cos(y_{ij}\lambda)}
{\lambda+\sqrt{\lambda^2+\gamma_g^2}}\,d\lambda\right],
```

where ``\gamma_g^2=j\omega\mu_0\sigma_g``. The original Carson
assumptions neglect earth displacement current, air displacement current in
the correction, and longitudinal propagation.


# Reference

J. R. Carson, "Wave propagation in overhead wires with ground return,"
*Bell System Technical Journal*, vol. 5, pp. 539-554, 1926.
"""
function earth_impedance(
        ::Formula{:carson1926}, ::Union{Val{:self}, Val{:mutual}}, ::Val{1}, ::Val{1},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    gamma_squared = state.gamma_medium_squared[2]
    integral = integrate(functor.options.integration.method,
        SpectralIntegral(
            Val(:cosine), lambda -> begin
                attenuation = sqrt(lambda^2 + gamma_squared)
                one(lambda) /
                (lambda + attenuation)
            end,
            (height = geometry.H, separation = geometry.y_ij),
            float(nominal(abs(state.gamma[2])))),
        functor.options.integration.options, workspace)
    πT = one(geometry.H) * π
    return state.jω * state.mu[1] / (2πT) *
           (log(geometry.D_ij / geometry.d_ij) + 2 * integral)
end



function computation_options(::FormulaMethod{<:Formula{:carson1926}, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{1}, Val{1}}}
    (integration = (method = :quad, options = (;)),)
end

function validate(binding::FormulaMethod{<:Formula{:carson1926}, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:bottommost})
    binding
end

"""
$(TYPEDSIGNATURES)

Evaluate this formulation's medium state: absolute permeability \\[H/m\\]
and transverse propagation constant \\[1/m\\]. Air retains its prescribed
permeability; soil follows the selected source's magnetic approximation.
"""
function constitutive(::Formula{:carson1926}, ::Val{:air}, jω, μ, σ, ε)
    return (mu=μ, gamma=propagation(Val(:lossless), jω, μ, σ, ε))
end

function constitutive(::Formula{:carson1926}, ::Val{:earth}, jω, μ, σ, ε)
    permeability = vacuum_permeability(μ)
    return (mu=permeability, gamma=propagation(Val(:conductive), jω, permeability, σ, ε))
end

:carson1926
