function assumptions(::Val{:Pollaczek1926})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Classical homogeneous-earth underground integral.

**Expression.** The underground term is

```math
Z_{e,ij}^{11}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
K_0(\\gamma_1d_{ij})-K_0(\\gamma_1D_{ij})+2\\int_0^\\infty
\\frac{e^{-H\\sqrt{\\lambda^2+\\gamma_1^2}}}
{\\lambda+\\sqrt{\\lambda^2+\\gamma_1^2}}
\\cos(y_{ij}\\lambda)d\\lambda\\right],
```

**Reference.** F. Pollaczek, “Über das Feld einer unendlich langen
wechselstromdurchflossenen Einfachleitung,” *Elektrische Nachrichtentechnik*,
3, 339–360, 1926.
"""
function description(::Formula{:Pollaczek1926})
    "Pollaczek homogeneous-earth underground impedance (1926)"
end

function Γ(
        ::Val{:Pollaczek1926}, jω, materials, layers
)
    zero(jω)
end

raw"""
Evaluate Pollaczek's homogeneous-earth underground impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[K_0(\gamma_1d_{ij})-
K_0(\gamma_1D_{ij})+2\int_0^\infty
\frac{e^{-(h_i+h_j)\sqrt{\lambda^2+\gamma_1^2}}}
{\lambda+\sqrt{\lambda^2+\gamma_1^2}}
\cos(y_{ij}\lambda)\,d\lambda\right].
```
"""
function earth_impedance(
        ::Val{:Pollaczek1926}, ::Union{Val{:self}, Val{:mutual}}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    gamma = state.gamma[2]
    integral = integrate(functor.options.integration.method,
        SpectralIntegral(Val(:cosine),
            lambda -> begin
                u_1 = sqrt(lambda^2 + gamma^2)
                exp(-geometry.H * gamma^2 / (u_1 + lambda)) /
                (lambda + u_1)
            end,
            (height = geometry.H, separation = geometry.y_ij),
            float(nominal(abs(state.gamma[2])))),
        functor.options.integration.options, workspace)
    direct = oftype(
        state.jω,
        special_besselk(0, gamma * geometry.d_ij) -
        special_besselk(0, gamma * geometry.D_ij)
    )
    πT = one(geometry.H) * π
    return state.jω * state.mu[1] / (2πT) * (direct + 2 * integral)
end

Formulation(::LineCableModelsCoaxial, selected::Formula{:Pollaczek1926}) = selected

function hooks(::FormulaMethod{:Pollaczek1926, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{2}, Val{2}}}
    return (configurable = (:Γ, :earth, :permeability, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:Pollaczek1926), Γ),
            air = FormulaMethod(Val(:lossless), propagation),
            earth = FormulaMethod(Val(:conductive), propagation),
            permeability = vacuum_permeability,
            contribution = nothing))
end

function computation_options(::FormulaMethod{:Pollaczek1926, typeof(earth_impedance),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{2}, Val{2}}}
    (integration = (method = :quad, options = (;)),)
end

function validate(binding::FormulaMethod{:Pollaczek1926, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:default})
    binding
end

:Pollaczek1926
