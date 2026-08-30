function routes(::Val{:Pollaczek1926})
    (
        self = pollaczek1926,
        mutual = pollaczek1926,
        overhead = carson1926,
        underground = pollaczek1926_underground,
        mixed = pollaczek1926_mixed,
        Γ = pollaczek1926_gamma
    )
end

function assumptions(::Val{:Pollaczek1926})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Pollaczek1926}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Pair-complete classical homogeneous-earth recipe:
Carson overhead, Pollaczek underground, and the exact mixed integral.

**Expression.** The underground and mixed terms are

```math
Z_{e,ij}^{11}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
K_0(\\gamma_1d_{ij})-K_0(\\gamma_1D_{ij})+2\\int_0^\\infty
\\frac{e^{-H\\sqrt{\\lambda^2+\\gamma_1^2}}}
{\\lambda+\\sqrt{\\lambda^2+\\gamma_1^2}}
\\cos(y_{ij}\\lambda)d\\lambda\\right],
```

```math
Z_{e,ij}^{01}=\\frac{j\\omega\\mu_0}{\\pi}\\int_0^\\infty
\\frac{\\mu_1e^{-\\lambda|h_i|-a_1|h_j|}}
{\\lambda\\mu_1+a_1\\mu_0}\\cos(y_{ij}\\lambda)d\\lambda,
\\qquad a_1=\\sqrt{\\lambda^2+\\gamma_1^2}.
```

**Reference.** F. Pollaczek, “Über das Feld einer unendlich langen
wechselstromdurchflossenen Einfachleitung,” *Elektrische Nachrichtentechnik*,
3, 339–360, 1926.
"""
function description(::Formula{:Pollaczek1926})
    "Pollaczek homogeneous-earth overhead, underground, and mixed impedance (1926)"
end

function pollaczek1926_gamma(jω, permeability, permittivity)
    (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Pollaczek1926})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Pollaczek1926), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

"""
$(TYPEDSIGNATURES)

Select Pollaczek's leaf impedance from the physical conductor placement.

The registered recipe contains the homogeneous-earth overhead, underground,
and overhead-underground interactions. Its public identity therefore does not
encode which leaf the pair requires. The `overhead`, `underground`, and
`mixed` routes remain individually replaceable when composing an experiment.
"""
function pollaczek1926(functor, pair)
    return pollaczek1926(functor, pair, _placement(pair))
end

function pollaczek1926(functor, pair, ::Val{:overhead})
    return functor.routes.overhead(functor, pair)
end

function pollaczek1926(functor, pair, ::Val{:underground})
    return functor.routes.underground(functor, pair)
end

function pollaczek1926(functor, pair, ::Val{:mixed})
    return functor.routes.mixed(functor, pair)
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
function pollaczek1926_underground(functor, pair)
    state = functor.state
    geometry = _geometry(pair)
    gamma = state.gamma[2]
    integral = _quadrature(state) do lambda
        u_1 = sqrt(lambda^2 + gamma^2)
        exp(-geometry.H * u_1) * cos(geometry.y_ij * lambda) /
        (lambda + u_1)
    end
    direct = _complex_result(
        state.jω,
        special_besselk(0, gamma * geometry.d_ij) -
        special_besselk(0, gamma * geometry.D_ij)
    )
    πT = one(geometry.H) * π
    return state.jω * state.mu[1] / (2πT) * (direct + 2 * integral)
end

raw"""
Evaluate Pollaczek's overhead-underground mutual impedance:

```math
Z_{e,ij}^{01}=\frac{j\omega\mu_0}{\pi}\int_0^\infty
\frac{\mu_1e^{-\lambda|h_i|-a_1|h_j|}}
{\lambda\mu_1+a_1\mu_0}\cos(y_{ij}\lambda)d\lambda,
\qquad a_1=\sqrt{\lambda^2+\gamma_1^2}.
```
"""
function pollaczek1926_mixed(functor, pair)
    state = functor.state
    air = pair.layers[1] == 1 ? 1 : 2
    earth = air == 1 ? 2 : 1
    h_air = abs(pair.heights[air])
    h_earth = abs(pair.heights[earth])
    integral = _quadrature(state) do lambda
        a_1 = sqrt(lambda^2 + state.gamma_medium_squared[2])
        state.mu[2] * exp(-lambda * h_air - a_1 * h_earth) /
        (lambda * state.mu[2] + a_1 * state.mu[1]) *
        cos(pair.separation * lambda)
    end
    πT = one(h_air) * π
    return state.jω * state.mu[1] / πT * integral
end

:Pollaczek1926
