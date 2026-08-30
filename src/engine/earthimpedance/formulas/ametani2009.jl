function routes(::Val{:Ametani2009})
    (
        self = ametani2009,
        mutual = ametani2009,
        overhead = carson1926,
        underground = pollaczek1926_underground,
        mixed = ametani2009_mixed,
        Γ = ametani2009_gamma
    )
end

function assumptions(::Val{:Ametani2009})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Ametani2009}) = Val(:zero)
description(::Formula{:Ametani2009}) = "Ametani"

ametani2009_gamma(jω, permeability, permittivity) = (Γ = zero(jω), squared = zero(jω))

function (formula::Formula{:Ametani2009})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Ametani2009), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

function ametani2009(functor, pair)
    return ametani2009(functor, pair, _placement(pair))
end

function ametani2009(functor, pair, ::Val{:overhead})
    return functor.routes.overhead(functor, pair)
end

function ametani2009(functor, pair, ::Val{:underground})
    return functor.routes.underground(functor, pair)
end

function ametani2009(functor, pair, ::Val{:mixed})
    return functor.routes.mixed(functor, pair)
end

raw"""
Evaluate Ametani's approximation for the mutual impedance between one
overhead and one buried conductor:

```math
Z_{e,ij}^{01}=\frac{j\omega\mu_0}{2\pi}
e^{-h_g/h_e}\ln\frac{S}{D},
```

```math
h_e=\frac1{\sqrt{j\omega\mu_0\sigma_g}},\qquad
S=\sqrt{(h_a+h_g+2h_e)^2+y_{ij}^2},\qquad
D=\sqrt{(h_a+h_g)^2+y_{ij}^2}.
```

Here ``h_a`` and ``h_g`` are positive height and burial-depth magnitudes.
The complete recipe retains Pollaczek's exact homogeneous same-medium leaves
and replaces only the mixed interaction.

# Reference

A. Ametani, "An investigation of earth-return impedance between overhead
and underground conductors and its approximation," *IEEE Transactions on
Electromagnetic Compatibility*, vol. 51, pp. 860-867, 2009.
DOI: 10.1109/TEMC.2009.2019953.
"""
function ametani2009_mixed(functor, pair)
    state = functor.state
    air = pair.layers[1] == 1 ? 1 : 2
    earth = air == 1 ? 2 : 1
    h_a = abs(pair.heights[air])
    h_g = abs(pair.heights[earth])
    h_e = inv(state.gamma[2])
    D = hypot(pair.separation, h_a + h_g)
    S = sqrt((h_a + h_g + 2h_e)^2 + pair.separation^2)
    πT = one(h_a) * π
    return state.jω * state.mu[1] / (2πT) *
           exp(-h_g / h_e) * log(S / D)
end

:Ametani2009
