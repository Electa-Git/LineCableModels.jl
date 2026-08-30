function routes(::Val{:Lucca1994})
    (
        self = lucca1994,
        mutual = lucca1994,
        overhead = carson1926,
        underground = pollaczek1926_underground,
        mixed = lucca1994_mixed,
        Γ = lucca1994_gamma
    )
end

function assumptions(::Val{:Lucca1994})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Lucca1994}) = Val(:zero)
description(::Formula{:Lucca1994}) = "Lucca"

lucca1994_gamma(jω, permeability, permittivity) = (Γ = zero(jω), squared = zero(jω))

function (formula::Formula{:Lucca1994})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Lucca1994), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

function lucca1994(functor, pair)
    return lucca1994(functor, pair, _placement(pair))
end

function lucca1994(functor, pair, ::Val{:overhead})
    return functor.routes.overhead(functor, pair)
end

function lucca1994(functor, pair, ::Val{:underground})
    return functor.routes.underground(functor, pair)
end

function lucca1994(functor, pair, ::Val{:mixed})
    return functor.routes.mixed(functor, pair)
end

raw"""
Evaluate Lucca's approximation for the mutual impedance between one overhead
and one buried conductor:

```math
Z_{e,ij}^{01}=\frac{j\omega\mu_0}{2\pi}\left[
\ln\frac{S}{D}-\frac23\left(\frac{h_e}{S^2}\right)^3
H(H^2-3y_{ij}^2)\right],
```

```math
h_e=\frac1{\sqrt{j\omega\mu_0\sigma_g}},\qquad
H=h_a+h_g+2h_e,\qquad
S=\sqrt{H^2+y_{ij}^2},\qquad
D=\sqrt{(h_a+h_g)^2+y_{ij}^2}.
```

The complete recipe retains Pollaczek's exact homogeneous same-medium leaves
and replaces only the mixed interaction.

# Reference

G. Lucca, "Mutual impedance between an overhead and a buried line with earth
return," *9th International Conference on Electromagnetic Compatibility*,
1994. DOI: 10.1049/cp:19940679.
"""
function lucca1994_mixed(functor, pair)
    state = functor.state
    air = pair.layers[1] == 1 ? 1 : 2
    earth = air == 1 ? 2 : 1
    h_a = abs(pair.heights[air])
    h_g = abs(pair.heights[earth])
    h_e = inv(state.gamma[2])
    H = h_a + h_g + 2h_e
    S_squared = H^2 + pair.separation^2
    S = sqrt(S_squared)
    D = hypot(pair.separation, h_a + h_g)
    correction = (2 * one(h_a) / 3) * (h_e / S_squared)^3 *
                 H * (H^2 - 3 * pair.separation^2)
    πT = one(h_a) * π
    return state.jω * state.mu[1] / (2πT) * (log(S / D) - correction)
end

:Lucca1994
