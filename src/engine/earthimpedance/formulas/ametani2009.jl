function routes(identifier::Val{:Ametani2009})
    (
        self = formula_method(identifier, earth_impedance, Val(:self)),
        mutual = formula_method(identifier, earth_impedance, Val(:mutual)),
        overhead = formula_method(identifier, earth_impedance, Val(:overhead)),
        underground = formula_method(
            identifier, earth_impedance, Val(:underground)
        ),
        mixed = formula_method(identifier, earth_impedance, Val(:mixed)),
        Γ = formula_method(identifier, propagation_constant)
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
"""
$(TYPEDSIGNATURES)

**Identification.** Pair-complete homogeneous-earth recipe using Carson for
overhead pairs, Pollaczek for underground pairs, and Ametani's approximation
for mixed overhead-underground pairs.

**Expression.** Its distinctive mixed term is

```math
Z_{e,ij}^{01}=\\frac{j\\omega\\mu_0}{2\\pi}e^{-h_g/h_e}\\ln\\frac{S}{D},
\\quad h_e=(j\\omega\\mu_0\\sigma_g)^{-1/2},
```

```math
S=\\sqrt{(h_a+h_g+2h_e)^2+y_{ij}^2},\\qquad
D=\\sqrt{(h_a+h_g)^2+y_{ij}^2}.
```

**Reference.** A. Ametani, T. Yoneda, Y. Baba, and N. Nagaoka, “An
Investigation of Earth-Return Impedance Between Overhead and Underground
Conductors and Its Approximation,” *IEEE Transactions on Electromagnetic
Compatibility*, 51, 860–867, 2009.
DOI: 10.1109/TEMC.2009.2019953.
"""
description(::Formula{:Ametani2009}) =
    "Ametani pair-complete homogeneous-earth impedance (2009)"

function propagation_constant(::Val{:Ametani2009}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Ametani2009})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Ametani2009), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

function earth_impedance(
        ::Val{:Ametani2009}, ::Val{:mutual}, functor, pair
)
    placement = _placement(pair)
    typeof(placement) === Val{:overhead} &&
        return functor.routes.overhead(functor, pair)
    typeof(placement) === Val{:underground} &&
        return functor.routes.underground(functor, pair)
    return functor.routes.mixed(functor, pair)
end

function earth_impedance(
        ::Val{:Ametani2009}, ::Val{:overhead}, functor, pair
)
    return earth_impedance(Val(:Carson1926), Val(:mutual), functor, pair)
end

function earth_impedance(
        ::Val{:Ametani2009}, ::Val{:underground}, functor, pair
)
    return earth_impedance(Val(:Pollaczek1926), Val(:underground), functor, pair)
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
function earth_impedance(
        ::Val{:Ametani2009}, ::Val{:mixed}, functor, pair
)
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
