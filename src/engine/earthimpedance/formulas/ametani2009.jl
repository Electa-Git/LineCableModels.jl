function routes(identifier::Val{:Ametani2009})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        overhead = FormulaMethod(identifier, earth_impedance, Val(:overhead)),
        underground = FormulaMethod(
            identifier, earth_impedance, Val(:underground)
        ),
        mixed = FormulaMethod(identifier, earth_impedance, Val(:mixed)),
        Γ = FormulaMethod(identifier, propagation_constant)
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
function Formula(::Val{:Ametani2009};approximation::Symbol=:image,kwargs...)
    approximation in (:image,:power_frequency) ||
        throw(ArgumentError(":Ametani2009 approximation must be :image or :power_frequency"))
    identifier=Val(:Ametani2009)
    defaults=routes(identifier)
    if approximation===:power_frequency
        defaults=merge(defaults,(;mixed=FormulaMethod(
            identifier,earth_impedance,Val(:power_frequency))))
    end
    overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) ||
        throw(ArgumentError("unknown routes for earth-impedance formula :Ametani2009"))
    return Formula(identifier,merge(defaults,overrides),assumptions(identifier))
end
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel horizontal isolated conductors with heights/depths ``h_1,h_2`` and transverse separation ``y``. Their radii are drawn in Fig. 1 but do not enter this mutual earth term. |
| Calculated quantities | Mutual earth-return impedance per unit length between one overhead and one buried conductor |
| Earth structure | Homogeneous conducting half-space below homogeneous air with a planar interface. |
| Model and approximation | The parent is the Pollaczek mixed integral as restated in (5)–(6). The authors first substitute ``\\sqrt{s^2+m^2}\\simeq s+m`` in the exponent, (17), then differentiate with respect to ``y``, change variables ``s=mt``, deform the contour under a stated no-pole assumption, and substitute ``t/(\\sqrt{t^2+1}+t)\\simeq(1-e^{-2t})/2``, (24). Integration and the subsequent antiderivative in ``y`` lead to (27). No controlled retained/discarded series order or uniform expansion parameter is assigned to these two substitutions. A further power-frequency approximation is recorded separately. |
| Main source | Akihiro Ametani, Tetsuzo Yoneda, Yoshihiro Baba, and Naoto Nagaoka (2009) |
| Citation key(s) | `:Ametani2009` |
| Evidence status | Original PDF equations checked against page images; source coordinate convention unresolved. |

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

## Power-frequency reduction

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel horizontal wires, horizontal separation ``y``, signed source coordinates ``h_1,h_2``; ``d^2=(h_1-h_2)^2+y^2`` in section II. |
| Calculated quantities | Low-frequency mutual earth-return impedance between overhead and buried conductors |
| Earth structure | Homogeneous half-space below air. |
| Model and approximation | Leading low-frequency reduction of the exponential image (27), requiring the complex penetration depth to be large relative to the conductor heights and separation. The numerical coefficient retains the parent's SI scale and full frequency factor; the source supplies no error remainder. |
| Main source | Ametani, Yoneda, Baba, and Nagaoka (2009), who identify this reduction with a previously known overhead approximation rather than claim priority for that older expression |
| Citation key(s) | `:Ametani2009` |
| Evidence status | Original PDF equations checked against page images; numerical low-frequency coefficient obtained from the printed parent (27). |

**Expression.** Select `Formula(:Ametani2009; approximation=:power_frequency)`
for the leading low-frequency limit of (27):

```math
Z_{e,ij}^{01}\\simeq\\frac{j\\omega\\mu_0}{2\\pi}
\\ln\\frac{2}{g d},\\qquad
g=\\sqrt{j\\omega\\mu_0\\sigma_g},\\qquad
d=\\sqrt{x^2+(h_a+h_g)^2}.
```

Both heights here are positive magnitudes. The condition is
``|g|\\max(h_a,h_g,d)\\ll1``, not a dimensional comparison with unity.
The parent fixes the frequency multiplier of the complete complex
coefficient. With the source's units and rounding, its positive-frequency
form is ``f\\{1+j[8.253+0.628\\ln(\\rho_g/(fd^2))]\\}``
in milliohms per kilometre. The evaluator retains the unrounded SI constants.

The printed equations (31)–(32) remain in the source record. No bare
dimensionless ``1/12`` is added to an impedance. Lucca's existing
full image correction supplies its own dimensioned contribution.
The default `approximation=:image` retains the exponential image;
same-medium routes remain unchanged in either selection.

**Reference.** [Ametani2009](@cite), equations (27), (30)–(32).
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

function earth_impedance(::Val{:Ametani2009},::Val{:power_frequency},functor,pair)
    _require(pair,Val(:mixed))
    state=functor.state
    d=hypot(pair.separation,sum(abs,pair.heights))
    isfinite(d) && d>0 && !iszero(state.jω) ||
        throw(DomainError(d,":Ametani2009 requires positive separation and nonzero frequency"))
    return _complex_result(state.jω,state.jω*state.mu[1]/(2*(one(d)*π))*
        log(2/(state.gamma[2]*d)))
end

:Ametani2009
