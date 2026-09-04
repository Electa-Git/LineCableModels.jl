function routes(identifier::Val{:Pollaczek1926})
    (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        overhead = FormulaMethod(
            identifier, earth_potential_coefficient, Val(:overhead)
        ),
        underground = FormulaMethod(
            identifier, earth_potential_coefficient, Val(:underground)
        ),
        mixed = FormulaMethod(identifier, earth_potential_coefficient, Val(:mixed)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Pollaczek1926})
    (
        air = _vacuum,
        earth = _vacuum,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Pollaczek1926}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Classical pair-complete potential-coefficient recipe:
electrostatic image coefficient overhead, the classical underground
coefficient, and zero mixed-media coefficient.

**Expression.**

```math
P_{e,ij}^{00}=\\frac{1}{2\\pi\\varepsilon_0}\\ln\\frac{D_{ij}}{d_{ij}},
```

```math
P_{e,ij}^{11}=\\frac{j\\omega}{2\\pi(\\sigma_1+j\\omega\\varepsilon_1)}
[K_0(\\gamma_0d_{ij})-K_0(\\gamma_0D_{ij})],\\qquad
P_{e,ij}^{01}=0.
```

**Reference.** F. Pollaczek, “Über das Feld einer unendlich langen
wechselstromdurchflossenen Einfachleitung,” *Elektrische Nachrichtentechnik*,
3, 339–360, 1926; potential-coefficient transcription follows Ametani et al.,
IET, 2021.
"""
function description(::Formula{:Pollaczek1926})
    "Pollaczek classical overhead, underground, and mixed potential coefficients (1926)"
end

function propagation_constant(
        ::Val{:Pollaczek1926}, jω, permeability, permittivity
)
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

Select Pollaczek's leaf potential coefficient from conductor placement.

The same registered recipe supplies same-medium self and mutual terms together
with the classical zero mixed-media coefficient. The `overhead`,
`underground`, and `mixed` leaves remain individually replaceable.
"""
function earth_potential_coefficient(
        ::Val{:Pollaczek1926}, ::Val{:mutual}, functor, pair
)
    placement = _placement(pair)
    typeof(placement) === Val{:overhead} &&
        return functor.routes.overhead(functor, pair)
    typeof(placement) === Val{:underground} &&
        return functor.routes.underground(functor, pair)
    return functor.routes.mixed(functor, pair)
end

raw"""
Evaluate the classical overhead electrostatic-image coefficient:

```math
P_{e,ij}^{00}=\frac{1}{2\pi\varepsilon_0}
\ln\frac{D_{ij}}{d_{ij}}.
```

For a self term, ``d_{ii}`` is the conductor outer radius carried by the
canonical earth pair.
"""
function earth_potential_coefficient(
        ::Val{:Pollaczek1926}, ::Val{:overhead}, functor, pair
)
    state = functor.state
    geometry = _geometry(pair)
    return log(geometry.D_ij / geometry.d_ij) / (2π * state.epsilon[1])
end

raw"""
Evaluate Pollaczek's classical homogeneous-earth underground potential
coefficient:

```math
P_{e,ij}^{11}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
\left[K_0(\gamma_0d_{ij})-K_0(\gamma_0D_{ij})\right],
\qquad \gamma_0=j\omega\sqrt{\mu_0\varepsilon_0}.
```

This is the classical reference used with the zero mixed-media coupling; the
earth conductivity remains in the potential-coefficient prefactor.
"""
function earth_potential_coefficient(
        ::Val{:Pollaczek1926}, ::Val{:underground}, functor, pair
)
    state = functor.state
    geometry = _geometry(pair)
    gamma_0 = state.gamma[2]
    direct = special_besselk(0, gamma_0 * geometry.d_ij) -
             special_besselk(0, gamma_0 * geometry.D_ij)
    kappa_1 = state.sigma[2] + state.jω * state.epsilon[2]
    return state.jω / (2π * kappa_1) * direct
end

raw"""
Evaluate Pollaczek's classical overhead-underground potential coefficient:

```math
P_{e,ij}^{01}=0,\qquad I_{ij}^{01}(\lambda)=0,\qquad A_{ij}^{01}=0.
```
"""
function earth_potential_coefficient(
        ::Val{:Pollaczek1926}, ::Val{:mixed}, functor, pair
)
    return zero(functor.state.jω)
end

:Pollaczek1926
