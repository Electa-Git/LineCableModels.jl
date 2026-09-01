function routes(identifier::Val{:Pettersson1994})
    return (
        self = formula_method(identifier, earth_potential_coefficient, Val(:self)),
        mutual = formula_method(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = formula_method(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Pettersson1994})
    (
        air = _full,
        earth = _full,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Pettersson1994}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Wideband complex-image potential coefficient for overhead
conductors above homogeneous earth. Pettersson's paper treats wires above, on,
and under ground; the registered route implements only its overhead equation.

**Expression.**

```math
P_{e,ij}=\\frac{1}{2\\pi\\varepsilon_0}\\ln\\frac{D_{ij}}{d_{ij}}+N_{e,ij},
```

```math
N_{e,ij}=\\frac{1}{(n^2+1)\\pi\\varepsilon_0}
\\ln\\frac{\\sqrt{[H+(n^2+1)/\\beta_\\gamma]^2+y_{ij}^2}}{D_{ij}},
\\quad n^2=\\varepsilon_{rg}+\\frac{\\sigma_g}{j\\omega\\varepsilon_0}.
```

**Reference.** P. Pettersson, “Image Representation of Wave Propagation on
Wires Above, On and Under Ground,” *IEEE Transactions on Power Delivery*, 9,
1049–1055, 1994. DOI: 10.1109/61.296290.
"""
function description(::Formula{:Pettersson1994})
    "Pettersson wideband overhead image potential coefficient (1994)"
end

function propagation_constant(
        ::Val{:Pettersson1994}, jω, permeability, permittivity
)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Pettersson1994})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Pettersson1994), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate Pettersson's wideband overhead earth potential coefficient:

```math
P_{e,ij}=\frac{1}{2\pi\varepsilon_0}\ln\frac{D_{ij}}{d_{ij}}
+N_{e,ij},
```

```math
N_{e,ij}=\frac{1}{(n^2+1)\pi\varepsilon_0}
\ln\frac{\sqrt{[H+(n^2+1)/\beta_\gamma]^2+y_{ij}^2}}{D_{ij}},
```

```math
n^2=\varepsilon_{rg}+\frac{\sigma_g}{j\omega\varepsilon_0},
\qquad
\beta_\gamma=\sqrt{\gamma_g^2-\gamma_0^2}.
```

The 2020 secondary transcription omits the perfect-ground coefficient
``1/(2\pi\varepsilon_0)`` in its admittance equation and a square in the
radicand of ``N_{e,ij}``; both are restored here so that the terms have the
published potential-coefficient dimensions and the image distance has units
of length.

# Reference

P. Pettersson, "Image representation of wave propagation on wires above,
on and under ground," *IEEE Transactions on Power Delivery*, vol. 9,
pp. 1049-1055, 1994. DOI: 10.1109/61.296290.
"""
function earth_potential_coefficient(
        ::Val{:Pettersson1994}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    γ0_squared, γg_squared = state.gamma_medium_squared
    βγ = sqrt(γg_squared - γ0_squared)
    n_squared = state.epsilon[2] / state.epsilon[1] +
                state.sigma[2] / (state.jω * state.epsilon[1])
    image = sqrt(
        (geometry.H + (n_squared + 1) / βγ)^2 + geometry.y_ij^2
    )
    πT = one(geometry.H) * π
    perfect = log(geometry.D_ij / geometry.d_ij) /
              (2πT * state.epsilon[1])
    correction = log(image / geometry.D_ij) /
                 ((n_squared + 1) * πT * state.epsilon[1])
    return perfect + correction
end

:Pettersson1994
