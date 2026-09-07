function routes(identifier::Val{:Pettersson1994})
    return (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
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

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Bare round thin wire represented by a line source; no coating. |
| Calculated quantities | Generalized p.u.l. self/mutual shunt admittance of a thin wire above, at, or below a lossy planar interface |
| Earth structure | Two homogeneous half-spaces with a flat boundary. |
| Model and approximation | Uses the same integrand replacement (9) as the series record, but ``b=n^2`` for ``Q``. On the interface the source uses (12)–(15), not the ``h>0`` distances above. |
| Main source | Pär Pettersson (1994 publication; 1993 conference manuscript) |
| Citation key(s) | `:Pettersson1994` |
| Evidence status | Original publication page images checked |

**Expression.**

```math
P_{e,ij}=\\frac{1}{2\\pi\\varepsilon_0}\\ln\\frac{D_{ij}}{d_{ij}}+N_{e,ij},
```

```math
N_{e,ij}=\\frac{1}{(n^2+1)\\pi\\varepsilon_0}
\\ln\\frac{\\sqrt{[H+(n^2+1)/\\beta_\\gamma]^2+y_{ij}^2}}{D_{ij}},
\\quad n^2=\\varepsilon_{rg}+\\frac{\\sigma_g}{j\\omega\\varepsilon_0}.
```

**Numerical interpretation.** The evaluator retains conduction and displacement
current in both media, uses the source-medium factor
``j\\omega/(\\sigma_1+j\\omega\\varepsilon_1)``, and selects the
printed sign of the Q image. Equations (10)–(11) apply separately above or
below ground; (14)–(15) apply to wires on the interface. Mixed and
interface/off-interface pairs require a common modal prescription and are
not provided by these closed images. The returned coefficient populates
the full potential matrix before admittance conversion.

**Reference.** P. Pettersson, “Image Representation of Wave Propagation on
Wires Above, On and Under Ground,” *IEEE Transactions on Power Delivery*, 9,
1049–1055, 1994. DOI: 10.1109/61.296290.
"""
function description(::Formula{:Pettersson1994})
    "Pettersson two-half-space and interface image potential coefficient (1994)"
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
    state=functor.state
    images=pettersson_images(state,pair)
    return _complex_result(state.jω,state.jω/
        (2*(one(pair.separation)*π)*images.kappa)*images.electric)
end

:Pettersson1994
