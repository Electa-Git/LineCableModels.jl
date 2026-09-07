function routes(identifier::Val{:Dubanton1969})
    return (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

assumptions(::Val{:Dubanton1969}) = (
    air = _lossless, earth = _conductive, permeability = vacuum_permeability
)
propagation(::Val{:Dubanton1969}) = Val(:zero)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Ideal, infinitely long parallel thin wires represented by radius, heights, and separation; conductor internal impedance and insulation are outside the expression. |
| Calculated quantities | Per-unit-length self and mutual ideal-conductor/ground-return loop impedances for overhead wires |
| Earth structure | Homogeneous conductive half-space under a plane interface. |
| Model and approximation | The authors transform Carson's correction integral and introduce the explicit kernel approximation in their equation (37), then obtain (3) and (4). Thus the closed logarithms are approximations to Carson, not exact full-wave expressions. Their proof of the mutual result initially requires ``β<1``; the wider practical range is supported by numerical testing. |
| Main source | Equations originally proposed by Dubanton and published by Gary; heuristic justification, analytical relation to Carson, and numerical error evaluation by A. Déri, G. Tevan, A. Semlyen, and A. Castanheira (1981) |
| Citation key(s) | Primary attribution: `:Dubanton1969`; equation-complete English derivation: `:Deri1981`; auxiliary complex-image derivation: `:Wait1969` |
| Evidence status | PDF page images checked for the English derivation; DOI agrees with its primary-PDF metadata |

**Expression.** With ``p=1/\\sqrt{j\\omega\\mu_0\\sigma}``,

```math
Z_{e,ii}=\\frac{j\\omega\\mu_0}{2\\pi}\\ln\\frac{2(h_i+p)}{r_i},
```

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}
\\ln\\frac{\\sqrt{(h_i+h_j+2p)^2+x_{ij}^2}}
{\\sqrt{(h_i-h_j)^2+x_{ij}^2}}.
```

**Reference.** [Dubanton1969](@cite); English derivation and error analysis:
[Deri1981](@cite), equations (3), (4), and (18). The former selections
`:Gary1976` and `:DeriSemlyen1981` resolve to this formula.
"""
description(::Formula{:Dubanton1969}) =
    "Dubanton complex ground-return-plane approximation (1969)"

function propagation_constant(::Val{:Dubanton1969}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Dubanton1969})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Dubanton1969), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

function earth_impedance(::Val{:Dubanton1969}, ::Val{:self}, functor, pair)
    return _image_plane_impedance(functor, pair, inv(functor.state.gamma[2]), Val(:self))
end

function earth_impedance(::Val{:Dubanton1969}, ::Val{:mutual}, functor, pair)
    return _image_plane_impedance(functor, pair, inv(functor.state.gamma[2]), Val(:mutual))
end

:Dubanton1969
