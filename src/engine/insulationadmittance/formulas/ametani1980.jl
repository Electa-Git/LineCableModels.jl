assumptions(::Val{:Ametani1980}) = (;)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Insulation admittance |
| Geometry | Concentric annuli ``r_2<r<r_3``, ``r_4<r<r_5``, and ``r_6<r<r_7``. |
| Calculated quantities | Core, sheath, and armor insulation potential coefficients and admittance-matrix assembly |
| Earth structure | Not applicable. |
| Model and approximation | Not an analytical approximation within the concentric, lossless dielectric model. Omission of dielectric loss is a constitutive restriction, not a corpus modification. |
| Main source | A. Ametani (1980) |
| Citation key(s) | `:Ametani1980` |
| Evidence status | PDF page images checked |

**Expression.** The lossless material relation and radial potential coefficient are

```math
\\kappa=j\\omega\\varepsilon_0\\varepsilon_r,
\\qquad
p=\\frac{\\ln(b/a)}{2\\pi\\varepsilon_0\\varepsilon_r}.
```

The shared radial assembler forms the source's nested potential matrix;
admittance is obtained by inverting that complete matrix. Finite material
conductivity is omitted by this selection, as required by the source model.

**Reference.** [Ametani1980](@cite), equations (4), (21)–(24).
The former Gustavsen2013 selection described this same lossless relation and
remains a constructor alias.
"""
description(::Formula{:Ametani1980}) =
    "Ametani lossless coaxial-insulation potential coefficients (1980)"

@inline function insulation_material(
        ::Val{:Ametani1980}, material::Material{T},
        frequency::T, temperature::T, values::NamedTuple
) where {T <: Real}
    ε0=one(T)*88541878128*(one(T)*10)^(-22)
    ω=2*(one(T)*π)*frequency
    return complex(zero(T),ω)*ε0*material.eps_r
end

:Ametani1980
