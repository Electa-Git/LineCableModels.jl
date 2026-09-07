assumptions(::Val{:Weeks1984}) = (;loss_tangent=0.0)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Insulation admittance |
| Geometry | Single coaxial core-to-sheath radial path. Inner semiconducting screen, main insulation, outer semiconducting screen. |
| Calculated quantities | Core-to-sheath p.u.l. admittance of inner screen, lossy insulation, and outer screen in series |
| Earth structure | None. |
| Model and approximation | The series equation is the coaxial layer model; the ``y\\simeq y_2`` reduction assumes both screen admittances greatly exceed the insulation admittance. |
| Main source | W. Weeks and Yi Diao (1984) |
| Citation key(s) | `:Weeks1984` |
| Evidence status | Original publication page images checked |

**Expression.** The selected main-insulation relation is

```math
\\kappa=\\omega\\varepsilon_0\\varepsilon_r D_f+
j\\omega\\varepsilon_0\\varepsilon_r,
\\qquad
y_2=\\frac{2\\pi\\kappa}{\\ln(a_2/a_1)}.
```

Set `loss_tangent` to the nonnegative dissipation factor `D_f`.
It represents the complete dielectric loss in this expression; material
DC conductivity is not added again. The default zero factor gives the
lossless limit. Semiconducting screens use their separately selected
constitutive relation. Their potential coefficients are added to the
insulation coefficient before inversion, giving
`inv(y) = inv(y1) + inv(y2) + inv(y3)`.

**Reference.** [Weeks1984](@cite), equations (9)–(10).
"""
description(::Formula{:Weeks1984}) =
    "Weeks–Diao loss-tangent insulation and radial series admittance (1984)"

@inline function insulation_material(
        ::Val{:Weeks1984}, material::Material{T},
        frequency::T, temperature::T, values::NamedTuple
) where {T <: Real}
    factor=values.loss_tangent
    factor isa Real && isfinite(factor) && factor>=zero(factor) ||
        throw(DomainError(factor,"dielectric loss tangent must be nonnegative and finite"))
    ε0=one(T)*88541878128*(one(T)*10)^(-22)
    ω=2*(one(T)*π)*frequency
    return complex(T(factor),one(T))*ω*ε0*material.eps_r
end

:Weeks1984
