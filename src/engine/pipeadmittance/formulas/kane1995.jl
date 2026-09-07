function routes(identifier::Val{:Kane1995})
    return (
        self=FormulaMethod(identifier,pipe_potential_coefficient,Val(:self)),
        mutual=FormulaMethod(identifier,pipe_potential_coefficient,Val(:mutual))
    )
end
assumptions(::Val{:Kane1995})=(;)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Insulation admittance |
| Geometry | Core ``i`` radius ``a_i`` and offset ``b_i``; core ``j`` offset ``b_j``; angle ``\\theta_{i,j}``; cylindrical shield inner radius ``c_1``. A dielectric occupies the shield hollow. No multiple insulation shells, separate screen thicknesses, or layer-combination rule is printed here. |
| Calculated quantities | ``L_{ii}C_{ii}=\\mu_0\\varepsilon_0\\varepsilon_{ri}`` and ``P_{ji}=L_{ji}/(\\mu_0\\varepsilon_0\\varepsilon_{ri})``, with their full source-provided geometric inductance definitions |
| Earth structure | Not applicable. |
| Model and approximation | Lossless identities coupled to the paper's filament-based geometric inductances. No dielectric asymptotic expansion, discarded-loss term derivation, truncation order, or scalar-to-matrix construction is provided. “Lossless” is the author's case identification, not a reviewer-specified conductivity or loss-tangent value inserted into a more general expression. |
| Main source | Kane, Ahmad and Auriol (1995), “Multiwire Shielded Cable Parameter Computation,” IEEE Transactions on Magnetics 31(3), pp. 1646–1649; reproduced lossless/classical relations, not established original priority |
| Citation key(s) | `:Kane1995` |
| Evidence status | Original publication equations checked against PDF page images; dielectric terminology and multi-core assembly interpretation unresolved |

**Expression.** In a common homogeneous lossless dielectric, equation (8)
gives the self potential coefficient as the reciprocal of its scalar
core-to-shield capacitance. Equation (16) supplies the mutual coefficient:

```math
P_{kk}=\\frac{G_{kk}}{2\\pi\\varepsilon},
\\qquad P_{km}=\\frac{G_{km}}{2\\pi\\varepsilon}.
```

The geometric logarithms are exactly those used by the separate pipe
magnetic-field term. They are evaluated once by the shared circular geometry
function. This relation does not identify a scalar reciprocal with an
element of an inverse capacitance matrix.

**Assembly.** The coaxial engine selects one homogeneous cavity material.
It assembles these coefficients with local radial dielectric coefficients
and the outer pipe reference, then inverts the complete potential matrix.
A conducting-dielectric extension multiplies the lossless coefficient by
jωε/κ, using the independently selected material admittivity κ.
That extension and the matrix inversion are engine operations, not additional
equations attributed to Kane.

**Reference.** [Kane1995](@cite), equations (8), (15), and (16).
"""
description(::Formula{:Kane1995}) =
    "Kane et al. homogeneous pipe-interior potential coefficients (1995)"

function (formula::Formula{:Kane1995})(radius::T,epsilon::T) where {T <: Real}
    isfinite(radius) && radius>zero(T) && isfinite(epsilon) && epsilon>zero(T) ||
        throw(DomainError((radius,epsilon),"pipe radius and permittivity must be positive and finite"))
    return Functor{:Kane1995,typeof(formula.routes),T}(formula.routes,radius,epsilon)
end

function pipe_potential_coefficient(::Val{:Kane1995},::Val{:self},functor,pair)
    pair.row==pair.column || throw(ArgumentError("pipe self coefficient requires equal indices"))
    geometry=PipeImpedance._geometry(pair,functor.radius)
    return geometry.geometric/(2*(one(functor.radius)*π)*functor.epsilon)
end

function pipe_potential_coefficient(::Val{:Kane1995},::Val{:mutual},functor,pair)
    pair.row!=pair.column || throw(ArgumentError("pipe mutual coefficient requires distinct indices"))
    geometry=PipeImpedance._geometry(pair,functor.radius)
    return geometry.geometric/(2*(one(functor.radius)*π)*functor.epsilon)
end

:Kane1995
