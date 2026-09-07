function routes(identifier::Val{:Kane1995})
    return (
        self=FormulaMethod(identifier,pipe_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,pipe_impedance,Val(:mutual)))
end
assumptions(::Val{:Kane1995})=(;)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | ``P`` round cores at arbitrary offsets inside a finite cylindrical shield; conductor insulation and bonding are excluded. |
| Calculated quantities | Filament–shield self and mutual impedances; inner shield surface impedance; solid-core skin term; other-core eddy-current contributions; complete self and mutual core–shield loop assembly |
| Earth structure | Not applicable. |
| Model and approximation | The model combines filament–shield fields, solid-core skin impedance, and selected pairwise proximity terms. The source gives no interaction order, discarded-term list, error bound, or recurrence equations. |
| Main source | M. Kane, A. Ahmad, and P. Auriol (1995). |
| Citation key(s) | `:Kane1995` |
| Evidence status | Original publication equations checked; excitation, definition, and citation conflicts remain unresolved. |

**Numerical scope.** The finite-wall response is shared with DaSilva2006.
The default adds Kane's core-to-core series through physical core data
supplied with `with_cores`. Core skin belongs to the individual coaxial
block and is not added again here.

**Expression.** The core contribution uses

```math
\\Delta Z_k(d)=\\frac{j\\omega\\mu_0}{2\\pi}
\\sum_{n=1}^{\\infty}
\\frac{2\\mu_{rk}(a_k/d)^{2n}}
{n(\\mu_{rk}-1)+z_k I_{n-1}(z_k)/I_n(z_k)},
\\qquad z_k=a_k\\sqrt{j\\omega\\mu_0\\mu_{rk}\\sigma_k}.
```

The self increment sums all other cores. The source's mutual increment
uses the receiving core. Finite-wall and geometric terms are the
shared pipe coefficients, not another copy of the core-skin response.

The reciprocal matrix route requires equal solid-core radii,
resistivities, and permeabilities. The separate `core_proximity` leaf
retains the source's directional unequal-core expression; no averaging
of its unequal mutual coefficients is performed.

The same pair term is restated by [DaSilva2006](@cite), (4), and
[Hoidalen2013](@cite), (25). It does not solve multiple scattering.

**Reference.** [Kane1995](@cite), (4)–(6),(9)–(11),(14),(17)–(18).
"""
description(::Formula{:Kane1995}) =
    "Kane finite-pipe response with pairwise solid-core proximity (1995)"

function (formula::Formula{:Kane1995})(radius::T,outer_radius::T,
        rho::T,mur::T,s::Complex{T}) where {T <: Real}
    parent=Formula(:DaSilva2006;rtol=formula.assumptions.rtol,
        max_terms=formula.assumptions.max_terms,proximity=formula.assumptions.proximity)
    full=parent(radius,outer_radius,rho,mur,s)
    return Functor{:Kane1995,typeof(formula.routes),typeof(full.state)}(formula.routes,full.state)
end

function pipe_impedance(::Val{:Kane1995},::Val{:self},functor,pair)
    pair.row==pair.column || throw(ArgumentError("pipe self route requires equal indices"))
    return pipe_impedance(Val(:DaSilva2006),Val(:cavity),functor,pair)
end

function pipe_impedance(::Val{:Kane1995},::Val{:mutual},functor,pair)
    pair.row!=pair.column || throw(ArgumentError("pipe mutual route requires distinct indices"))
    return pipe_impedance(Val(:DaSilva2006),Val(:cavity),functor,pair)
end

:Kane1995
