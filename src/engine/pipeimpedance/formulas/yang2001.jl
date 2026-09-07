function routes(identifier::Val{:Yang2001})
    return (
        self=FormulaMethod(identifier,pipe_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,pipe_impedance,Val(:mutual))
    )
end
assumptions(::Val{:Yang2001}) = (;)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Multiple cores at arbitrary locations inside a circular enclosure. Circular pipe radii ``R_1,R_2`` and internal cable axes. |
| Calculated quantities | Finite-thickness pipe input, self, and mutual impedance components |
| Earth structure | External ground contribution ``Z_0`` appended separately. |
| Model and approximation | Hybrid finite-pipe harmonic expansion; numerical truncation of ``n`` is required. |
| Main source | Yixin Yang, J. Ma, and F. P. Dawalibi (2001) |
| Citation key(s) | `:Yang2001` |
| Evidence status | Original-page image verified |

**Expression.** The cavity coefficient combines the finite-wall inner-surface
impedance with infinite-wall positive-order harmonics:

```math
H_{p,km}=Z_{pi}+\\frac{j\\omega\\mu_0}{2\\pi}
\\left[G_{km}+\\sum_{n=1}^{\\infty}
\\frac{2\\mu_{rp}\\Re(w_{km}^{n})}
{n(1+\\mu_{rp})+x_1K_{n-1}(x_1)/K_n(x_1)}\\right],
\\qquad x_1=a\\sqrt{j\\omega\\mu_0\\mu_{rp}/\\rho}.
```

Here ``G_{km}`` is the circular cavity logarithm and
``w_{km}=z_k\\overline z_m/a^2``.
The finite-wall inner, outer, and transfer terms use Schelkunoff1934;
the outer and transfer terms are required for the unreduced pipe terminal
matrix, as restated in [DeSilva2019](@cite), (A3)–(A6).

**Numerical interpretation.** Equation (7)'s extra permeability symbol is
resolved by the dimensionally consistent finite-wall coefficient (A6) of
DeSilva2019. Dimensionless geometric harmonics carry one common
``j\\omega\\mu_0/(2\\pi)`` factor, not the repeated factor obtained by
substituting the printed (12) literally into (9) and (11).
The source record preserves those printed equations.

This hybrid is not the fully finite-wall eddy-current solution. The inner
radius enters the harmonic response, but the outer radius does not.
For thin magnetic walls at low frequency it therefore differs from
DaSilva2006. No wall-thickness error bound is supplied by the source.
Individual core skin and core proximity are not part of this cavity term.

The equivalent appendix record is
[DeSilva2019](@cite).
Its equation (A2) contains ``K_{n-1}/K_n``, not ``I_{n-1}/K_n``.
The source-record transcription has been corrected against the page image.

**Reference.** [Yang2001](@cite), equations (7)–(12);
[DeSilva2019](@cite), equations (A1)–(A6).
"""
description(::Formula{:Yang2001}) =
    "Yang finite-wall surfaces with infinite-wall pipe harmonics (2001)"

function (formula::Formula{:Yang2001})(
        radius::T,outer_radius::T,rho::T,mur::T,s::Complex{T}
) where {T <: Real}
    full=Formula(:DaSilva2006;
        rtol=formula.assumptions.rtol,max_terms=formula.assumptions.max_terms,
        proximity=formula.assumptions.proximity
    )(radius,outer_radius,rho,mur,s)
    state=full.state
    return Functor{:Yang2001,typeof(formula.routes),typeof(state)}(formula.routes,state)
end

function pipe_impedance(::Val{:Yang2001},::Val{:self},functor,pair)
    pair.row==pair.column || throw(ArgumentError("pipe self route requires equal indices"))
    return pipe_impedance(Val(:Yang2001),Val(:cavity),functor,pair)
end

function pipe_impedance(::Val{:Yang2001},::Val{:mutual},functor,pair)
    pair.row!=pair.column || throw(ArgumentError("pipe mutual route requires distinct indices"))
    return pipe_impedance(Val(:Yang2001),Val(:cavity),functor,pair)
end

function pipe_impedance(::Val{:Yang2001},::Val{:cavity},functor,pair)
    state=functor.state
    geometry=_geometry(pair,state.radius)
    harmonic=_harmonic_sum(state,geometry.w) do n
        _infinite_pipe_harmonic(n,state.x1,state.mur)
    end
    return state.inner+state.prefactor*(geometry.geometric+harmonic)+_core_increment(functor,pair)
end

:Yang2001
