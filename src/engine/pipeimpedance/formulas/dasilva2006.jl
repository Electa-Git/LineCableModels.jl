function routes(identifier::Val{:DaSilva2006})
    return (
        self=FormulaMethod(identifier,pipe_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,pipe_impedance,Val(:mutual))
    )
end
assumptions(::Val{:DaSilva2006}) = (;)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Round solid cores at arbitrary offsets inside a finite-wall circular pipe; insulation parameters and bonding transformations are excluded. |
| Calculated quantities | Self ``Z(i,i)`` and mutual ``Z(j,i)`` impedances for proposed Method 3; finite-wall pipe skin and pipe-mediated coupling, and self core skin; no core-to-core eddy-current increment |
| Earth structure | Not applicable. |
| Model and approximation | Method 3 retains full finite-wall Bessel coefficients and omits core-to-core proximity. No expansion parameter or error bound is specified. |
| Main source | D. da Silva, G. Fernández, and R. A. Rivas (2006), Method 3. |
| Citation key(s) | `:DaSilva2006` |
| Evidence status | English equation pages checked; the Spanish pages corroborate the equations. |

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Round solid cores at arbitrary offsets inside an infinite-wall circular pipe; no outer pipe radius, insulation radius, or bonding network is included. |
| Calculated quantities | Self ``Z(i,i)`` and mutual ``Z(j,i)`` loop impedances for proposed Method 2; core skin and pairwise proximity, infinite-wall pipe surface and pipe-mediated coupling contributions |
| Earth structure | Not applicable; no earth boundary or layers in the extracted model. |
| Model and approximation | Method 2 combines core skin, core-to-core proximity, and infinite-wall pipe coupling. The positive-order series remain infinite, and the source gives no truncation error bound. |
| Main source | D. da Silva, G. Fernández, and R. A. Rivas (2006), Method 2. |
| Citation key(s) | `:DaSilva2006` |
| Evidence status | English equation pages checked; the Spanish pages corroborate them. The summation bound and root branch remain unresolved. |

**Expression.** The coefficient supplied to the shared pipe-interior block is

```math
H_{p,km}=Z_{pi}+\\frac{j\\omega\\mu_0}{2\\pi}
\\left[G_{km}+\\sum_{n=1}^{\\infty}a_n
\\Re\\left\\{\\left(\\frac{z_k\\overline z_m}{R_1^2}\\right)^n\\right\\}\\right],
```

where ``z_k=x_k+jy_k`` denotes an axis position relative to the pipe,
``G_{kk}=\\ln[(R_1^2-|z_k|^2)/(R_1r_k)]``, and
``G_{km}=\\ln[|R_1^2-z_k\\overline z_m|/(R_1|z_k-z_m|)]``.
The coefficients ``a_n`` are the finite-wall boundary response of
source equations (11)–(14), after extracting ``j\\omega\\mu_0/(2\\pi)``.

**Numerical conventions.** The positive pipe-return inner-surface response
``Z_{pi}`` is evaluated by Schelkunoff1934, also printed in
[Hoidalen2013](@cite), equation (8). This fixes the zeroth-order return-current
orientation; directly copying the opposite denominator sign of the 2006
``A_0,B_0`` pair would give a negative DC resistance. The source record
retains the published transcription.

The self core-skin term in source equation (7) is deliberately excluded from
``H_p``: it belongs to the individual coaxial block and must be added
once. Pipe outer-surface, through-wall transfer, jacket, and exterior-earth
terms also remain separate. Positive-order pipe harmonics are retained;
core-to-core eddy-current increments are not added to Method 3.

Scaled Bessel functions remove exponential growth. Wider arithmetic is used
when finite boundary ratios would otherwise overflow or underflow. The
`max_terms` limit reports failure rather than silently truncating a
nonconverged sum.

**Method 2.** `wall=:infinite` selects the infinite-wall source and defaults
to `proximity=:Kane1995`. It requires an infinite outer-radius argument
and physical core data through `with_cores`. The self sum includes all
other cores, as specified by the parent Kane (9), not an index-dependent
omission from the printed `P-1` upper bound. Outer-surface and transfer
requests are rejected: this selection supplies core-to-pipe coefficients,
not a finite-wall exterior assembly. The individual core skin term must
still be added once.

For finite walls, `proximity=:Hoidalen2013` adds the separate corrected
core increment to Method 3's pipe response. That combination is explicitly
labelled by both selections and is not attributed to Method 3 alone.

**Reference.** [DaSilva2006](@cite), Method 3, equations (5)–(14);
[Fortin2005](@cite), equations (6)–(7), for the equivalent finite-wall
vector-potential boundary problem;
[Hoidalen2013](@cite), equations (8)–(11), for the finite-wall surface and
harmonic restatement.
"""
description(::Formula{:DaSilva2006}) =
    "Da Silva et al. finite-pipe cavity without core-to-core proximity (2006)"

function (formula::Formula{:DaSilva2006})(
        radius::T,outer_radius::T,rho::T,mur::T,s::Complex{T}
) where {T <: Real}
    if formula.assumptions.wall===:infinite
        isinf(outer_radius) && outer_radius>0 ||
            throw(ArgumentError("infinite-wall coefficients require outer_radius=Inf and cannot supply a finite exterior matrix"))
        _validate_pipe(radius,2radius,rho,mur,s)
        piT=one(T)*π; mu0=T(4)*piT/T(10)^7
        x=sqrt(s*mu0*mur/rho)*radius
        prefactor=s*mu0/(2piT)
        inner=Complex{T}(prefactor*mur*special_besselkx(0,x)/(x*special_besselkx(1,x)))
        state=(;radius,outer_radius,rho,mur,s,x1=x,prefactor,inner,
            outer=nothing,transfer=nothing,wall=:infinite,
            proximity=formula.assumptions.proximity,
            rtol=max(T(formula.assumptions.rtol),eps(T)),max_terms=formula.assumptions.max_terms)
        return Functor{:DaSilva2006,typeof(formula.routes),typeof(state)}(formula.routes,state)
    end
    _validate_pipe(radius,outer_radius,rho,mur,s)
    πT=one(T)*π; μ0=T(4)*πT/T(10)^7
    m=sqrt(s*μ0*mur/rho)
    wall=InternalImpedance.Formula(:Schelkunoff1934)(radius,outer_radius,rho,mur,s)
    state=(;
        radius,outer_radius,rho,mur,s,x1=m*radius,x2=m*outer_radius,
        wall=:finite,proximity=formula.assumptions.proximity,
        prefactor=s*μ0/(2πT),inner=wall(Val(:inner)),
        outer=wall(Val(:outer)),transfer=wall(Val(:mutual)),
        rtol=max(T(formula.assumptions.rtol),eps(T)),
        max_terms=formula.assumptions.max_terms
    )
    return Functor{:DaSilva2006,typeof(formula.routes),typeof(state)}(formula.routes,state)
end

function pipe_impedance(::Val{:DaSilva2006},::Val{:self},functor,pair)
    pair.row==pair.column || throw(ArgumentError("pipe self route requires equal indices"))
    return pipe_impedance(Val(:DaSilva2006),Val(:cavity),functor,pair)
end

function pipe_impedance(::Val{:DaSilva2006},::Val{:mutual},functor,pair)
    pair.row!=pair.column || throw(ArgumentError("pipe mutual route requires distinct indices"))
    return pipe_impedance(Val(:DaSilva2006),Val(:cavity),functor,pair)
end

function pipe_impedance(::Val{:DaSilva2006},::Val{:cavity},functor,pair)
    state=functor.state
    geometry=_geometry(pair,state.radius)
    harmonic=state.wall===:infinite ? _harmonic_sum(
        n->_infinite_pipe_harmonic(n,state.x1,state.mur),state,geometry.w) :
        _harmonic_sum(state,geometry.w)
    return state.inner+state.prefactor*(geometry.geometric+harmonic)+_core_increment(functor,pair)
end

:DaSilva2006
