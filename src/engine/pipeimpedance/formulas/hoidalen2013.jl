function routes(identifier::Val{:Hoidalen2013})
    return (
        self=FormulaMethod(identifier,pipe_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,pipe_impedance,Val(:mutual))
    )
end
assumptions(::Val{:Hoidalen2013}) = (;)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Finite circular pipe and equal-core geometry using radii ``r_{p1},r_{p2},r_{p3},r_{1k},r_{2k}``, offsets ``d_k``, and separations ``d_{km}``. |
| Calculated quantities | Low-frequency pipe inner-surface impedance, pipe-mediated core interaction, surface-connection combination, and source-provided cable-loop/mode limits |
| Earth structure | Not applicable to pipe/internal and mode-cancelled terms. Not stated for the unspecified ground-return input in the total-loop formula. |
| Model and approximation | Low-frequency limits of the finite-pipe Bessel model combined through the printed cable and modal relations. Equation (35) also assumes a thin wall; the source gives no discarded order or error remainder. |
| Main source | Høidalen's low-frequency analysis of the finite-pipe formulas attributed to Kane et al. (1995), da Silva et al. (2006), and earlier cable assembly |
| Citation key(s) | `:Hoidalen2013` |
| Evidence status | Author manuscript equations checked; printed limit notation documented. |

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Infinite-wall parent has finite inner radius ``r_{p1}``, no outer radius in (14)–(16). Core offsets ``d_k,d_m`` and angle ``\\theta_{km}`` define coupling. The separately labelled mixed-model diagnostic reuses finite-wall outer/transfer terms and pipe DC resistance. |
| Calculated quantities | Low-frequency logarithmic inner-pipe surface term; pipe-mediated core-coupling correction; author's incomplete diagnostic expression for mixed finite/infinite pipe assembly |
| Earth structure | Not applicable. |
| Model and approximation | Høidalen takes the low-frequency behavior of the infinite-wall Bessel parent (2)–(3), identifying its ``K_0`` logarithmic behavior and retaining the additional ``C_1`` term (16). The source does not supply a complete term-by-term expansion, remainder, or a uniform finite-thickness limit. Taking a low-frequency limit after an infinite-wall assumption is not equated here with taking the low-frequency limit at fixed physical wall thickness. |
| Main source | Høidalen's low-frequency analysis of the infinite-wall model attributed to Brown–Rocamora (1976) and Ametani (1980) |
| Citation key(s) | `:Hoidalen2013` |
| Evidence status | Manuscript page images checked for (14)–(16) and their dependencies; the source leaves diagnostic ``X`` unexpanded. |

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Round solid cores separated by ``d_{km}`` in the tested pipe-cable geometry. |
| Calculated quantities | Core-to-core proximity contribution to per-length series impedance; published differential/common-mode multipliers for the symmetrical three-core case |
| Earth structure | Not applicable; the correction is independent of earth layering or burial depth. |
| Model and approximation | The source subtracts the low-frequency inductive term from Kane's expression and adds the ``k=1`` Dwight term. The sum over ``n`` remains infinite, and no discarded-series error bound is supplied. |
| Main source | Hans Kr. Høidalen, “Analysis of Pipe-Type Cable Impedance Formulations at Low Frequencies,” 2013, DOI `10.1109/TPWRD.2013.2272343` |
| Citation key(s) | `:Hoidalen2013` |
| Evidence status | Equation (36) and its dependencies checked. |

**Expression.** The low-frequency surface terms are
``Z_{pi}=R_{dc}+j\\omega L_i``,
``Z_{po}=R_{dc}+j\\omega L_o``, and
``Z_{pm}=R_{dc}+j\\omega L_t``. With ``q=\\ln(b/a)``
and ``c=\\mu_0\\mu_{rp}/(2\\pi)``,

```math
\\begin{aligned}
L_t/c&=(q/\\sinh^2 q-\\coth q)/4, \\\\
L_i/c&=L_t/c+q/2+(q\\coth q-1)/2, \\\\
L_o/c&=L_t/c+q/2-(q\\coth q-1)/2.
\\end{aligned}
```

These expressions give equation (12) and the low-frequency expansion of
(20)–(21), with ``L_i+L_o-2L_t=cq`` as in (22).
The DC resistance is ``R_{dc}=\\rho/[\\pi(b-a)(b+a)]``.

The cavity coefficient uses the same circular geometric term as DaSilva2006,
with the finite-wall static harmonic response derived from (9)–(11):

```math
a_{n,0}=\\frac{2\\mu_{rp}}{n(1+\\mu_{rp})}
\\frac{1-r(a/b)^{2n}}{1-r^2(a/b)^{2n}},
\\qquad r=\\frac{\\mu_{rp}-1}{\\mu_{rp}+1}.
```

**Numerical interpretation.** At ``\\mu_{rp}=1``, the harmonic sum is
``-\\ln|1-w|`` and agrees with the printed equation (13).
For a magnetic finite wall, the thickness-dependent expression above is
the limit of the full Bessel parent. Equation (13) instead gives its
infinite-outer-radius limit; that expression is not substituted for the
finite wall. The printed equation is retained in the source record.

This is an explicitly selected low-frequency approximation, not an automatic
switch from the full skin-effect solution. The source supplies no numerical
frequency threshold. Core skin, core-to-core proximity, external insulation,
and earth return remain separate. Source total-loop and symmetric-mode
relations follow only under their stated geometry and material restrictions.
Thin-wall surface combinations use a local series to avoid cancellation.

**Infinite-wall selection.** `wall=:infinite` evaluates (14) and the
low-frequency harmonic expansion of (3), with an infinite outer-radius
argument. It has no finite outer or transfer surface. With
``L=\\ln(2/x_1)-\\gamma``, the first positive-order parent coefficient is
``2\\mu_{rp}/(1+\\mu_{rp}+x_1^2L)``. It is the static coefficient minus
the positive correction printed in (16). The numerical evaluation therefore
subtracts that correction, while the survey preserves the plus sign in
printed (15). This is a low-frequency approximation within an already
infinite wall; no finite-wall DC limit is implied.

**Core proximity.** `core_proximity(Val(:Hoidalen2013), ...)` supplies (36)
independently of the wall approximation. A pipe selection with
`proximity=:Hoidalen2013` adds it through the physical core data. Matrix
assembly is restricted to the source's symmetrical three equal
nonmagnetic cores: the increment has diagonal entries twice the pair
term and off-diagonal entries equal to that term. Its eigenvalues are
the published differential and common factors one and four.
The static logarithm is subtracted analytically term by term, so the
low-frequency difference is not lost by subtracting nearly equal sums.

**Reference.** [Hoidalen2013](@cite), equations (8)–(13), (20)–(22).
"""
description(::Formula{:Hoidalen2013}) =
    "Høidalen finite-pipe low-frequency surface and cavity terms (2013)"

function (formula::Formula{:Hoidalen2013})(
        radius::T,outer_radius::T,rho::T,mur::T,s::Complex{T}
) where {T <: Real}
    if formula.assumptions.wall===:infinite
        isinf(outer_radius) && outer_radius>0 ||
            throw(ArgumentError("infinite-wall coefficients require outer_radius=Inf and cannot supply a finite exterior matrix"))
        _validate_pipe(radius,2radius,rho,mur,s)
        piT=one(T)*π; mu0=T(4)*piT/T(10)^7
        x=sqrt(s*mu0*mur/rho)*radius
        logarithm=log(2/x)-T(Base.MathConstants.eulergamma)
        prefactor=s*mu0/(2piT)
        state=(;radius,outer_radius,rho,mur,s,x1=x,logarithm,prefactor,
            inner=prefactor*mur*logarithm,outer=nothing,transfer=nothing,wall=:infinite,
            proximity=formula.assumptions.proximity,
            rtol=max(T(formula.assumptions.rtol),eps(T)),max_terms=formula.assumptions.max_terms)
        return Functor{:Hoidalen2013,typeof(formula.routes),typeof(state)}(formula.routes,state)
    end
    _validate_pipe(radius,outer_radius,rho,mur,s)
    πT=one(T)*π; μ0=T(4)*πT/T(10)^7
    q=log1p((outer_radius-radius)/radius)
    if T <: Union{Float32,Float64} && q<T(1)/100
        q2=q*q
        lt=q*(-one(T)/6+q2*(one(T)/45+q2*(-one(T)/315+q2*2/4725)))
        half_difference=q2*(one(T)/6+q2*(-one(T)/90+q2*(one(T)/945-q2/9450)))
    else
        cot=inv(tanh(q))
        lt=(q/sinh(q)^2-cot)/4
        half_difference=(q*cot-one(T))/2
    end
    mean=lt+q/2
    scale=μ0*mur/(2πT)
    dc=rho/(πT*(outer_radius-radius)*(outer_radius+radius))
    state=(;
        radius,outer_radius,rho,mur,s,log_ratio=q,
        wall=:finite,proximity=formula.assumptions.proximity,
        prefactor=s*μ0/(2πT),
        inner=dc+s*scale*(mean+half_difference),
        outer=dc+s*scale*(mean-half_difference),
        transfer=dc+s*scale*lt,
        rtol=max(T(formula.assumptions.rtol),eps(T)),
        max_terms=formula.assumptions.max_terms
    )
    return Functor{:Hoidalen2013,typeof(formula.routes),typeof(state)}(formula.routes,state)
end

function pipe_impedance(::Val{:Hoidalen2013},::Val{:self},functor,pair)
    pair.row==pair.column || throw(ArgumentError("pipe self route requires equal indices"))
    return pipe_impedance(Val(:Hoidalen2013),Val(:cavity),functor,pair)
end

function pipe_impedance(::Val{:Hoidalen2013},::Val{:mutual},functor,pair)
    pair.row!=pair.column || throw(ArgumentError("pipe mutual route requires distinct indices"))
    return pipe_impedance(Val(:Hoidalen2013),Val(:cavity),functor,pair)
end

function pipe_impedance(::Val{:Hoidalen2013},::Val{:cavity},functor,pair)
    state=functor.state
    geometry=_geometry(pair,state.radius)
    if state.wall===:infinite
        # From parent (3), the first harmonic is smaller than its static value.
        # Thus the positive Delta printed in (16) is subtracted from (13).
        scale=2state.mur/(1+state.mur)
        correction=scale*real(geometry.w)*state.logarithm*state.x1^2/
            (1+state.mur+state.logarithm*state.x1^2)
        harmonic=-scale*log(abs(1-geometry.w))-correction
        return state.inner+state.prefactor*(geometry.geometric+harmonic)+_core_increment(functor,pair)
    end
    # Nonmagnetic static harmonics cancel the geometric image logarithm.
    # Evaluate the resulting free-space logarithm without a truncated sum.
    if isone(state.mur)
        distance=pair.row==pair.column ? pair.radii[1] :
            hypot(pair.positions[1][1]-pair.positions[2][1],
                  pair.positions[1][2]-pair.positions[2][2])
        field=log(state.radius/distance)
    else
        field=geometry.geometric+_harmonic_sum(state,geometry.w) do n
            _static_pipe_harmonic(n,state.log_ratio,state.mur)
        end
    end
    return state.inner+state.prefactor*field+_core_increment(functor,pair)
end

:Hoidalen2013
