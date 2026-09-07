function routes(identifier::Val{:Uribe2008})
    return (
        self=FormulaMethod(identifier,earth_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,earth_impedance,Val(:mutual)),
        overhead=FormulaMethod(identifier,earth_impedance,Val(:overhead)),
        underground=FormulaMethod(identifier,earth_impedance,Val(:underground)),
        mixed=FormulaMethod(identifier,earth_impedance,Val(:mixed)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:Uribe2008})=(
    air=_lossless,earth=_conductive,permeability=vacuum_permeability,
    approximation=:ccitt
)
propagation(::Val{:Uribe2008})=Val(:zero)

function Formula(::Val{:Uribe2008};approximation::Symbol=:ccitt,kwargs...)
    approximation in (:ccitt,:wedepohl) ||
        throw(ArgumentError(":Uribe2008 approximation must be :ccitt or :wedepohl"))
    defaults=routes(Val(:Uribe2008)); overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) ||
        throw(ArgumentError("unknown routes for earth-impedance formula :Uribe2008"))
    selected=merge(defaults,overrides)
    values=merge(assumptions(Val(:Uribe2008)),(;approximation))
    return Formula{:Uribe2008,typeof(selected),typeof(values)}(selected,values)
end

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel overhead and buried line axes. |
| Calculated quantities | CCITT-recommended mixed mutual ground-impedance approximation and Uribe's normalized form |
| Earth structure | Homogeneous conductive soil half-space below air. |
| Model and approximation | Source-labelled recommended approximation; derivation, retained order, and universal range unresolved. |
| Main source | CCITT via G. Lucca and H. W. Dommel–J. Sawada, secondary attribution through F. A. Uribe (2008) |
| Citation key(s) | `:Uribe2008` |
| Evidence status | Secondary transcription checked against page images against Uribe; credited recommendation/report not equation-verified |

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel overhead and buried line axes. |
| Calculated quantities | Wedepohl-attributed low-frequency mixed approximation and Uribe's normalized form |
| Earth structure | Homogeneous conductive soil half-space below air. |
| Model and approximation | Low-frequency retained-term approximation of the Wedepohl–Wilcox series; exact retained order and unambiguous frequency bound unresolved. |
| Main source | L. M. Wedepohl and D. J. Wilcox (1973), secondary mixed-form attribution through F. A. Uribe (2008) |
| Citation key(s) | `:Uribe2008` |
| Evidence status | Secondary mixed-form transcription checked against page images against Uribe; exact original mixed locator unresolved |

**Expression.** The default `approximation=:ccitt` evaluates the
dimensional source equation (6c). With
``g=\\sqrt{j\\omega\\mu_0\\sigma_g}``, positive overhead height
``h_a``, positive burial depth ``h_g``, and direct distance
``d=\\sqrt{x^2+(h_a+h_g)^2}``, it is

```math
Z_C=\\frac{j\\omega\\mu_0}{2\\pi}
\\left[\\ln\\frac{1.851}{gd}+\\frac23g(h_a-h_g)\\right].
```

The `approximation=:wedepohl` selection evaluates the dimensional
structure of (6e),

```math
Z_W=\\frac{j\\omega\\mu_0}{2\\pi}
\\left[-\\ln\\frac{e^{\\gamma_E}gd}{2}
+\\frac12-\\frac43g(h_a+h_g)\\right].
```

The multiplicative Euler factor is ``e^{\\gamma_E}\\simeq1.78107``,
as fixed by the small-argument Bessel normalization of the credited
Wedepohl–Wilcox parent. The coefficient ``4/3`` is retained from
Uribe's dimensional (6e), not replaced by the ``2/3`` in (6f).
Consequently this mixed restatement is not registered as another copy
of the 1973 underground mutual expression. It reuses that logarithmic
evaluator with twice the height-sum argument.

**Assembly and limits.** Only the mixed mutual coefficient changes.
Overhead and underground pairs use the existing exact Pollaczek/Carson
leaves. Coordinate signs enter the CCITT height sum; interchanging
source and observation preserves reciprocity. The square root has positive
real part, with negative-frequency conjugation.

Both alternatives are conductivity-only approximations in two
nonmagnetic half-spaces. Source-normalized (6d) and (6f) are not treated
as extra formulas or used to rescale the dimensional coefficients.
The source supplies no universal error bound; these selections must
not be read as uniformly accurate substitutes for the mixed integral.

**Reference.** [Uribe2008](@cite), equations (6c)–(6f), as the
secondary mixed-form source. [Wedepohl1973](@cite) supplies the
Bessel-series logarithmic normalization, not independent verification
of Uribe's mixed height coefficient.
"""
description(::Formula{:Uribe2008}) =
    "Uribe's CCITT and Wedepohl mixed approximations (2008)"

function propagation_constant(::Val{:Uribe2008},s,mu,epsilon)
    return (Γ=zero(s),squared=zero(s))
end

function (formula::Formula{:Uribe2008})(rho,epsilon,mu,s,Γ,segments=nothing)
    _check(rho,epsilon,mu)
    length(rho)==2 || throw(ArgumentError(":Uribe2008 requires one earth half-space"))
    isinf(rho[1]) && isfinite(rho[2]) && rho[2]>0 ||
        throw(DomainError(rho,":Uribe2008 requires lossless air and conducting earth"))
    all(x->isapprox(x,vacuum_permeability(x)),mu) ||
        throw(DomainError(mu,":Uribe2008 assumes nonmagnetic media"))
    isfinite(s) && !iszero(s) && iszero(real(s)) ||
        throw(DomainError(s,":Uribe2008 requires nonzero real frequency"))
    return _homogeneous_functor(Val(:Uribe2008),formula,rho,epsilon,mu,s,Γ,segments)
end

function earth_impedance(::Val{:Uribe2008},::Val{:mutual},functor,pair)
    placement=_placement(pair)
    placement===Val(:overhead) && return functor.routes.overhead(functor,pair)
    placement===Val(:underground) && return functor.routes.underground(functor,pair)
    return functor.routes.mixed(functor,pair)
end

function earth_impedance(::Val{:Uribe2008},::Val{:overhead},functor,pair)
    return earth_impedance(Val(:Carson1926),Val(:mutual),functor,pair)
end

function earth_impedance(::Val{:Uribe2008},::Val{:underground},functor,pair)
    return earth_impedance(Val(:Pollaczek1926),Val(:underground),functor,pair)
end

function earth_impedance(::Val{:Uribe2008},::Val{:mixed},functor,pair)
    _require(pair,Val(:mixed))
    state=functor.state
    air=pair.layers[1]==1 ? 1 : 2
    soil=air==1 ? 2 : 1
    ha=abs(pair.heights[air]); hg=abs(pair.heights[soil])
    d=hypot(pair.separation,ha+hg)
    all(isfinite,(ha,hg,d)) && ha>0 && hg>0 && d>0 ||
        throw(DomainError((ha,hg,d),":Uribe2008 requires positive height and burial depth"))
    if state.formula.assumptions.approximation===:wedepohl
        return earth_impedance(Val(:Wedepohl1973),Val(:kernel),functor,d,2(ha+hg))
    end
    g=state.gamma[2]; T=typeof(ha)
    bracket=log(T(1851)/(T(1000)*g*d))+2g*(ha-hg)/3
    return _complex_result(state.jω,state.jω*state.mu[1]/(2T(π))*bracket)
end

:Uribe2008
