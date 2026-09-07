function routes(identifier::Val{:Lima2012})
    return (
        self=FormulaMethod(identifier,earth_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,earth_impedance,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:Lima2012}) = (air=_lossless,earth=_full,permeability=vacuum_permeability)
propagation(::Val{:Lima2012}) = Val(:zero)

function Formula(identifier::Val{:Lima2012};displacement_current::Bool=true,kwargs...)
    defaults=routes(identifier); overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) ||
        throw(ArgumentError("unknown Lima2012 routes"))
    selected=merge(defaults,overrides)
    values=merge(assumptions(identifier),(earth=displacement_current ? _full : _conductive,))
    return Formula{:Lima2012,typeof(selected),typeof(values)}(selected,values)
end

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filamentary external conductor positions at heights ``h_i,h_j``; horizontal separation ``d_{ij}``. |
| Calculated quantities | Exact closed-form evaluation of the quasi-TEM Carson self and mutual ground-return impedance using Struve and Bessel functions |
| Earth structure | Homogeneous half-space. |
| Model and approximation | Not an analytical approximation of the parent Carson integral (1). The source applies an integral transformation and the integral definition of the Struve function to obtain (7). The physical quasi-TEM/Carson model remains restricted. |
| Main source | A. C. S. de Lima and C. Portela (2012) |
| Citation key(s) | `:Lima2012` |
| Evidence status | Original publication page images checked. |

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | ``n`` parallel insulated single-core circular cables, outermost radius ``R\\ll h_i``; self sets ``x=R`` and ``\\ell=2h_i``. |
| Calculated quantities | Self and mutual closed-form asymptotic ground-return impedance for underground cables |
| Earth structure | Homogeneous earth. |
| Model and approximation | The source rewrites the parent residual ``J_m`` as (13)–(15), then expands the ``I_4`` integrand about ``\\xi\\to\\infty`` and retains the leading terms to obtain (16). The stated expansion is asymptotic and divergent; no convergence order or remainder bound is supplied. |
| Main source | A. C. S. de Lima and C. Portela (2012) |
| Citation key(s) | `:Lima2012` |
| Evidence status | Original publication page images checked |

**Expression.** With u₁,₂=η(H∓jx), the overhead coefficient is

```math
Z_{oh}=\\frac{j\\omega\\mu_0}{2\\pi}
\\left[\\ln(D/d)+F(u_1)+F(u_2)\\right],\\qquad
F(u)=\\frac{\\pi}{2u}[\\mathbf H_1(u)-Y_1(u)]-u^{-2}.
```

The distinct underground approximation is

```math
Z_{ug}=\\frac{j\\omega\\mu_0}{2\\pi}
\\left[K_0(\\eta d)+\\frac{H^2-x^2}{D^2}K_2(\\eta D)
-\\frac{2(H^2-x^2)}{\\eta^2D^4}(1+\\eta H)e^{-\\eta H}\\right].
```


**Selection.** Overhead pairs use the exact Struve evaluation of the
Carson-type integral; buried pairs use the distinct asymptotic equation (16).
The latter is restricted to horizontal separation smaller than half the
sum of burial depths.

For unequal overhead heights, both complex arguments use the height sum
from source integral (1). This preserves that parent and reciprocity.
The printed unequal-height arguments in (7) are not substituted for the
parent integral. Equal heights reproduce the printed closed form directly.
The engine adds the separate ideal-ground logarithm to return a complete
external coefficient.

**Numerical evaluation.** The overhead parent integral is evaluated directly
when small arguments, large arguments, or unsupported special-function
arithmetic would make the Struve difference inaccurate. The buried expression
shares the algebraically regularized Sunde image terms with DeConti2024,
but does not include the latter's Padé residual. These are distinct
approximations despite their common terms.

The former Theodoulidis2015 selection remains a conduction-only constructor
alias for the same overhead closed-form identity. Lima2012 retains the
source's full soil admittivity by default.

**Reference.** [Lima2012](@cite), equations (1), (7)–(8), and (16).
"""
description(::Formula{:Lima2012}) =
    "Lima–Portela overhead closed form and underground asymptote (2012)"

function propagation_constant(::Val{:Lima2012},jω,permeability,permittivity)
    return (Γ=zero(jω),squared=zero(jω))
end

function (formula::Formula{:Lima2012})(rho,epsilon,mu,jω,Γ,segments=nothing)
    return _homogeneous_functor(Val(:Lima2012),formula,rho,epsilon,mu,jω,Γ,segments)
end

function struve_h1(::Val{:Lima2012}, z::Complex{T}) where {T <: Real}
    term = 2z^2 / (3 * (one(T) * π))
    value = term
    tolerance = max(eps(T), T(1e-15))
    for k in 0:9999
        term *= -(z / 2)^2 / ((k + T(3) / 2) * (k + T(5) / 2))
        updated = value + term
        abs(term) <= tolerance * max(one(T), abs(updated)) && return updated
        value = updated
    end
    throw(ErrorException("complex Struve H1 series did not converge"))
end

@inline function closed_form_term(identifier::Val{:Lima2012}, u)
    return π / (2u) * (struve_h1(identifier, u) - bessely(1, u)) - inv(u^2)
end


function earth_impedance(tag::Val{:Lima2012},::Val{:self},functor,pair)
    return earth_impedance(tag,Val(:mutual),functor,pair)
end

function earth_impedance(tag::Val{:Lima2012},::Val{:mutual},functor,pair)
    state=functor.state; geometry=_geometry(pair)
    T=typeof(geometry.H); πT=one(T)*π
    if _placement(pair)===Val(:underground)
        geometry.y_ij<geometry.H/2 || throw(DomainError(geometry.y_ij,
            "Lima2012 underground asymptote requires separation < half the depth sum"))
        value=_sunde_image_terms(state.gamma[2],geometry.d_ij,geometry.H,geometry.y_ij)
        return _complex_result(state.jω,state.jω*state.mu[1]/(2πT)*value)
    end
    _require(pair,Val(:overhead))
    self=pair.row==pair.column
    lateral=self ? zero(T) : geometry.y_ij
    γ=state.gamma[2]
    u1=γ*(geometry.H-complex(zero(T),lateral))
    u2=γ*(geometry.H+complex(zero(T),lateral))
    correction=if T===BigFloat || min(abs(u1),abs(u2))<T(0.1) ||
        max(abs(u1),abs(u2))>T(15)
        integral=_carson_integral(state,geometry.H,lateral)
        2integral
    else
        closed_form_term(tag,u1)+closed_form_term(tag,u2)
    end
    ideal=self ? log(geometry.H/geometry.y_ij) : log(geometry.D_ij/geometry.d_ij)
    return _complex_result(state.jω,state.jω*state.mu[1]/(2πT)*(ideal+correction))
end

:Lima2012
