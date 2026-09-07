function routes(identifier::Val{:Kikuchi1957})
    return (
        self=FormulaMethod(identifier,earth_potential_coefficient,Val(:self)),
        mutual=FormulaMethod(identifier,earth_potential_coefficient,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant))
end
assumptions(::Val{:Kikuchi1957})=(air=_full,earth=_full,permeability=vacuum_permeability)
propagation(::Val{:Kikuchi1957})=Val(:explicit)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite thin circular overhead conductor of radius ``a`` and height ``h``; no insulation. The internal current distribution is approximated by that of an isolated conductor. |
| Calculated quantities | Scalar potentials ``V_1(x,y)`` in air and ``V_2(x,y)`` in earth due to one overhead wire; source-defined self terminal voltage ``V=V_1(0,h-a)``; charge/current relation retained separately |
| Earth structure | One flat homogeneous half-space ``y<0`` below air ``y>0``; surface ``y=0``. |
| Model and approximation | Hankel-plus-integral representation under the thin-wire, homogeneous-medium, and principal-mode assumptions. The small-Hankel and ``\\Gamma\\simeq jk_1`` reductions are auxiliary evaluators. |
| Main source | H. Kikuchi (1957); the 1955 Japanese paper is an earlier-language witness, and priority remains unresolved. |
| Citation key(s) | `:Kikuchi1957` |
| Evidence status | All 13 original pages checked; the auxiliary reductions are not independent physical formulations. |

**Expression.** With source charge continuity
``Q=\\Gamma_{source}I/(j\\omega)``, the returned coefficient is
``P=V/Q``. The engine wavenumber is ``k_x=j\\Gamma_{source}``,
so ``\\chi_n^2=\\gamma_n^2+k_x^2`` and the source's outgoing
Hankel argument is ``\\lambda_n=j\\chi_n``. Thus

```math
P=\\frac{1}{2\\pi\\varepsilon_1}
\\left[K_0(\\chi_1a)-K_0(\\chi_1[2h-a])
+2\\int_0^\\infty
\\frac{\\gamma_1^2e^{-(2h-a)q_1}}
{\\gamma_2^2q_1+\\gamma_1^2q_2}\\,du\\right],
\\qquad q_n=\\sqrt{u^2+\\chi_n^2}.
```

The equality ``K_0(\\chi)=j\\pi H_0^{(1)}(j\\chi)/2``
uses the same first-kind Hankel solution, not a replacement kind.
The default air-reference prescription implements the separate
source reduction (5.9); an explicit engine wavenumber retains the
general scalar potential without solving its modal dispersion relation.

The self boundary is the lower wire surface, so the image distance
and integral height are ``2h-a``, not ``2h``. Only the single-wire
terminal coefficient is registered. The separate observation leaf
also evaluates (5.2)/(5.4) below the interface; it is not a mixed-cable
matrix formula. No entrywise reciprocal is inserted into admittance.

**Reference.** [Kikuchi1957](@cite), (3.5),(5.1)–(5.4),(5.9), section 9;
[DLMF 10.27.8](https://dlmf.nist.gov/10.27.E8) for the Hankel connection.
"""
description(::Formula{:Kikuchi1957})="Kikuchi overhead-wire terminal scalar potential (1957)"
propagation_constant(::Val{:Kikuchi1957},s,mu,epsilon)=
    (Γ=sqrt(-s*s*mu*epsilon),squared=-s*s*mu*epsilon)

function (formula::Formula{:Kikuchi1957})(rho,epsilon,mu,s,Γ,segments=nothing)
    _check(rho,epsilon,mu)
    isinf(rho[1]) || throw(DomainError(rho[1],"Kikuchi scalar source requires lossless air"))
    all(x->isapprox(x,vacuum_permeability(x)),mu) ||
        throw(DomainError(mu,"Kikuchi scalar potential requires nonmagnetic media"))
    return _scalar_potential_functor(Val(:Kikuchi1957),formula,rho,epsilon,mu,s,Γ,segments)
end

function earth_potential_coefficient(::Val{:Kikuchi1957},::Val{:self},functor,pair)
    _require(pair,Val(:overhead))
    pair.row==pair.column && pair.heights[1]==pair.heights[2] ||
        throw(ArgumentError("Kikuchi terminal coefficient requires self geometry"))
    h=pair.heights[1]; a=pair.separation
    0<a<h || throw(DomainError((a,h),"the complete wire section must lie above the interface"))
    return _scalar_potential_coefficient(functor.state,1,a,2h-a,2h-a,zero(h))
end

function earth_potential_coefficient(::Val{:Kikuchi1957},::Val{:mutual},functor,pair)
    throw(ArgumentError("Kikuchi's observation potential is not a registered multiwire mutual coefficient"))
end

function earth_potential_coefficient(::Val{:Kikuchi1957},::Val{:observation},
        functor,h,x,y)
    all(isfinite,(h,x,y)) && h>0 || throw(DomainError((h,x,y),"invalid observation coordinates"))
    state=functor.state
    if y>=0
        return _scalar_potential_coefficient(state,1,hypot(x,y-h),
            hypot(x,y+h),y+h,abs(x))
    end
    bulk=state.gamma_medium_squared
    integral=_scalar_potential_integral(state,h-y,1) do u
        q1=spectral_root(u^2+bulk[1]+state.gamma_squared,state.jω)
        q2=spectral_root(u^2+bulk[2]+state.gamma_squared,state.jω)
        bulk[1]*exp(y*q2-h*q1)*cos(x*u)/(bulk[2]*q1+bulk[1]*q2)
    end
    return _complex_result(state.jω,integral/((one(h)*π)*state.epsilon[1]))
end

:Kikuchi1957
