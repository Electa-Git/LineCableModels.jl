function routes(identifier::Val{:Moghram1998})
    return (
        self=FormulaMethod(identifier,earth_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,earth_impedance,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:Moghram1998})=(air=_conductive,earth=_conductive,permeability=_material)
propagation(::Val{:Moghram1998})=Val(:zero)
media(::Formula{:Moghram1998})=Val(:stratified)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Self and mutual overhead conductors. Filament source with radius used for self evaluation. |
| Calculated quantities | Overhead-line self and mutual impedances over three-layer earth, with two-layer and homogeneous limits |
| Earth structure | Three conducting layers below air; reductions to two/one. |
| Model and approximation | Displacement and longitudinal variation are neglected; within that model the DCFT layer solution is not a fitted approximation. |
| Main source | I. S. Moghram (1998), extending Wedepohl–Wasley |
| Citation key(s) | `:Moghram1998` |
| Evidence status | Original publication equations and text checked |

**Expression.** Equations (14)–(15) reduce to the weighted surface field
ratio obtained by successively matching the two finite earth layers.
With ``q_k^2=\\lambda^2+j\\omega\\mu_k\\sigma_k``, start with
``b_4=q_4/\\mu_4`` and apply

```math
b_k=\\frac{q_k}{\\mu_k}
\\frac{b_{k+1}+(q_k/\\mu_k)\\tanh(q_kt_k)}
{(q_k/\\mu_k)+b_{k+1}\\tanh(q_kt_k)},\\qquad k=3,2.
```

The kernel in (16)–(17) is ``1/(\\lambda+\\mu_1b_2)``.
Equation (16) uses zero lateral spacing in its ground correction, with
the radius only in ``\\ln(2h/r)``. Equation (17) retains the mutual
cosine and direct/image distances.

The implementation shares the boundary calculation with the two-layer
formula. It accepts the source's one-, two-, and three-layer cases;
no arbitrary-layer claim is attributed to this publication. Displacement
current is omitted, and layer permeabilities remain independent.

**Reference.** [Moghram1998](@cite), equations (14)–(19), p. 447.
"""
description(::Formula{:Moghram1998})="Moghram permeable three-layer overhead impedance (1998)"
propagation_constant(::Val{:Moghram1998},jω,permeability,permittivity)=(Γ=zero(jω),squared=zero(jω))
function (formula::Formula{:Moghram1998})(rho,epsilon,mu,jω,Γ,segments,thickness)
    2<=length(rho)<=4 || throw(DimensionMismatch(
        "Moghram1998 supports one to three earth layers"))
    return _stratified_functor(
        Val(:Moghram1998),formula,rho,epsilon,mu,jω,Γ,segments,thickness)
end
function (formula::Formula{:Moghram1998})(rho,epsilon,mu,jω,Γ,segments=nothing)
    length(rho)==2 || throw(DimensionMismatch(
        "Moghram1998 requires thicknesses for multiple earth layers"))
    return formula(rho,epsilon,mu,jω,Γ,segments,fill(oftype(first(rho),Inf),2))
end
function earth_impedance(::Val{:Moghram1998},::Val{:mutual},functor,pair)
    _require(pair,Val(:overhead))
    return _layered_overhead_coefficient(functor.state,pair)
end

:Moghram1998
