function routes(identifier::Val{:Lee2014})
    return (
        self=FormulaMethod(identifier,earth_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,earth_impedance,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:Lee2014})=(air=_full,earth=_full,permeability=vacuum_permeability)
propagation(::Val{:Lee2014})=Val(:explicit)
media(::Formula{:Lee2014})=Val(:stratified)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Mutual pair; self follows standard substitution. Filamentary external model; no insulation. |
| Calculated quantities | Arbitrary-horizontal-layer overhead earth-return mutual impedance and accelerated numerical evaluator |
| Earth structure | Arbitrary ``N`` horizontal layers. |
| Model and approximation | Spline interpolation and asymptotic-tail truncation are numerical approximations to the parent integral. |
| Main source | Jae-bok Lee, Jun Zou, Benliang Li, and Munno Ju (2014) |
| Citation key(s) | `:Lee2014` |
| Evidence status | Original publication page images checked |

**Expression.** The source's nonmagnetic multilayer parent uses
uₘ²=λ²+γₘ²−γ₀², with the air-reference subtraction retained.
Equations (4)–(9) determine the reflection coefficients, and equation (2)
integrates their surface response. The ideal-ground logarithm from (1) is
added to produce the complete external coefficient.

**Shared evaluation.** The boundary-matching equations are evaluated by
the existing multilayer kernel with the air longitudinal wavenumber supplied
explicitly. This differs from the default zero-longitudinal prescription
of Tsiamitros2008. No second boundary solver is introduced.

Adaptive quadrature evaluates the unapproximated parent integral.
The paper's spline moments and asymptotic tail are alternative numerical
evaluations, not extra physical terms; their finite truncation is not
reproduced by this selection.

**Reference.** [Lee2014](@cite), equations (1)–(9), printed pp. 1381–1382.
"""
description(::Formula{:Lee2014}) =
    "Lee et al. air-referenced nonmagnetic multilayer parent impedance (2014)"

function propagation_constant(::Val{:Lee2014},jω,permeability,permittivity)
    squared=oftype(jω,-jω^2*permeability*permittivity)
    return (Γ=sqrt(squared),squared)
end

function (formula::Formula{:Lee2014})(rho,epsilon,mu,jω,Γ,segments=nothing)
    length(rho)==2 || throw(DimensionMismatch("Lee2014 requires thicknesses for multiple earth layers"))
    return formula(rho,epsilon,mu,jω,Γ,segments,fill(oftype(first(rho),Inf),2))
end

function (formula::Formula{:Lee2014})(rho,epsilon,mu,jω,Γ,segments,thickness)
    length(rho)>=2 || throw(DimensionMismatch("Lee2014 requires air and at least one earth layer"))
    reference = isnothing(Γ) ? formula.routes.Γ(
        jω, _permeability(mu, 1, formula.assumptions.permeability), epsilon[1]
    ).Γ : Γ
    return _stratified_functor(
        Val(:Lee2014),formula,rho,epsilon,mu,jω,reference,segments,thickness
    )
end

function earth_impedance(::Val{:Lee2014},::Val{:mutual},functor,pair)
    _require(pair,Val(:overhead))
    return earth_impedance(Val(:Tsiamitros2008),Val(:overhead),functor,pair)
end

:Lee2014
