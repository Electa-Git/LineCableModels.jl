function routes(identifier::Val{:Wedepohl1966})
    return (
        self=FormulaMethod(identifier,earth_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,earth_impedance,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:Wedepohl1966})=(air=_conductive,earth=_conductive,permeability=_material)
propagation(::Val{:Wedepohl1966})=Val(:explicit)
media(::Formula{:Wedepohl1966})=Val(:stratified)
function Formula(::Val{:Wedepohl1966};displacement_current::Bool=false,kwargs...)
    identifier=Val(:Wedepohl1966); defaults=routes(identifier)
    if displacement_current
        defaults=merge(defaults,(Γ=FormulaMethod(
            identifier,propagation_constant,Val(:air_reference)),))
    end
    overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) || throw(ArgumentError(
        "unknown routes for earth-impedance formula :Wedepohl1966"))
    values=displacement_current ? (air=_full,earth=_full,permeability=_material) :
        assumptions(identifier)
    return Formula(identifier,merge(defaults,overrides),values)
end

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Mutual pair at heights ``h_r,y`` and horizontal coordinates ``s_r,x``. Filamentary external conductors; no insulation term. |
| Calculated quantities | Mutual p.u.l. series-impedance contribution of a multilayer earth; conduction-only and displacement-current variants |
| Earth structure | Finite upper layer of thickness ``d`` above a lower half-space. |
| Model and approximation | Equation (8) neglects displacement current. Equation (9) restores it after prescribing free-space longitudinal propagation; it is not a solved full-wave modal result. |
| Main source | L. M. Wedepohl and R. G. Wasley (1966) |
| Citation key(s) | `:Wedepohl1966` |
| Evidence status | Original publication page images checked |

**Expression.** The source's hyperbolic expression can be written without
growing functions as
``A=q_2[(\\mu_2/\\mu_3)q_3+q_2\\tanh(dq_2)]/
[q_2+(\\mu_2/\\mu_3)q_3\\tanh(dq_2)]``.
The ground correction is

```math
Z_g=\\frac{j\\omega\\mu_1}{\\pi}\\int_0^\\infty
\\frac{e^{-\\lambda H}\\cos(\\lambda x)}
{\\lambda+(\\mu_1/\\mu_2)A}\\,d\\lambda.
```

The default uses equation (8), with conduction-only roots. Select
`displacement_current=true` for equation (9), which uses bulk soil roots
minus the air constant squared. Independent layer permeabilities are
retained in both cases. The complete coefficient adds the ideal-ground
logarithm.

The same weighted boundary relation is shared with Sunde's nonmagnetic
case and Moghram's three-layer extension. The self coefficient uses the
thin-wire limit: zero lateral spacing in the ground correction and the
wire radius in the ideal-ground logarithm. This assembly convention is
not a separately printed finite-radius formula in Wedepohl's mutual record.

**Reference.** [Wedepohl1966](@cite), equations (8)–(9).
"""
description(::Formula{:Wedepohl1966})="Wedepohl–Wasley permeable two-layer overhead impedance (1966)"
propagation_constant(::Val{:Wedepohl1966},jω,permeability,permittivity)=(Γ=zero(jω),squared=zero(jω))
function propagation_constant(::Val{:Wedepohl1966},::Val{:air_reference},jω,permeability,permittivity)
    squared=oftype(jω,-jω^2*permeability*permittivity)
    return (Γ=sqrt(squared),squared)
end
function (formula::Formula{:Wedepohl1966})(rho,epsilon,mu,jω,Γ,segments,thickness)
    length(rho)==3 || throw(DimensionMismatch("Wedepohl1966 requires exactly two earth layers"))
    reference=formula.routes.Γ(jω,mu[1],epsilon[1])
    selected=isnothing(Γ) ? reference.Γ : Γ
    selected≈reference.Γ || throw(ArgumentError(
        "Wedepohl1966 uses the propagation reference prescribed by its selected variant"))
    return _stratified_functor(
        Val(:Wedepohl1966),formula,rho,epsilon,mu,jω,selected,segments,thickness)
end
function earth_impedance(::Val{:Wedepohl1966},::Val{:mutual},functor,pair)
    _require(pair,Val(:overhead))
    return _layered_overhead_coefficient(functor.state,pair)
end

:Wedepohl1966
