function routes(identifier::Val{:Xue2021})
    return (
        self=FormulaMethod(identifier,earth_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,earth_impedance,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:Xue2021})=(air=_full,earth=_full,permeability=_material)
propagation(::Val{:Xue2021})=Val(:explicit)
media(::Formula{:Xue2021})=Val(:stratified)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinite parallel round conductors, heights ``h_i,h_j``, horizontal spacing ``y_{ij}``; self uses radius in ``D_1``. |
| Calculated quantities | Self and mutual overhead earth-return impedance above a four-layer earth |
| Earth structure | Four horizontal layers below air; fourth is semi-infinite. |
| Model and approximation | No fitted approximation is applied to (1); ``F`` is the exact four-layer interface factor within the source's line model. Its explicit algebraic dependencies ``A_1,A_2,T_1,\\ldots,T_6`` are printed in (A.1)–(A.15). |
| Main source | H. Xue, J. Mahseredjian, A. Ametani, J. Morales, and I. Kocar (2021) |
| Citation key(s) | `:Xue2021` |
| Evidence status | Accepted-publication page image verified |

**Expression.** Equations (1)–(3) and Appendix A are evaluated by the
existing permeability-weighted layer response, with
``a_n^2=\\lambda^2+\\gamma_n^2-\\gamma_0^2``.
The source's ``F`` is twice the shared magnetic surface kernel.
This is an exact algebraic elimination of the four-layer boundary
conditions, not the equivalent homogeneous earth approximation.

The final publication numbers Appendix A as (10)–(24).
Its depths are cumulative positive interface depths; engine inputs are
individual layer thicknesses. The two representations are compared
independently, including unequal layer permeabilities.

**Scope.** Exactly four earth layers, three finite thicknesses, lossless
air and overhead conductors are required. A supplied longitudinal
wavenumber must match the source air reference. Self uses the physical
radius only in the direct logarithm and zero lateral image separation.

**Reference.** [Xue2021](@cite), equations (1)–(3), (10)–(24).
"""
description(::Formula{:Xue2021}) =
    "Xue et al. exact four-layer overhead impedance (2021)"

function propagation_constant(::Val{:Xue2021},s,mu,epsilon)
    squared=oftype(s,-s^2*mu*epsilon)
    return (Γ=sqrt(squared),squared)
end

function earth_impedance(::Val{:Xue2021},::Val{:validate},rho,epsilon,mu,s,Γ,thickness)
    length(rho)==length(epsilon)==length(mu)==length(thickness)==5 ||
        throw(DimensionMismatch(":Xue2021 exact formula requires air and four earth layers"))
    isinf(rho[1]) && rho[1]>0 ||
        throw(DomainError(rho[1],":Xue2021 exact formula requires lossless air"))
    all(x->x>0 && (isfinite(x)||isinf(x)),rho) &&
        all(x->isfinite(x)&&x>0,epsilon) && all(x->isfinite(x)&&x>0,mu) ||
        throw(DomainError((rho,epsilon,mu),":Xue2021 exact formula requires positive material constants"))
    all(x->isfinite(x)&&x>=0,thickness[2:4]) &&
        thickness[1]==Inf && thickness[5]==Inf ||
        throw(DomainError(thickness,":Xue2021 requires three finite nonnegative layer thicknesses"))
    isfinite(s) && iszero(real(s)) && !iszero(s) ||
        throw(DomainError(s,":Xue2021 requires nonzero real frequency"))
    reference=propagation_constant(Val(:Xue2021),s,mu[1],epsilon[1])
    Γ===nothing || (isfinite(Γ)&&isapprox(Γ^2,reference.squared)) ||
        throw(ArgumentError(":Xue2021 exact formula fixes the air-reference longitudinal wavenumber"))
    return reference
end

function (formula::Formula{:Xue2021})(rho,epsilon,mu,s,Γ,segments=nothing)
    throw(DimensionMismatch(":Xue2021 exact formula requires the four-layer thickness vector"))
end

function (formula::Formula{:Xue2021})(rho,epsilon,mu,s,Γ,segments,thickness)
    reference=earth_impedance(Val(:Xue2021),Val(:validate),rho,epsilon,mu,s,Γ,thickness)
    functor=_stratified_functor(Val(:Xue2021),formula,rho,epsilon,mu,s,reference.Γ,segments,thickness)
    state=merge(functor.state,(gamma_squared=-functor.state.gamma_medium_squared[1],))
    return Functor{:Xue2021,typeof(formula.routes),typeof(state)}(formula.routes,state)
end

function earth_impedance(::Val{:Xue2021},::Val{:mutual},functor,pair)
    _require(pair,Val(:overhead))
    minimum(pair.heights)>0 || throw(DomainError(pair.heights,":Xue2021 requires positive heights"))
    return _layered_overhead_coefficient(functor.state,pair)
end

:Xue2021
