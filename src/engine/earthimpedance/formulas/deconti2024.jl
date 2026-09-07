function routes(identifier::Val{:DeConti2024})
    return (
        self=FormulaMethod(identifier,earth_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,earth_impedance,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:DeConti2024}) = (air=_lossless,earth=_full,permeability=_material)
propagation(::Val{:DeConti2024}) = Val(:zero)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Circular insulated cables represented externally by their cable coordinates; ``h_m,h_n`` are burial depths and ``r`` is horizontal separation (or the source's self-distance substitution). |
| Calculated quantities | Self and mutual per-unit-length ground-return impedance of underground cables; closed-form 1/1 Padé approximation of Sunde's residual integral |
| Earth structure | Homogeneous earth; the displayed Sunde parent has no upper-medium spectral term. |
| Model and approximation | Starting from Sunde's exact-in-model rearrangement (3)–(4), the source expands about ``t=1`` and replaces only ``e^{-t\\gamma D}`` by the 1/1 Padé form ``[(2-\\gamma D(t-1))/(2+\\gamma D(t-1))]e^{-\\gamma D}`` in (6). Integrating that rational replacement gives (7)–(10), which is inserted into (5). No additional term is discarded in the printed construction. |
| Main source | A. De Conti and A. C. S. de Lima (2024) |
| Citation key(s) | `:DeConti2024` |
| Evidence status | Original publication page images checked |

**Expression.** Equations (5) and (7)–(10) supply the Padé approximation.
The earth diffusion constant retains physical permeability and full
admittivity; the impedance prefactor remains the source's vacuum permeability.

**Numerical evaluation.** At small arguments, the singular leading terms of
the Bessel and exponential contributions are cancelled algebraically.
The convergent remainder series is summed to working precision.
For the Padé residual, a finite angular integral evaluates the same rational
approximant when its closed expression would subtract large nearly equal
terms. This changes the evaluation, not the approximated exponential.
At larger arguments, the printed arctangent expression is used directly.
No air-wave term or external-admittance formula is inferred.

**Reference.** [DeConti2024](@cite), equations (5)–(10).
"""
description(::Formula{:DeConti2024}) =
    "De Conti–de Lima Padé underground impedance (2024)"

function propagation_constant(::Val{:DeConti2024},jω,permeability,permittivity)
    return (Γ=zero(jω),squared=zero(jω))
end

function (formula::Formula{:DeConti2024})(rho,epsilon,mu,jω,Γ,segments=nothing)
    return _homogeneous_functor(
        Val(:DeConti2024),formula,rho,epsilon,mu,jω,Γ,segments
    )
end


function closed_form_term(tag::Val{:DeConti2024},::Val{:pade_residual},z,H,r,D)
    T=typeof(H)
    iszero(r) && return zero(z)
    if abs(z)<one(T)
        angle=atan(r,H)
        integrand(θ)=-cos(2θ)*(2-z*(cos(θ)-1))/(2+z*(cos(θ)-1))
        return exp(-z)*quadgk(integrand,zero(T),angle;rtol=eps(T)^(T(3)/4))[1]
    end
    u=r/(H+D); q=sqrt(1-z)
    i1=(H/D-8/z)*r/D
    i2=16*(2-z)/z^2*atan(u)
    i3=-4*(8-8z+z^2)/z^2*atan(u*q)/q
    return exp(-z)*(i1+i2+i3)
end

function closed_form_term(tag::Val{:DeConti2024},γ,d,H,r)
    D=hypot(H,r); z=γ*D
    residual=closed_form_term(tag,Val(:pade_residual),z,H,r,D)
    return _sunde_image_terms(γ,d,H,r)-2r*H/D^2*residual
end

function earth_impedance(tag::Val{:DeConti2024},::Val{:mutual},functor,pair)
    _require(pair,Val(:underground))
    state=functor.state; geometry=_geometry(pair)
    value=closed_form_term(tag,state.gamma[2],geometry.d_ij,geometry.H,geometry.y_ij)
    μ0=vacuum_permeability(one(geometry.H))
    return _complex_result(state.jω,state.jω*μ0/(2*(one(geometry.H)*π))*value)
end

:DeConti2024
