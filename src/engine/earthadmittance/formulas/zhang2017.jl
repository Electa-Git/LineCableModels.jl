function routes(identifier::Val{:Zhang2017})
    return (
        self=FormulaMethod(identifier,earth_potential_coefficient,Val(:self)),
        mutual=FormulaMethod(identifier,earth_potential_coefficient,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:Zhang2017}) = (
    air=_vacuum,earth=_full,permeability=vacuum_permeability,
    evaluation=:integral,source_radius=nothing
)
propagation(::Val{:Zhang2017}) = Val(:zero)

function Formula(::Val{:Zhang2017};evaluation::Symbol=:integral,source_radius=nothing,kwargs...)
    evaluation in (:integral,:asymptotic_tail) ||
        throw(ArgumentError(":Zhang2017 evaluation must be :integral or :asymptotic_tail"))
    source_radius===nothing ||
        (source_radius isa Real && isfinite(source_radius) && source_radius>0) ||
        throw(DomainError(source_radius,":Zhang2017 source_radius must be positive and finite"))
    defaults=routes(Val(:Zhang2017)); overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) ||
        throw(ArgumentError("unknown routes for earth-admittance formula :Zhang2017"))
    selected=merge(defaults,overrides)
    values=merge(assumptions(Val(:Zhang2017)),(;evaluation,source_radius))
    return Formula{:Zhang2017,typeof(selected),typeof(values)}(selected,values)
end

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | One buried insulated conductor at depth ``h``; equation (11) samples at the metal-core radius ``r``, distinct from the outer insulation radius ``R``. |
| Calculated quantities | Ground-admittance integral and asymptotic-extraction numerical evaluator |
| Earth structure | Homogeneous earth below air. |
| Model and approximation | The evaluator splits the spectral integral at a source-defined threshold, integrates spline moments on the finite interval, and evaluates the extracted tail with exponential integrals. It retains earth displacement current. |
| Main source | Boyuan Zhang, Jun Zou, Xuelong Du, Jaebok Lee, and Mun-No Ju (2017) |
| Citation key(s) | Primary: `:Zhang2017`; parent integral: `:Papadopoulos2010b`; Vance reduction: `:Vance1978` |
| Evidence status | Original PDF equations (11)–(18) and Appendix A checked |

**Expression.** The right-hand side of printed equation (11) is evaluated
as a potential coefficient, with ground admittance obtained from
``Y_g=j\\omega/P_g``. With ``s=j\\omega``,
``a=s\\mu_0\\sigma_g``, ``q=\\sqrt{\\lambda^2+a}``,
``c=\\gamma_0^2/\\gamma_1^2``, and ``B=\\sqrt{4h^2+r^2}``,
the exact integral is evaluated as

```math
P_g=\\frac{s}{2\\pi(\\sigma_g+s\\varepsilon_g)}
\\left[
K_0(\\sqrt a\\,r)-K_0(\\sqrt a\\,B)
+2\\int_0^\\infty\\frac{e^{-2hq}\\cos(\\lambda r)}
{\\lambda+cq}\\,d\\lambda
\\right].
```

The Bessel difference integrates the nondecaying part of (12) exactly.
The source fixes different lossless axial references in air and earth:
``u_0=\\lambda``, ``u_1=q``. The interface weights retain full
bulk propagation constants. This is not the single axial reference used
by the default Papadopoulos2010b selection.

**Evaluation.** The default `evaluation=:integral` uses the unapproximated
kernel. `evaluation=:asymptotic_tail` integrates the exact kernel up to
``T=10|\\sqrt a|`` and evaluates (A.6)–(A.8) above it.
A moment of order ``n`` includes ``T^{1-n}E_n(Tz)``;
in particular, cubic inverse powers require ``T^{-2}``.
Adaptive quadrature evaluates the finite interval; the paper's separate
spline-moment implementation and timing are not reproduced.

**Geometry and scope.** Only buried self coefficients are defined.
The pair radius is the sampling radius by default. Standard coaxial
assembly supplies the outer earth-interface radius; `source_radius=r`
selects the printed metal-core sampling explicitly for a coated source.
The source's choice and the matrix interface choice are therefore not
silently equated. Internal insulation terms remain separate. Two
nonmagnetic half-spaces, vacuum air, positive earth conductivity, and
nonzero frequency are required. A nonzero common axial input is rejected.
The normal self route also rejects a negative real ground admittance.
Such values occur in some high-frequency cases of this prescribed kernel;
they are not made passive by changing the source coefficients.

**Reference.** [Zhang2017](@cite), equations (11)–(18) and Appendix A.
The source's logarithmic/Vance alternative uses the existing
Petrache2005 impedance and Vance1978 scalar conversion; it is not added
as another copy of those formulas.
"""
description(::Formula{:Zhang2017}) =
    "Zhang et al. buried-wire potential integral and asymptotic tail (2017)"

function propagation_constant(::Val{:Zhang2017},s,mu,epsilon)
    return (Γ=zero(s),squared=zero(s))
end

function (formula::Formula{:Zhang2017})(
        rho::AbstractVector{T},epsilon::AbstractVector{T},mu::AbstractVector{T},
        s::Complex{T},Γ,segments=nothing
) where {T <: Real}
    _check(rho,epsilon,mu)
    length(rho)==2 || throw(ArgumentError(":Zhang2017 requires one earth half-space"))
    isinf(rho[1]) && isfinite(rho[2]) && rho[2]>zero(T) ||
        throw(DomainError(rho,":Zhang2017 requires vacuum air and conducting earth"))
    all(x->isfinite(x)&&x>zero(T),epsilon) &&
        isapprox(epsilon[1],vacuum_permittivity(epsilon[1])) ||
        throw(DomainError(epsilon,":Zhang2017 requires vacuum air and positive permittivity"))
    all(x->isapprox(x,vacuum_permeability(x)),mu) ||
        throw(DomainError(mu,":Zhang2017 assumes nonmagnetic media"))
    isfinite(s) && !iszero(s) && iszero(real(s)) ||
        throw(DomainError(s,":Zhang2017 requires nonzero real frequency"))
    return _homogeneous_functor(Val(:Zhang2017),formula,rho,epsilon,mu,s,Γ,segments)
end

function earth_potential_coefficient(
        ::Val{:Zhang2017},::Val{:tail_moment},n::Int,z::Complex{T}
) where {T <: AbstractFloat}
    n in (1,3) || throw(ArgumentError("the extracted tail uses moments of order 1 and 3"))
    e=exp(-z)
    n==1 && return e*scaled_expint_negative(-z)
    if abs(z)<=2
        first=e*scaled_expint_negative(-z)
        second=e-z*first
        return (e-z*second)/2
    end
    # Laplace representation avoids cancellation in upward E_n recurrence.
    value=quadgk(v->v*v*exp(-v)/(v+z),zero(T),T(Inf);
        rtol=max(sqrt(eps(T)),T(1e-30)))[1]
    return e*value/2
end

function earth_potential_coefficient(::Val{:Zhang2017},::Val{:self},functor,pair)
    value=earth_potential_coefficient(Val(:Zhang2017),Val(:coefficient),functor,pair)
    state=functor.state
    sign(imag(state.jω))*imag(value)>=-state.tolerance*abs(value) ||
        throw(DomainError(state.jω,
            ":Zhang2017 gives negative ground conductance here; select a different potential model"))
    return value
end

function earth_potential_coefficient(::Val{:Zhang2017},::Val{:coefficient},functor,pair)
    _require(pair,Val(:underground))
    pair.row==pair.column || throw(ArgumentError(":Zhang2017 defines only a scalar self coefficient"))
    state=functor.state; T=typeof(real(state.jω))
    h=abs(first(pair.heights))
    selected=state.formula.assumptions.source_radius
    r=selected===nothing ? T(pair.separation) : T(selected)
    isfinite(h) && h>0 && isfinite(r) && r>0 ||
        throw(DomainError((h,r),":Zhang2017 requires positive depth and sampling radius"))
    a=state.jω*state.mu[2]*state.sigma[2]
    g=sqrt(a); bulk=state.gamma_medium_squared
    c=bulk[1]/bulk[2]; H=2h
    if state.formula.assumptions.evaluation===:integral
        integrand=t->begin
            lambda=t/H; q=sqrt(lambda*lambda+a)
            exp(-H*q)*cos(lambda*r)/(lambda+c*q)/H
        end
        transition=min(one(T),abs(c*g)*H)
        reflection=quadgk(integrand,zero(T),transition,one(T),T(Inf);
            rtol=state.tolerance)[1]
        direct=special_besselk(0,g*r)-special_besselk(0,g*hypot(H,r))
        total=direct+2reflection
    else
        threshold=10abs(g)
        integrand=t->begin
            lambda=threshold*t; q=sqrt(lambda*lambda+a)
            decay=exp(-H*q)
            (-expm1(-H*q)/q+2decay/(lambda+c*q))*cos(lambda*r)*threshold
        end
        transition=min(one(T),abs(c*g)/threshold)
        finite=quadgk(integrand,zero(T),transition,one(T);
            rtol=state.tolerance)[1]
        moment(n,z)=real(earth_potential_coefficient(
            Val(:Zhang2017),Val(:tail_moment),n,threshold*z))*threshold^(1-n)
        pure=complex(zero(T),r); damped=complex(H,r)
        R=(bulk[2]-bulk[1])/(bulk[2]+bulk[1])
        kx1_squared=-state.jω^2*state.mu[2]*state.epsilon[2]
        tail=moment(1,pure)-a*moment(3,pure)/2+a*moment(3,damped)/2+
            R*(moment(1,damped)+(bulk[2]-kx1_squared)*moment(3,damped)/2)
        total=finite+tail
    end
    return _complex_result(state.jω,state.jω*total/
        (2T(π)*(state.sigma[2]+state.jω*state.epsilon[2])))
end

function earth_potential_coefficient(::Val{:Zhang2017},::Val{:mutual},functor,pair)
    throw(ArgumentError(":Zhang2017 supplies no mutual potential coefficient"))
end

:Zhang2017
