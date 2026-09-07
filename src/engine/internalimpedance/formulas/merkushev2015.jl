function routes(identifier::Val{:Merkushev2015})
    return (
        inner=FormulaMethod(identifier,internal_impedance,Val(:inner)),
        outer=FormulaMethod(identifier,internal_impedance,Val(:outer)),
        mutual=FormulaMethod(identifier,internal_impedance,Val(:mutual)))
end
assumptions(::Val{:Merkushev2015})=(;)
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | One central steel wire and six aluminum strands, homogenized as an anisotropic layer; strand radius ``R`` and helix pitch ``h``. |
| Calculated quantities | Low-frequency p.u.l. internal impedance of a single-layer aluminum conductor helically stranded over a magnetic steel core |
| Earth structure | None. |
| Model and approximation | Six discrete strands are homogenized into a continuous anisotropic surface layer; aluminum skin effect is ignored. The core solution retains cylindrical skin effect through Bessel functions. |
| Main source | A. G. Merkushev and I. A. Elagin (2015) |
| Citation key(s) | `:Merkushev2015` |
| Evidence status | Original publication page images checked |

**Expression.** The six-strand anisotropic layer supplies
``\\Sigma_z=6\\pi R^2\\sigma_{Al}Q``. Set
``w^2=j\\omega\\mu_0\\mu_{St}\\sigma_{St}R^2``,
``q=I_1(w)/(wI_0(w))``, and
``C=2\\Sigma_c/\\Sigma_z``. The source's distinct fractions in (2)
give the equivalent numerical expression

```math
Z_{int}=\\frac{C+\\theta^2w^2q}
{2\\Sigma_c[1+(C+\\theta^2w^2q)q]},
\\qquad\\Sigma_c=\\pi R^2\\sigma_{St},\\quad\\theta=2\\pi R/h.
```

The zero-frequency limit is ``1/(\\Sigma_c+\\Sigma_z)``.
The geometry is one central wire and six equal-radius helical strands,
not a solid equivalent annulus. The coaxial adapter retains both
materials and the declared pitch. It does not solve a three-dimensional
helical field or nonlinear magnetic saturation.

The necessary low-frequency condition
``|\\omega|\\mu_0\\sigma_{Al}R^2<1`` is enforced; the source's
stronger asymptotic condition is much less than one. This check is
not a claimed error bound. Only the scalar outer term is available.

**Reference.** [Merkushev2015](@cite), (1)–(2), p. 402.
"""
description(::Formula{:Merkushev2015})="Merkushev linear single-layer ACSR impedance (2015)"

function (formula::Formula{:Merkushev2015})(R::T,h::T,rho_core::T,
        rho_strand::T,mu_core::T,s::Complex{T}) where {T <: AbstractFloat}
    all(isfinite,(R,rho_core,rho_strand,mu_core,s)) &&
        R>0 && h>0 && !isnan(h) && rho_core>0 && rho_strand>0 && mu_core>0 &&
        iszero(real(s)) || throw(DomainError((R,h,rho_core,rho_strand,mu_core,s),
            "Merkushev requires positive physical inputs and a finite real frequency"))
    mu0=vacuum_permeability(R)
    parameter=abs(s)*mu0*R^2/rho_strand
    parameter<1 || throw(DomainError(parameter,
        "Merkushev's aluminium-skin-effect expansion requires |omega|*mu0*sigmaAl*R^2 much smaller than one"))
    theta=isinf(h) ? zero(T) : 2*(one(T)*π)*R/h
    Q=if iszero(theta)
        one(T)
    else
        integral=quadgk(one(T),T(3);rtol=max(64eps(T),T(1e-12))) do r
            r*acos(clamp((r^2+3)/(4r),-one(T),one(T)))/(1+theta^2*r^2)
        end
        2/(one(T)*π)*integral[1]
    end
    core=(one(T)*π)*R^2/rho_core
    layer=6*(one(T)*π)*R^2/rho_strand*Q
    w2=s*mu0*mu_core*R^2/rho_core
    w=sqrt(w2)
    q=iszero(w) ? complex(one(T)/2,zero(T)) :
        special_besselix(1,w)/(w*special_besselix(0,w))
    C=2core/layer
    U=C+theta^2*w2*q
    outer=Complex{T}(U/(2core*(1+U*q)))
    state=(;R,h,rho_core,rho_strand,mu_core,s,Q,outer)
    return Functor{:Merkushev2015,typeof(formula.routes),typeof(state)}(formula.routes,state)
end

function (formula::Formula{:Merkushev2015})(r_in,r_ex,rho,mu,s)
    throw(ArgumentError("Merkushev requires the retained seven-wire geometry, pitch, and two physical materials"))
end
internal_impedance(::Val{:Merkushev2015},::Val{:outer},state)=state.outer
internal_impedance(::Val{:Merkushev2015},::Union{Val{:inner},Val{:mutual}},state)=zero(state.s)
(functor::Functor{:Merkushev2015})(::Val{:outer})=functor.routes.outer(functor.state)
(functor::Functor{:Merkushev2015})(::Val{:inner})=functor.routes.inner(functor.state)
(functor::Functor{:Merkushev2015})(::Val{:mutual})=functor.routes.mutual(functor.state)

:Merkushev2015
