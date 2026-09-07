function routes(identifier::Val{:Maaouni2001})
    return (
        self=FormulaMethod(identifier,earth_potential_coefficient,Val(:self)),
        mutual=FormulaMethod(identifier,earth_potential_coefficient,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:Maaouni2001}) = (
    air=_full,earth=_full,permeability=vacuum_permeability
)
propagation(::Val{:Maaouni2001}) = Val(:explicit)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Parallel thin bare overhead wires; scalar self/mutual coefficients; no coating. |
| Calculated quantities | Thin-wire potential coefficient and qTEM admittance near a lossy interface |
| Earth structure | Homogeneous earth below lossless air. |
| Model and approximation | qTEM limit of the air transverse constant, followed by the exponential approximation (23) used in (25)–(26); no independent full-wave modal solve. |
| Main source | A. Maaouni, A. Amri, and A. Zouhir (2001) |
| Citation key(s) | `:Maaouni2001` |
| Evidence status | Original publication page images checked |

**Expression.** The qTEM potential coefficient is
``P_{mn}=[\\ln(D_{mn}/d_{mn})+G_{mn}]/(2\\pi\\varepsilon_{air})``,
where ``G`` is source equation (26) with its helper (25).
The source uses ``e^{-j\\omega t}``; the evaluator uses its conjugate
under ``e^{j\\omega t}``.

Writing ``N=n^2``, ``A=\\sqrt{1-N}``,
``b=-j/\\sqrt{1+N}`` at positive frequency,
``z_\\pm=k_0(H\\pm jX)``, and
``Q(w)=e^{-w}E_1(-w)``, cancellation of the intermediate ``P(b,z)``
terms gives the numerically evaluated expression

```math
G=\\frac14\\sum_{z\\in\\{z_+,z_-\\}}
\\left[
2\\ln\\left(1+\\frac{2}{Az}\\right)
+Q\\left(b[z+2/A]\\right)+Q\\left(-b[z+2/A]\\right)
-\\frac{N-1}{N+1}\\{Q(bz)+Q(-bz)\\}
\\right].
```

This is algebraically equation (26), not replacement by the unapproximated
Wise1948 integral. Both ``\\Re(Az_\\pm)>0`` conditions are enforced.
The phase-domain potential matrix is assembled before inversion.
For self terms, ``d=r``, ``D=2h``, and ``X=0``.

**Assumptions.** Two nonmagnetic half-spaces, lossless air, conducting earth,
and overhead thin wires. The axial wavenumber is fixed to the real air
wavenumber by the qTEM reduction; a different explicit longitudinal input
is rejected. No full-wave mode, buried pair, or mixed pair is inferred.
Wider arithmetic retains cancellation in the closed expression without
changing the caller's global precision.

**Reference.** [Maaouni2001](@cite), equations (2), (7), (23), (25)–(26).
The unapproximated qTEM potential integral is [Wise1948](@cite).
"""
description(::Formula{:Maaouni2001}) =
    "Maaouni et al. overhead qTEM analytical potential approximation (2001)"

function propagation_constant(::Val{:Maaouni2001},s,mu,epsilon)
    squared=-s^2*mu*epsilon
    return (Γ=sqrt(squared),squared)
end

function (formula::Formula{:Maaouni2001})(
        rho::AbstractVector{T},epsilon::AbstractVector{T},mu::AbstractVector{T},
        s::Complex{T},Γ,segments=nothing
) where {T <: Real}
    _check(rho,epsilon,mu)
    length(rho)==2 || throw(ArgumentError(":Maaouni2001 requires two homogeneous half-spaces"))
    isinf(rho[1]) && isfinite(rho[2]) && rho[2]>zero(T) ||
        throw(DomainError(rho,":Maaouni2001 requires lossless air and conducting earth"))
    all(x->isfinite(x)&&x>zero(T),epsilon) ||
        throw(DomainError(epsilon,":Maaouni2001 requires positive finite permittivities"))
    all(x->isapprox(x,vacuum_permeability(x)),mu) ||
        throw(DomainError(mu,":Maaouni2001 assumes nonmagnetic media"))
    isfinite(s) && !iszero(s) && iszero(real(s)) ||
        throw(DomainError(s,":Maaouni2001 requires nonzero real frequency"))
    reference=formula.routes.Γ(s,mu[1],epsilon[1]).Γ
    Γ===nothing || isapprox(Γ,reference) ||
        throw(ArgumentError(":Maaouni2001 fixes the axial wavenumber to the air reference"))
    return _homogeneous_functor(Val(:Maaouni2001),formula,rho,epsilon,mu,s,reference,segments)
end

function earth_potential_coefficient(
        ::Val{:Maaouni2001},::Val{:kernel},
        N::Complex{T},k0::T,H::T,x::T,sense::T
) where {T <: AbstractFloat}
    if T <: Union{Float32,Float64} && precision(BigFloat)>=128
        value=earth_potential_coefficient(Val(:Maaouni2001),Val(:kernel),
            Complex{BigFloat}(N),BigFloat(k0),BigFloat(H),BigFloat(x),BigFloat(sense))
        return Complex{T}(value)
    end
    A=sqrt(1-N); b=complex(zero(T),-sense)/sqrt(1+N)
    H*real(A)>abs(x*imag(A)) || throw(DomainError((H,x),
        ":Maaouni2001 requires positive real parts for both transformed separation arguments"))
    result=zero(N)
    for z in (k0*complex(H,x),k0*complex(H,-x))
        shift=z+2/A
        result+=2log1p(2/(A*z))+
            scaled_expint_negative(b*shift)+scaled_expint_negative(-b*shift)-
            (N-1)/(N+1)*(scaled_expint_negative(b*z)+scaled_expint_negative(-b*z))
    end
    return result/4
end

function earth_potential_coefficient(::Val{:Maaouni2001},::Val{:mutual},functor,pair)
    _require(pair,Val(:overhead))
    state=functor.state
    H=sum(pair.heights)
    self=pair.row==pair.column
    x=self ? zero(pair.separation) : abs(pair.separation)
    d=self ? pair.separation : hypot(x,pair.heights[1]-pair.heights[2])
    all(isfinite,(H,x,d)) && minimum(pair.heights)>0 && d>0 ||
        throw(DomainError((pair.heights,pair.separation),":Maaouni2001 requires positive overhead geometry"))
    N=state.epsilon[2]/state.epsilon[1]+state.sigma[2]/(state.jω*state.epsilon[1])
    k0=abs(imag(state.jω))*sqrt(state.mu[1]*state.epsilon[1])
    G=earth_potential_coefficient(Val(:Maaouni2001),Val(:kernel),N,k0,H,x,sign(imag(state.jω)))
    πT=one(H)*π
    return _complex_result(state.jω,(log(hypot(H,x)/d)+G)/(2πT*state.epsilon[1]))
end

:Maaouni2001
