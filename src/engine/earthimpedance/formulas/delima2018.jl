function routes(identifier::Val{:DeLima2018})
    return (
        self=FormulaMethod(identifier,earth_impedance,Val(:self)),
        mutual=FormulaMethod(identifier,earth_impedance,Val(:mutual)),
        Γ=FormulaMethod(identifier,propagation_constant)
    )
end
assumptions(::Val{:DeLima2018})=(
    air=_full,earth=_full,permeability=vacuum_permeability,approximation=:image
)
propagation(::Val{:DeLima2018})=Val(:explicit)

function Formula(::Val{:DeLima2018};approximation::Symbol=:image,kwargs...)
    approximation in (:image,:quasi_full_wave) ||
        throw(ArgumentError(":DeLima2018 approximation must be :image or :quasi_full_wave"))
    identifier=Val(:DeLima2018); defaults=routes(identifier); overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) ||
        throw(ArgumentError("unknown routes for earth-impedance formula :DeLima2018"))
    selected=merge(defaults,overrides)
    values=merge(assumptions(identifier),(;approximation))
    return Formula{:DeLima2018,typeof(selected),typeof(values)}(selected,values)
end


"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinitely long circular conductor of radius ``r`` at signed region-1 coordinate ``h``; buried application is bare. |
| Calculated quantities | Full-wave parent per-unit-length impedance and noniterative quasi-full-wave evaluation for one overhead or bare buried conductor parallel to a planar interface |
| Earth structure | Two homogeneous half-spaces separated by a plane interface. |
| Model and approximation | The qFW operation is the source's substitution of an image-approximation propagation constant ``\\bar\\gamma`` into ``u_i`` and ``\\eta_1``. It avoids Newton iteration of the full-wave modal equation (7). It is distinct from qTEM, which sets ``\\gamma\\approx\\gamma_1`` and uses a leading logarithmic term, and from closed-form image evaluation of the Sommerfeld integrals. |
| Main source | A. C. S. de Lima, A. P. C. Magalhães, P. E. D. Rocha, R. A. Meyberg, and M. T. C. de Barros (2018) |
| Citation key(s) | `:DeLima2018` |
| Evidence status | Original publication page images checked |

**Expression.** Appendix (23) gives the external image coefficient
``Z_{image}=j\\omega\\mu_0[\\ln(2h/r)+\\bar S_1-\\bar S_2-\\bar S_4]/(2\\pi)``.
For a prescribed source propagation constant ``\\gamma_p``, equation (12)
gives
``Z_e=j\\omega\\mu_0[\\Lambda+S_1-(\\gamma_p^2/\\gamma_1^2)(S_2+S_4)]/(2\\pi)``.
The kernel uses ``u_i^2=\\lambda^2+\\gamma_i^2+k_x^2``,
``k_x=j\\gamma_p``, and ``\\eta_1^2=\\gamma_1^2+k_x^2``
from definition (2), not the inconsistent radial expression in (15).
The exact zero-radial Bessel difference is ``\\ln(D/r)``.

**Image expression.** The default `approximation=:image` evaluates
Appendix D, equations (23)–(24), which are printed in the same paper.
The separate conductor impedance is supplied by the internal family.
For the medium containing the wire, ``\\kappa_1=\\sigma_1+j\\omega\\varepsilon_1``,
``N=\\gamma_2^2/\\gamma_1^2``, ``\\beta=\\sqrt{\\gamma_2^2-\\gamma_1^2}``,
``D=\\sqrt{4h^2+r^2}``, and ``A=(N+1)/(\\beta D)``:

```math
\\begin{aligned}
\\bar S_1&=\\ln\\left(1+\\frac{2}{\\beta D}\\right),\\\\
\\bar S_2&=\\frac{2}{N+1}\\ln(1+A),\\\\
\\bar S_4&=\\frac{2\\ln2}{N+1}
+\\frac{2N}{N+1}\\ln\\left(1+\\frac{1}{1+2A}\\right).
\\end{aligned}
```

The last expression is algebraically (24), arranged to retain the
cancellation at large earth-to-air material ratios.

**Propagation input.** Select `approximation=:quasi_full_wave`
with an explicit backend input `k_x=j*gamma_image`. The image estimate is
``\\gamma_{image}=\\sqrt{(z_i+Z_{image})j\\omega/P_{image}}``.
The source's internal impedance ``z_i`` must be supplied explicitly
when forming that estimate; it is not silently set to zero. The helper
`EarthImpedance.propagation_constant(Val(:DeLima2018), Val(:image),
parameters, z_i, jomega)` returns the backend wavenumber and its square.
Here `parameters` contains the image `impedance` and `potential`
coefficients returned by the shared `Val(:parameters)` route.

The qFW selection evaluates the integral coefficients at that prescribed
wavenumber. It does not solve the full-wave modal equation, select between
fast and TL modes, or claim that the image estimate satisfies the exact
modal equation. A missing qFW propagation input is rejected.


**Scope.** Only a single bare thin wire above or below a planar interface
is covered; media are exchanged for the buried case. Both half-spaces
are homogeneous and nonmagnetic. Mutual coefficients are not supplied.
The source's conductor-loss approximation requires
``|\\gamma_c|\\gg|\\gamma_p|``. Internal impedance remains separate
in normal matrix assembly and is not counted twice.

**Reference.** [DeLima2018](@cite), equations (2)–(4), (7), (10)–(16),
and Appendix D, equations (23)–(24).
"""
description(::Formula{:DeLima2018}) =
    "De Lima et al. single-wire image and prescribed qFW impedance (2018)"

function propagation_constant(::Val{:DeLima2018},s,mu,epsilon)
    return (Γ=zero(s),squared=zero(s))
end

function propagation_constant(::Val{:DeLima2018},::Val{:image},parameters,internal,s)
    physical=sqrt((internal+parameters.impedance)*s/parameters.potential)
    real(physical)<0 && (physical=-physical)
    kx=complex(zero(real(s)),one(real(s)))*physical
    return (Γ=kx,squared=kx^2)
end

function (formula::Formula{:DeLima2018})(rho,epsilon,mu,s,Γ,segments=nothing)
    earth_impedance(Val(:DeLima2018),Val(:validate),formula,rho,epsilon,mu,s,Γ)
    return _homogeneous_functor(Val(:DeLima2018),formula,rho,epsilon,mu,s,Γ,segments)
end

function earth_impedance(::Val{:DeLima2018},::Val{:validate},formula,rho,epsilon,mu,s,Γ)
    length(rho)==length(epsilon)==length(mu)==2 ||
        throw(DimensionMismatch(":DeLima2018 requires two homogeneous half-spaces"))
    all(x->(isinf(x)&&x>0)||(isfinite(x)&&x>0),rho) &&
        all(x->isfinite(x)&&x>0,epsilon) ||
        throw(DomainError((rho,epsilon),":DeLima2018 requires positive material constants"))
    all(x->isapprox(x,vacuum_permeability(x)),mu) ||
        throw(DomainError(mu,":DeLima2018 assumes nonmagnetic media"))
    isfinite(s) && !iszero(s) && iszero(real(s)) ||
        throw(DomainError(s,":DeLima2018 requires nonzero real frequency"))
    if formula.assumptions.approximation===:quasi_full_wave
        Γ===nothing && throw(ArgumentError(":DeLima2018 qFW requires an explicit image-derived wavenumber"))
    else
        Γ===nothing || iszero(Γ) ||
            throw(ArgumentError(":DeLima2018 image coefficients do not use an independent wavenumber"))
    end
    return nothing
end

function earth_impedance(::Val{:DeLima2018},::Val{:self},functor,pair)
    return earth_impedance(Val(:DeLima2018),Val(:parameters),functor,pair).impedance
end

function earth_impedance(::Val{:DeLima2018},::Val{:mutual},functor,pair)
    throw(ArgumentError(":DeLima2018 supplies only single-wire self coefficients"))
end

function earth_impedance(::Val{:DeLima2018},::Val{:parameters},functor,pair)
    pair.row==pair.column && pair.layers[1]==pair.layers[2] &&
        pair.layers[1] in (1,2) ||
        throw(ArgumentError(":DeLima2018 supplies only single-wire self coefficients"))
    state=functor.state; s=state.jω
    source=pair.layers[1]; other=source==1 ? 2 : 1
    h=abs(pair.heights[1]); r=pair.separation; D=hypot(2h,r)
    all(isfinite,(h,r)) && h>0 && r>0 ||
        throw(DomainError((h,r),":DeLima2018 requires positive interface distance and radius"))
    g1=state.gamma_medium_squared[source]; g2=state.gamma_medium_squared[other]
    kappa=state.sigma[source]+s*state.epsilon[source]
    πT=one(h)*π
    if state.formula.assumptions.approximation===:image
        beta=sqrt(g2-g1)
        iszero(beta) && throw(DomainError(beta,":DeLima2018 image form requires nonzero medium contrast"))
        N=g2/g1; A=(N+1)/(beta*D)
        S1=log1p(2/(beta*D))
        S2=2log1p(A)/(N+1)
        S4=2log(one(h)*2)/(N+1)+2N/(N+1)*log1p(inv(1+2A))
        ideal=log(2h/r)
        impedance=s*state.mu[source]/(2πT)*(ideal+S1-S2-S4)
        potential=s/(2πT*kappa)*(ideal-S4)
    else
        source_squared=g1+state.gamma_squared
        other_squared=g2+state.gamma_squared
        radial=spectral_root(source_squared,s)
        Lambda=iszero(radial) ? log(D/r) :
            special_besselk(0,radial*r)-special_besselk(0,radial*D)
        ratio=-state.gamma_squared/g1; N=g2/g1
        terms=quadgk(zero(h),oftype(h,Inf);rtol=state.tolerance) do t
            lambda=t/h
            q1=spectral_root(lambda^2+source_squared,s)
            q2=spectral_root(lambda^2+other_squared,s)
            decay=exp(-h*q1)
            quotient=iszero(q1) ? complex(h) : -expm1(-h*q1)/q1
            common=N*q1+q2
            electric=2q2*decay*quotient/common
            combined=(N==1 && ratio==1) ? zero(s) :
                2decay^2*((N-ratio)*q1+(1-ratio)*q2)/((q1+q2)*common)
            weight=cos(lambda*r)/h
            [(combined-ratio*electric)*weight,electric*weight]
        end[1]
        impedance=s*state.mu[source]/(2πT)*(Lambda+terms[1])
        potential=s/(2πT*kappa)*(Lambda-terms[2])
    end
    return (
        impedance=_complex_result(s,impedance),
        potential=_complex_result(s,potential)
    )
end

:DeLima2018
