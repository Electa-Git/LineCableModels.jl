"""
$(TYPEDEF)

Bound the absolute remainder of a homogeneous-half-space spectral kernel on
``λ=t\\exp(iθ)``. For ``t≥L`` the outgoing roots satisfy
``|a_m-λ|≤β_m/t``, where
``β_m=|κ_m²|/(1+\\sqrt{1-|κ_m²|/L²})``. Triangle inequalities then bound the
material denominators and root ratios. Each exponential envelope is integrated
using ``E_1(dL)≤\\exp(-dL)/(dL)``; the bound is valid only after its root and
denominator conditions hold. All coefficients are copied from the physical
state. Uncertain inputs retain full-range numerical verification.

$(TYPEDFIELDS)
"""
struct EarthSpectralTail{Kind, P, Q, T} <: AbstractSpectralTailBound
    "Squared-root magnitudes |κ²| \\[1/m²\\]."
    k2::NTuple{2,T}
    "Permeability magnitudes \\[H/m\\]."
    mu::NTuple{2,T}
    "Magnitude of the summed permeability \\[H/m\\]."
    musum::T
    "Admittivity magnitudes \\[S/m\\]."
    sh::NTuple{2,T}
    "Magnitude of the summed admittivity \\[S/m\\]."
    shsum::T
    "Receiver and source root-distance coefficients \\[m\\]."
    heights::NTuple{2,T}
    "Explicit exponential weight length, including path padding \\[m\\]."
    height::T
    "Horizontal separation and Bessel radius \\[m\\]."
    extent::NTuple{2,T}
    "Contour angle \\[rad\\]."
    angle::T
    "Column/receiver logarithmic normalization \\[dimensionless\\]."
    logscale::T
    "Receiver I₀ magnitude for a combined voltage path \\[dimensionless\\]."
    i0::T
end

spectral_cacheable_tail(tail) = tail===nothing
spectral_cacheable_tail(::EarthSpectralTail) = true

function earth_tail(::Val{Kind}, ::Val{P}, ::Val{Q}, u, hp, hq, height,
        separation, radius, angle, logscale, i0 = 1) where {Kind,P,Q}
    values=(u.k2..., u.mu..., u.sh..., hp, hq, height, separation, radius, logscale, i0)
    # A nominal bound cannot bound an uncertain derivative. Do not discard a
    # physical correlation merely to make the tail metadata concrete.
    all(x->typeof(real(x))<:Union{AbstractFloat,Integer}, values) || return nothing
    T=promote_type(map(x->typeof(float(real(x))), values)...)
    return EarthSpectralTail{Kind,P,Q,T}(Tuple(T.(abs.(u.k2))), Tuple(T.(abs.(u.mu))), T(abs(sum(u.mu))),
        Tuple(T.(abs.(u.sh))), T(abs(sum(u.sh))), (T(hp),T(hq)), T(height),
        (T(separation),T(radius)), T(angle), T(real(logscale)), T(abs(i0)))
end

function spectral_kernel_tail(kernel::EarthSpectrum{Kind,P,Q}, w, angle) where {Kind,P,Q}
    g=kernel.geometry
    return earth_tail(Val(Kind), Val(P), Val(Q), kernel.state, g.hp, g.hq,
        w.height, w.separation, get(w,:radius,zero(w.height)), angle, g.logscale)
end

function spectral_kernel_tail(kernel::EarthRadialSpectrum{Kind}, w, angle) where {Kind}
    return earth_tail(Val(Kind), Val(2), Val(2), kernel.state, zero(w.height),
        w.height, w.height, w.separation, zero(w.height), angle, kernel.logscale)
end

function spectral_kernel_tail(kernel::EarthPathVoltageSpectrum{P,Q,Reference}, w, angle) where {P,Q,Reference}
    g=kernel.geometry
    kind=Reference===:deep ? Val(:path_deep) : Val(:path_interface)
    return earth_tail(kind, Val(P), Val(Q), kernel.state, g.hp, g.hq,
        w.height+g.padding, w.separation, g.radius, angle, g.logscale, 1+g.i0minus)
end

function spectral_kernel_tail(kernel::RadializedEarthSpectrum, w, angle)
    w.height==kernel.height || return nothing
    return spectral_kernel_tail(kernel.kernel, (height=w.height,separation=w.separation), angle)
end

function (tail::EarthSpectralTail{Kind,P,Q,T})(limit::Real) where {Kind,P,Q,T}
    L=T(limit)
    L>0 && isfinite(L) || return T(Inf)
    eta=map(x->x/L^2, tail.k2)
    maximum(eta)<=T(1)/4 || return T(Inf)
    beta=map((x,e)->x/(1+sqrt(1-e)), tail.k2, eta)
    delta=map(x->x/L^2, beta)
    c,s=cos(tail.angle),sin(tail.angle)
    # This also makes the near-λ square root the outgoing one.
    maximum(delta)<c/2 || return T(Inf)
    dm=tail.musum-tail.mu[2]*delta[1]-tail.mu[1]*delta[2]
    ds=tail.shsum-tail.sh[2]*delta[1]-tail.sh[1]*delta[2]
    dm>0 && ds>0 || return T(Inf)
    lo=map(x->1-x,delta); hi=map(x->1+x,delta)
    hp,hq=tail.heights
    y,r=tail.extent
    correction=tail.logscale+(hp*beta[P]+hq*beta[Q])/L
    d=tail.height*c-(y+r)*s
    d>0 || return T(Inf)
    coefficient=if Kind===:Z
        prod(tail.mu)/dm
    elseif Kind===:phi
        (tail.mu[1]*hi[1]+tail.mu[2]*hi[2])/(dm*ds)
    elseif Kind===:voltage
        hi[3-P]/(lo[P]*ds)
    elseif Kind===:endpoint
        sum(tail.k2)/(lo[1]*lo[2]*ds*L^2)
    elseif Kind===:surface
        hi[1]/(lo[2]*ds)
    elseif Kind===:finite_direct
        inv(2tail.sh[2]*lo[2])
    elseif Kind===:finite_image
        (tail.sh[2]*hi[1]+tail.sh[1]*hi[2])/(2tail.sh[2]*lo[2]*ds)
    else
        # Bound the receiver endpoint and the interface/deep endpoint
        # separately; they have different geometric decay lengths.
        ratio=hi[3-P]/lo[P]
        first_rate=(tail.height+hp)*c-y*s
        first=earth_exponential_tail(correction, tail.i0*ratio/ds, first_rate, L)
        second_ratio=ratio+(Kind===:path_deep ? hi[1]/lo[2] : zero(T))
        second=earth_exponential_tail(tail.logscale+hq*beta[Q]/L, second_ratio/ds, d, L)
        return first+second
    end
    return earth_exponential_tail(correction, coefficient, d, L)
end

function earth_exponential_tail(logscale, coefficient, decay, limit)
    iszero(coefficient) && return zero(limit)
    value=exp(logscale+log(coefficient)-decay*limit-log(decay*limit))
    # Retain a positive representable upper envelope when exp underflows.
    return max(nextfloat(zero(limit)), value*(1+32eps(typeof(limit))))
end
