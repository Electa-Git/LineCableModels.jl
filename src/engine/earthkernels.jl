"Return electrical conductivity, including the open-circuit limit."
@inline conductivity(rho) = isinf(rho) ? zero(rho) : inv(rho)

@inline special_besselix(order::Integer, value) = SpecialFunctions.besselix(order, value)
@inline special_besselkx(order::Integer, value) = SpecialFunctions.besselkx(order, value)
@inline special_besselk(order::Integer, value) = SpecialFunctions.besselk(order, value)
@inline special_besseljx(order::Integer, value) = SpecialFunctions.besseljx(order, value)
@inline special_besselix(order::Integer,
    value::Complex{Float32}) = ComplexF32(SpecialFunctions.besselix(Float32(order), value))
@inline special_besselkx(order::Integer,
    value::Complex{Float32}) = ComplexF32(SpecialFunctions.besselkx(Float32(order), value))
@inline special_besselk(order::Integer,
    value::Complex{Float32}) = ComplexF32(SpecialFunctions.besselk(Float32(order), value))
@inline special_besseljx(order::Integer,
    value::Complex{Float32}) = ComplexF32(SpecialFunctions.besseljx(Float32(order), value))
function special_besseljx(order::Integer, value::Complex{BigFloat})
    order>=0 || throw(DomainError(order, "Bessel order must be nonnegative"))
    abs(value)>max(128, precision(BigFloat), order^2) ||
        return (-im)^order*special_besselix(order, im*value)
    real(value)<0 && return (-1)^order*special_besseljx(order, -value)
    # DLMF 10.17.1, 10.17.5–6: combine the two scaled Hankel expansions.
    # This avoids allocating an angular grid proportional to a DE tail node.
    plus=one(value)
    minus=one(value)
    tp=one(value)
    tm=one(value)
    for n in 1:100_000
        factor=(4BigFloat(order)^2-BigFloat(2n-1)^2)/(8n*value)
        tp*=im*factor
        tm*=-im*factor
        plus+=tp
        minus+=tm
        max(abs(tp), abs(tm))<=eps(BigFloat)*max(abs(plus), abs(minus)) && break
        n==100_000 &&
            throw(ErrorException("scaled BigFloat Bessel asymptotic did not converge"))
    end
    phase=cis(BigFloat(π)*(BigFloat(order)/2+BigFloat(1)/4))
    growth=abs(imag(value))
    return (exp(im*value-growth)*conj(phase)*plus +
            exp(-im*value-growth)*phase*minus)/sqrt(2BigFloat(π)*value)
end

# SpecialFunctions omits complex BigFloat Bessel functions. The local methods
# retain the caller's working precision for unsupported argument types.
function special_besselix(order::Integer, value::Complex{BigFloat})
    order >= 0 || throw(DomainError(order, "Bessel order must be nonnegative"))
    abs(value)>max(128, precision(BigFloat), order^2) &&
        return (-im)^order*special_besseljx(order, im*value)
    if abs(value)>8
        # DLMF 10.32.3, with the scaling inside the exponential. Seed the
        # angular phase scale rather than relying on two aliased rules.
        count=max(8, ceil(Int, abs(imag(value))+order))
        points=collect(range(zero(BigFloat), BigFloat(π); length = count+1))
        f=θ->exp(value*cos(θ)-abs(real(value)))*cos(order*θ)
        result, error=quadgk(f, points; rtol = sqrt(eps(BigFloat)))
        return result/BigFloat(π)
    end
    half = value / BigFloat(2)
    term = half^order / BigFloat(factorial(big(order)))
    result = term
    for index in 1:100_000
        term *= half^2 / (BigFloat(index) * BigFloat(index + order))
        next = result + term
        if next == result || abs(term) <= eps(BigFloat) * max(abs(next), one(BigFloat))
            return exp(-abs(real(value))) * next
        end
        result = next
    end
    throw(ErrorException("complex BigFloat besseli series did not converge"))
end

function special_besselk(order::Integer, value::Complex{BigFloat})
    return exp(-value)*special_besselkx(order, value)
end

function special_besselkx(order::Integer, value::Complex{BigFloat})
    order>=0 || throw(DomainError(order, "Bessel order must be nonnegative"))
    !iszero(value)&&real(value)>=0 || throw(DomainError(value,
        "complex BigFloat K requires a nonzero argument on the outgoing right half-plane"))
    # DLMF 10.32.8, w=z(t−1), rotate the w contour to the positive real
    # axis, then w=u². This remains exponentially decaying in the lossless
    # imaginary-argument limit; exp(-z*cosh(t)) does not.
    exponent=BigFloat(order)-BigFloat(1)/2
    f=u->exp(-u*u)*u^(2order)*(1+u*u/(2value))^exponent
    feature=min(sqrt(abs(value)), BigFloat(1)/2)
    result,
    error=quadgk(f, zero(BigFloat), feature, one(BigFloat), BigFloat(Inf);
        rtol = sqrt(eps(BigFloat)))
    return 2sqrt(BigFloat(π)/(2value))/SpecialFunctions.gamma(BigFloat(order)+BigFloat(1)/2)*result
end

"Evaluate the shared earth-return modified-Bessel difference."
@inline function bessel_difference(gamma_s, inner, outer)
    maximum_argument = max(abs(gamma_s) * inner, abs(gamma_s) * outer)
    isapprox(nominal(maximum_argument), 0; atol = 1.0e-6) &&
        return log(outer / inner)
    return special_besselk(0, gamma_s * inner) -
           special_besselk(0, gamma_s * outer)
end

"""
$(TYPEDEF)

Store the exterior circles used by the complete earth-current constraint.
Rows of the resulting matrices are receivers; columns are sources.

$(TYPEDFIELDS)
"""
struct EarthReturnGeometry{T <: Real}
    "Horizontal conductor centres \\[m\\]."
    horizontal::Vector{T}
    "Signed conductor heights above the interface \\[m\\]."
    height::Vector{T}
    "Exterior conductor or cable-port radii \\[m\\]."
    radius::Vector{T}

    function EarthReturnGeometry{T}(horizontal, height, radius) where {T <: Real}
        length(horizontal) == length(height) == length(radius) > 0 ||
            throw(DimensionMismatch("exterior centres and radii must align"))
        all(isfinite, horizontal) && all(isfinite, height) && all(isfinite, radius) ||
            throw(DomainError((horizontal, height, radius), "exterior geometry must be finite"))
        for p in eachindex(radius)
            zero(T) < radius[p] < abs(height[p]) || throw(DomainError(radius[p],
                "each exterior circumference must lie wholly in one half-space"))
            for q in 1:(p - 1)
                hypot(horizontal[p]-horizontal[q], height[p]-height[q]) >
                radius[p]+radius[q] ||
                    throw(DomainError((p, q), "exterior circumferences must not overlap"))
            end
        end
        return new{T}(collect(horizontal), collect(height), collect(radius))
    end
end

function EarthReturnGeometry(horizontal::AbstractVector{T}, height::AbstractVector{T},
        radius::AbstractVector{T}) where {T <: Real}
    return EarthReturnGeometry{T}(horizontal, height, radius)
end

"Select the outgoing square root for the exp(jωt) convention."
@inline function outgoing_root(z)
    root = sqrt(complex(z))
    return real(root) < 0 || (iszero(real(root)) && imag(root) < 0) ? -root : root
end

@inline root_difference(k2, a, λ) = iszero(k2) ? zero(a) : k2/(a+λ)

"Evaluate I₀(z)−1 without subtracting two nearly equal numbers."
function bessel_i0m1(z)
    if abs(nominal(z)) > 0.5
        return special_besselix(0, z)*exp(abs(real(z)))-one(z)
    end
    term=z*z/4
    result=term
    for n in 2:1000
        term *= (z/(2n))^2
        next=result+term
        next == result && return next
        result=next
    end
    return result
end

@inline bessel_j0m1(z) = bessel_i0m1(im*z)

# Combine receiver and reference endpoints before numerical integration.
# The interface kernel is other_root*(I0*exp(-hp*ap)-J0)/(ap*ds).
# For an air receiver at deep reference, add a0*J0/(ag*ds).
struct EarthPathVoltageSpectrum{P, Q, Reference, S, G}
    state::S
    geometry::G
end

function earth_weighted_spectrum(kernel::EarthPathVoltageSpectrum{P, Q, Reference},
        λ, logweight) where {P, Q, Reference}
    u=kernel.state
    g=kernel.geometry
    a0=outgoing_root(λ^2+u.k2[1])
    ag=outgoing_root(λ^2+u.k2[2])
    ap=P==1 ? a0 : ag
    aq=Q==1 ? a0 : ag
    other=P==1 ? ag : a0
    decay=exp(g.logscale-g.hq*root_difference(u.k2[Q], aq, λ)+logweight)
    j0, difference=if iszero(ap)
        argument=g.radius*λ
        j0=exp(-g.padding*λ+abs(imag(argument)))*special_besseljx(0, argument)
        (j0, -g.hp*j0)
    elseif abs(nominal(g.radius*λ))<0.5
        jminus=bessel_j0m1(g.radius*λ)
        envelope=exp(-g.padding*λ)
        (envelope*(one(jminus)+jminus),
            envelope*(g.i0minus*exp(-g.hp*ap)+expm1(-g.hp*ap)-jminus)/ap)
    else
        argument=g.radius*λ
        envelope=exp(-g.padding*λ+abs(imag(argument)))
        scaled=iszero(spectral_magnitude(envelope)) ? zero(envelope) :
               envelope*special_besseljx(0, argument)
        (scaled, ((one(g.i0minus)+g.i0minus)*exp(-g.padding*λ-g.hp*ap)-scaled)/ap)
    end
    numerator=other*difference
    Reference===:deep && (numerator+=a0/ag*j0)
    return decay*numerator/(u.sh[2]*a0+u.sh[1]*ag)
end
(kernel::EarthPathVoltageSpectrum)(λ) = earth_weighted_spectrum(kernel, λ, zero(λ))
function earth_combined_weight(kernel::EarthPathVoltageSpectrum)
    g=kernel.geometry
    return nominal(g.logscale+g.hq*maximum(abs, kernel.state.k))>300
end

"Evaluate the removable I₁(z)/(z I₀(z)) limit."
@inline function bessel_current_ratio(z)
    if abs(nominal(z)) < 1e-3
        z2=z*z
        return one(z)/2-z2/16+z2*z2/96-11z2^3/6144
    end
    return special_besselix(1, z)/(z*special_besselix(0, z))
end

"""
$(TYPEDSIGNATURES)

Prepare the manuscript medium roots and boundary-current factors at fixed
frequency. Source columns are scaled by exp(abs(real(κ r))) before assembly;
this common scaling cancels from the final matrix solves and prevents overflow
in products of growing I₀ and decaying exterior fields.

# Returns

- Typed medium data and per-conductor factors for K, H and L.
"""
function unified_earth_state(state, geometry::EarthReturnGeometry, arrays = nothing)
    s=state.jω
    sh=ntuple(m->state.sigma[m]+s*state.epsilon[m], 2)
    mu=ntuple(m->state.mu[m], 2)
    k2=ntuple(m->state.gamma_medium_squared[m]-state.Γ^2, 2)
    k=map(outgoing_root, k2)
    n=length(geometry.radius)
    if arrays===nothing
        T=eltype(geometry.radius)
        arrays=(x = Vector{Complex{T}}(undef, n), scaling = Vector{T}(undef, n),
            A = Vector{Complex{T}}(undef, n), F = Vector{Complex{T}}(undef, n))
    end
    x, scaling, A, F=arrays.x, arrays.scaling, arrays.A, arrays.F
    for p in eachindex(x)
        medium=geometry.height[p]>0 ? 1 : 2
        x[p]=k[medium]*geometry.radius[p]
        scaling[p]=abs(real(x[p]))
        A[p]=special_besselix(0, x[p])
        F[p]=2*(one(s)*π)*sh[medium]*geometry.radius[p]^2*bessel_current_ratio(x[p])
    end
    return (; s, Γ = state.Γ, sh, mu, k2, k, x, scaling, A, F)
end

# Each term excludes exp(-height*λ), cos(yλ), and optional J₀(rλ).
# Kind and the two ordered media dispatch outside the spectral loop.
struct EarthSpectrum{Kind, Receiver, Source, S, G}
    state::S
    geometry::G
end

@inline function earth_weighted_spectrum(
        kernel::EarthSpectrum{
            Kind, P, Q}, λ, logweight) where {Kind, P, Q}
    u=kernel.state
    g=kernel.geometry
    a0=outgoing_root(λ^2+u.k2[1])
    ag=outgoing_root(λ^2+u.k2[2])
    a=(a0, ag)
    dm=u.mu[2]*a0+u.mu[1]*ag
    ds=u.sh[2]*a0+u.sh[1]*ag
    # Rationalized root differences preserve the large-λ exponential.
    decay=exp(g.logscale-g.hp*root_difference(u.k2[P], a[P], λ) -
              g.hq*root_difference(u.k2[Q], a[Q], λ)+logweight)
    if Kind === :Z
        return decay*u.mu[1]*u.mu[2]/dm
    elseif Kind === :phi
        return decay*(u.mu[1]*a0+u.mu[2]*ag)/(dm*ds)
    elseif Kind === :voltage
        return decay*(P==2 ? a0/ag : ag/a0)/ds
    elseif Kind === :endpoint
        # a0/ag-ag/a0 = (κ0²-κg²)/(a0*ag), evaluated without cancellation.
        return decay*(u.k2[1]-u.k2[2])/(a0*ag*ds)
    elseif Kind === :surface
        return decay*a0/(ag*ds)
    elseif Kind === :finite_direct
        return decay/(2u.sh[2]*ag)
    elseif Kind === :finite_image
        return decay*(u.sh[2]*a0-u.sh[1]*ag)/(2u.sh[2]*ag*ds)
    end
    throw(ArgumentError("unknown earth spectral kernel"))
end

@inline (kernel::EarthSpectrum)(λ) = earth_weighted_spectrum(kernel, λ, zero(λ))
function earth_combined_weight(kernel::EarthSpectrum)
    g=kernel.geometry
    return nominal(g.logscale+(g.hp+g.hq)*maximum(abs, kernel.state.k))>300
end

function earth_spectrum(
        ::Val{Kind}, ::Val{P}, ::Val{Q}, state, hp, hq, logscale) where {Kind, P, Q}
    g=(; hp, hq, logscale)
    return EarthSpectrum{Kind, P, Q, typeof(state), typeof(g)}(state, g)
end

# Evaluate the averaged direct-image term with the source-column scaling.
function earth_direct(state, geometry, p, q, medium)
    u=state
    r=geometry.radius[p]
    h=abs(geometry.height[p])+abs(geometry.height[q])
    d=hypot(geometry.horizontal[p]-geometry.horizontal[q], h)
    D=p==q ? r :
      hypot(geometry.horizontal[p]-geometry.horizontal[q],
        geometry.height[p]-geometry.height[q])
    k=u.k[medium]
    xp=u.x[p]
    sp=u.scaling[p]
    sq=u.scaling[q]
    if abs(nominal(k*d)) < 0.5
        difference=iszero(k) ? log(d/D) : special_besselk(0, k*D)-special_besselk(0, k*d)
        if p==q
            correction=iszero(k) ? zero(k) : bessel_i0m1(xp)*special_besselk(0, k*d)
            return exp(sq)*(difference-correction)
        end
        return u.A[p]*exp(sp+sq)*difference
    end
    if p==q
        direct=special_besselkx(0, k*r)*exp(sq-k*r)
        image=u.A[p]*special_besselkx(0, k*d)*exp(sp+sq-k*d)
        return direct-image
    end
    return u.A[p]*(special_besselkx(0, k*D)*exp(sp+sq-k*D) -
                   special_besselkx(0, k*d)*exp(sp+sq-k*d))
end

# Same-earth terms have an exact Sommerfeld/Bessel image transform in u=ag.
# Fitting the bounded residual in u avoids fitting the artificial exp(h*λ)
# growth introduced when exp(-h*ag) is split into a cosine-weighted λ kernel.
struct EarthRadialSpectrum{Kind, S, T}
    state::S
    logscale::T
end

# Exact change of analytic weight: f(λ)e^(-hλ) becomes
# [u*f(λ)*e^(h(u-λ))] e^(-hu)/u, u²=λ²+q². This removes
# the air-root singularity from the residual fitted by complex images.
struct RadializedEarthSpectrum{K, Q, H}
    kernel::K
    q::Q
    height::H
end
function (kernel::RadializedEarthSpectrum)(u)
    λ=outgoing_root((u-kernel.q)*(u+kernel.q))
    return radial_kernel_value(kernel, u, λ)
end
function radial_kernel_value(kernel::RadializedEarthSpectrum, u, λ)
    difference=iszero(kernel.q) ? zero(u) : kernel.q^2/(u+λ)
    return u*kernel.kernel(λ)*exp(kernel.height*difference)
end
@inline function (kernel::EarthRadialSpectrum{Kind})(ag) where {Kind}
    u=kernel.state
    a0=outgoing_root((ag-u.k[2])*(ag+u.k[2])+u.k2[1])
    return earth_radial_value(kernel, ag, a0)
end
@inline function radial_kernel_value(kernel::EarthRadialSpectrum, ag, λ)
    return earth_radial_value(kernel, ag, outgoing_root(λ^2+kernel.state.k2[1]))
end
@inline function earth_radial_value(kernel::EarthRadialSpectrum{Kind}, ag, a0) where {Kind}
    u=kernel.state
    dm=u.mu[2]*a0+u.mu[1]*ag
    ds=u.sh[2]*a0+u.sh[1]*ag
    factor=exp(kernel.logscale)
    if Kind===:Z
        return factor*u.mu[1]*u.mu[2]*ag/dm
    elseif Kind===:phi
        return factor*ag*(u.mu[1]*a0+u.mu[2]*ag)/(dm*ds)
    else
        return factor*a0/ds
    end
end

# AMOS cannot reduce enormous complex arguments generated by very short images.
# The scaled large-argument K₀ expansion avoids that unnecessary phase reduction.
function image_besselk0x(z)
    abs(nominal(z))<1000 && return special_besselkx(0, z)
    term=one(z)
    value=term
    R=typeof(float(nominal(real(z))))
    for n in 1:1000
        term *= -((2n-1)^2)/(8n*z)
        next=value+term
        abs(nominal(term))<=eps(R)*abs(nominal(next)) &&
            return sqrt((one(z)*π)/(2z))*next
        value=next
    end
    throw(ErrorException("scaled image K₀ expansion did not converge"))
end

# Sampling decisions use nominal branch positions, including degenerate
# equal-media denominators. Physical kernel evaluations retain uncertainty.
function earth_denominator_pole(state)
    sh=map(nominal, state.sh)
    k2=map(nominal, state.k2)
    denominator=sh[2]^2-sh[1]^2
    iszero(denominator) && return nothing
    return outgoing_root((sh[1]^2*k2[2]-sh[2]^2*k2[1])/denominator)
end
function earth_valid_pole(state, candidate, ::Type{R}) where {R}
    sh=map(nominal, state.sh)
    k2=map(nominal, state.k2)
    a=sh[2]*outgoing_root(candidate^2+k2[1])
    b=sh[1]*outgoing_root(candidate^2+k2[2])
    return abs(a+b)<=sqrt(eps(R))*(abs(a)+abs(b))
end

# The medium and denominator scales are independent; retain both, even when
# their ratio spans many decades. Bridge them without prescribing a fixed grid.
function earth_spectral_features(
        state, height, separation, radius, angle, workspace = nothing)
    R=typeof(float(nominal(real(state.s))))
    scales=spectral_buffer(workspace, R, :scales)
    for k in state.k
        iszero(nominal(k)) || push!(scales, R(abs(nominal(k))))
    end
    for (p, q) in ((1, 2), (2, 1))
        !iszero(nominal(state.sh[q])) &&
            push!(scales, R(abs(nominal(state.k[q]*state.sh[p]/state.sh[q]))))
    end
    decay=R(nominal(height))*cos(angle)-R(nominal(separation+radius))*sin(angle)
    decay>0 || throw(DomainError(angle, "spectral contour has a growing geometric weight"))
    push!(scales, inv(decay))
    filter!(x->isfinite(x)&&x>zero(R), scales)
    spectral_sort_unique!(scales)
    points=spectral_buffer(workspace, R, :features)
    push!(points, zero(R))
    for scale in scales
        append!(points, (scale/2, scale, 2scale))
    end
    rotation=cis(angle)
    for k in state.k, sign in (-1, 1)

        earth_feature_neighbourhood!(points, nominal(sign*im*k)/rotation)
    end
    pole=haskey(state, :k2) ? earth_denominator_pole(state) : nothing
    if pole!==nothing
        for candidate in (pole, -pole)
            earth_valid_pole(state, candidate, R) || continue
            earth_feature_neighbourhood!(points, nominal(candidate)/rotation)
        end
    end
    spectral_sort_unique!(points)
    # Extra intervals only bridge gaps between physically declared scales.
    seeds=spectral_buffer(workspace, R, :seeds)
    append!(seeds, points)
    for i in 2:(length(seeds) - 1)
        x=4seeds[i]
        while x<seeds[i + 1]
            push!(points, x)
            x*=4
        end
    end
    spectral_sort_unique!(points)
    return SpectralFeatures{R, Nothing}(points, nothing)
end

function earth_feature_neighbourhood!(points::Vector{R}, location) where {R}
    centre=R(real(location))
    width=R(abs(imag(location)))
    centre>0 || return points
    width=max(width, 64eps(R)*centre)
    for offset in (-4, -1, 0, 1, 4)
        point=centre+offset*width
        point>0&&isfinite(point) && push!(points, point)
    end
    return points
end

# A prescribed Γ can move a branch point into the first quadrant. Keep the
# contour below it; the usual passive Γ=0 roots do not restrict the upper ray.
function earth_contour_angle(state, proposed)
    result=proposed
    for k in state.k, sign in (-1, 1)

        branch=nominal(sign*im*k)
        real(branch)>0&&imag(branch)>0 || continue
        result=min(result, atan(imag(branch), real(branch))/2)
    end
    # Squaring the electric denominator supplies pole candidates. Only a
    # candidate satisfying the original, unsquared denominator is a pole.
    pole=earth_denominator_pole(state)
    if pole!==nothing
        for candidate in (pole, -pole)
            real(candidate)>0&&imag(candidate)>0 || continue
            earth_valid_pole(state, candidate, typeof(proposed)) || continue
            result=min(result, atan(imag(candidate), real(candidate))/2)
        end
    end
    return typeof(proposed)(result)
end

function retained_earth_features(state, geometry, workspace = nothing)
    angle=min(pi/4, atan(float(nominal(geometry.H/(2geometry.y_ij)))))
    medium=(s = state.jω, k = state.gamma, k2 = state.gamma_medium_squared,
        sh = ntuple(m->state.sigma[m]+state.jω*state.epsilon[m], 2))
    features=earth_spectral_features(
        medium, geometry.H, geometry.y_ij, zero(geometry.H), angle, workspace)
    contrast=abs(nominal(sqrt(state.gamma_medium_squared[2]-state.gamma_medium_squared[1])))
    iszero(contrast) || append!(features.points, (contrast/2, contrast, 2contrast))
    return features
end
