"""
Select the outgoing square root for the exp(jωt) convention.
"""
@inline function outgoing_root(z)
    root = sqrt(complex(z))
    return real(root) < 0 || (iszero(real(root)) && imag(root) < 0) ? -root : root
end

@inline root_difference(k2, a, λ) = iszero(k2) ? zero(a) : k2/(a+λ)

# Keep the decay and oscillation in the same exponential on a complex ray.
@inline earth_cosine(height, separation, λ) = (exp(-(height+im*separation)*λ)+exp(-(height-im*separation)*λ))/2

function earth_spectral_value(kernel, λ, height, separation, radius)
    if iszero(radius)
        return earth_combined_weight(kernel) ?
               (earth_weighted_spectrum(kernel, λ, -(height+im*separation)*λ) +
                earth_weighted_spectrum(kernel, λ, -(height-im*separation)*λ))/2 :
               kernel(λ)*earth_cosine(height, separation, λ)
    end
    growth=abs(imag(radius*λ))
    weighted=if earth_combined_weight(kernel)
        (earth_weighted_spectrum(kernel, λ, -(height+im*separation)*λ+growth) +
         earth_weighted_spectrum(kernel, λ, -(height-im*separation)*λ+growth))/2
    else
        kernel(λ)*(exp(-(height+im*separation)*λ+growth) +
                   exp(-(height-im*separation)*λ+growth))/2
    end
    isfinite(weighted) ||
        throw(DomainError((λ, weighted), "nonfinite weighted earth kernel"))
    iszero(numerical_magnitude(weighted)) && return weighted
    return weighted*special_besseljx(0, radius*λ)
end

"""
Evaluate I₀(z)−1 without subtracting two nearly equal numbers.
"""
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

# Combine the air receiver and its interface endpoint before integration.
# The kernel is ag*(I0*exp(-hp*a0)-J0)/(a0*ds).
struct AirVoltageSpectrum{Q, S, G}
    state::S
    geometry::G
end

function earth_weighted_spectrum(kernel::AirVoltageSpectrum{Q},
        λ, logweight) where {Q}
    u=kernel.state
    g=kernel.geometry
    a0=outgoing_root(λ^2+u.k2[1])
    ag=outgoing_root(λ^2+u.k2[2])
    ap=a0
    aq=Q==1 ? a0 : ag
    decay=exp(g.logscale-g.hq*root_difference(u.k2[Q], aq, λ)+logweight)
    difference=if iszero(ap)
        argument=g.radius*λ
        j0=exp(-g.padding*λ+abs(imag(argument)))*special_besseljx(0, argument)
        -g.hp*j0
    elseif abs(nominal(g.radius*λ))<0.5
        jminus=bessel_j0m1(g.radius*λ)
        envelope=exp(-g.padding*λ)
        envelope*(g.i0minus*exp(-g.hp*ap)+expm1(-g.hp*ap)-jminus)/ap
    else
        argument=g.radius*λ
        envelope=exp(-g.padding*λ+abs(imag(argument)))
        scaled=iszero(numerical_magnitude(envelope)) ? zero(envelope) :
               envelope*special_besseljx(0, argument)
        ((one(g.i0minus)+g.i0minus)*exp(-g.padding*λ-g.hp*ap)-scaled)/ap
    end
    return decay*ag*difference/(u.sh[2]*a0+u.sh[1]*ag)
end
(kernel::AirVoltageSpectrum)(λ) = earth_weighted_spectrum(kernel, λ, zero(λ))
function earth_combined_weight(kernel::AirVoltageSpectrum)
    g=kernel.geometry
    return nominal(g.logscale+g.hq*maximum(abs, kernel.state.k))>300
end

"""
Evaluate the removable I₁(z)/(z I₀(z)) limit.
"""
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
function _unified_state!(arrays, state, geometry)
    s=state.jω
    sh=ntuple(m->state.sigma[m]+s*state.epsilon[m], 2)
    mu=ntuple(m->state.mu[m], 2)
    gamma=ntuple(m->sqrt(s*mu[m]*sh[m]), 2)
    k2=ntuple(m->gamma[m]^2-state.Γ^2, 2)
    k=map(outgoing_root, k2)
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
    elseif Kind === :air_reference
        return -decay*ag/(a0*ds)
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

# Subdivision decisions use nominal branch positions, including degenerate
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
function _unified_points!(arrays, state, height, separation, radius, angle)
    R=typeof(float(nominal(real(state.s))))
    scales=empty!(arrays.scales)
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
    unique!(sort!(scales; alg = Base.Sort.QuickSort))
    points=empty!(arrays.points)
    push!(points, zero(R))
    for scale in scales
        append!(points, (scale/2, scale, 2scale))
    end
    rotation=cis(angle)
    for k in state.k, sign in (-1, 1)

        _unified_neighbourhood!(points, nominal(sign*im*k)/rotation)
    end
    pole=earth_denominator_pole(state)
    if pole!==nothing
        for candidate in (pole, -pole)
            earth_valid_pole(state, candidate, R) || continue
            _unified_neighbourhood!(points, nominal(candidate)/rotation)
        end
    end
    unique!(sort!(points; alg = Base.Sort.QuickSort))
    # Extra intervals only bridge gaps between physically declared scales.
    seeds=empty!(arrays.seeds)
    append!(seeds, points)
    for i in 2:(length(seeds) - 1)
        x=4seeds[i]
        while x<seeds[i + 1]
            push!(points, x)
            x*=4
        end
    end
    return points
end

function _unified_neighbourhood!(points::Vector{R}, location) where {R}
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
