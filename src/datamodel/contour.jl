"""
$(TYPEDEF)

Represent a straight curve between two points in the cross-section plane.

$(TYPEDFIELDS)
"""
struct Segment{T <: Real}
    "Coordinates at `t = 0` in \\[m\\]."
    a::NTuple{2, T}
    "Coordinates at `t = 1` in \\[m\\]."
    b::NTuple{2, T}
    function Segment{T}(a::NTuple{2, T}, b::NTuple{2, T}) where {T <: Real}
        return validate(new{T}(a, b))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a segment from the points `a` and `b` \\[m\\], in the type that their
coordinates promote to with `Float64`.

# Errors

- Throws `DomainError` when a coordinate is not finite or the points coincide.
"""
function Segment(a::Tuple{Real, Real}, b::Tuple{Real, Real})
    T = promote_type(Float64, typeof.(a)..., typeof.(b)...)
    return Segment{T}(convert.(T, a), convert.(T, b))
end

function validate(segment::Segment)
    all(isfinite, (segment.a..., segment.b...)) || throw(DomainError(
        (segment.a, segment.b), "Segment coordinates must be finite"))
    segment.a != segment.b || throw(DomainError(
        (segment.a, segment.b), "Segment endpoints must differ"))
    return segment
end

"""
$(TYPEDEF)

Represent a circular curve from the x-axis of its placement over a signed angular span.
A positive span runs counterclockwise and a negative span clockwise.

$(TYPEDFIELDS)
"""
struct Arc{T <: Real}
    "Center \\[m\\] and initial direction \\[rad\\]."
    at::Placement{T}
    "Radius \\[m\\]."
    r::T
    "Signed angular span \\[rad\\]."
    span::T
    function Arc{T}(at::Placement{T}, r::T, span::T) where {T <: Real}
        return validate(new{T}(at, r, span))
    end
end

"""
$(TYPEDSIGNATURES)

Construct an arc at `at` with radius `r` \\[m\\] and signed angular `span` \\[rad\\].
The placement translates in \\[m\\] and rotates in \\[rad\\].
The arguments promote with `Float64`.

# Errors

- Throws `DomainError` when `r` is not positive and finite, or when `span` is not finite
  or does not satisfy ``0 < |\\mathrm{span}| \\leq 2\\pi``.
"""
function Arc(at::Placement, r::Real, span::Real)
    T = promote_type(Float64, eltype(at), typeof(r), typeof(span))
    return Arc{T}(convert(Placement{T}, at), convert(T, r), convert(T, span))
end

function validate(arc::Arc)
    isfinite(arc.r) && arc.r > zero(arc.r) ||
        throw(DomainError(arc.r, "Arc.r must be positive and finite"))
    isfinite(arc.span) && zero(arc.span) < abs(arc.span) <= 2π ||
        throw(DomainError(arc.span, "Arc.span must be finite with 0 < |span| ≤ 2π"))
    return arc
end

"""
$(TYPEDEF)

Represent a placed elliptic curve with parametric angle ``θ=φ0+s\\,span``.
A positive span runs counterclockwise. The normal offset moves the curve outward.

$(TYPEDFIELDS)
"""
struct EllipticArc{T <: Real}
    "Center \\[m\\] and direction of semi-axis `a` \\[rad\\]."
    at::Placement{T}
    "Semi-axis along the local x-axis \\[m\\]."
    a::T
    "Semi-axis along the local y-axis \\[m\\]."
    b::T
    "Parametric angle at the start \\[rad\\]."
    φ0::T
    "Signed parametric span \\[rad\\]."
    span::T
    "Distance along the outward normal \\[m\\], zero for an ellipse's own contour."
    offset::T
    @doc """
    $(TYPEDSIGNATURES)

    Construct an elliptic arc through `validate`. Its semi-axes and normal `offset` are
    in \\[m\\], and its parametric angles are in \\[rad\\].
    """
    function EllipticArc{T}(
            at::Placement{T}, a::T, b::T, φ0::T, span::T, offset::T) where {T <: Real}
        return validate(new{T}(at, a, b, φ0, span, offset))
    end
end

"""
$(TYPEDSIGNATURES)

Construct an elliptic arc from its placement, semi-axes \\[m\\], initial parametric
angle and signed span \\[rad\\], and normal `offset` \\[m\\] (zero by default),
promoted with `Float64`.
"""
function EllipticArc(
        at::Placement, a::Real, b::Real, φ0::Real, span::Real, offset::Real = 0)
    T = promote_type(
        Float64, eltype(at), typeof(a), typeof(b), typeof(φ0), typeof(span), typeof(offset))
    return EllipticArc{T}(convert(Placement{T}, at), convert(T, a), convert(T, b),
        convert(T, φ0), convert(T, span), convert(T, offset))
end

"""
$(TYPEDSIGNATURES)

Check that the semi-axes are positive and finite \\[m\\], the initial angle finite and
``0<|span|≤2π`` \\[rad\\]. The normal offset is finite and above ``-ρ_{min}`` \\[m\\]. The
radius ``ρ_{min}`` is the least of ``ρ(θ)=(a^2\\sin^2θ+b^2\\cos^2θ)^{3/2}/(ab)`` over the
range of the arc. Return `arc`, or throw `DomainError`.
"""
function validate(arc::EllipticArc)
    isfinite(nominal(arc.a)) && nominal(arc.a) > 0 ||
        throw(DomainError(arc.a, "EllipticArc.a must be positive and finite"))
    isfinite(nominal(arc.b)) && nominal(arc.b) > 0 ||
        throw(DomainError(arc.b, "EllipticArc.b must be positive and finite"))
    isfinite(nominal(arc.φ0)) ||
        throw(DomainError(arc.φ0, "EllipticArc.φ0 must be finite"))
    isfinite(nominal(arc.span)) && 0 < abs(nominal(arc.span)) <= 2π ||
        throw(DomainError(arc.span, "EllipticArc.span must be finite with 0 < |span| ≤ 2π"))
    θ₁, θ₂ = minmax(arc.φ0, arc.φ0 + arc.span)
    ρ_min = min(hypot(arc.a * sin(θ₁), arc.b * cos(θ₁))^3,
        hypot(arc.a * sin(θ₂), arc.b * cos(θ₂))^3) / (arc.a * arc.b)
    if nominal(arc.a) > nominal(arc.b)
        if ceil(nominal(θ₁ / π)) <= floor(nominal(θ₂ / π))
            ρ_min = arc.b^2 / arc.a
        end
    elseif nominal(arc.a) < nominal(arc.b)
        if ceil(nominal((θ₁ - π / 2) / π)) <= floor(nominal((θ₂ - π / 2) / π))
            ρ_min = arc.a^2 / arc.b
        end
    else
        ρ_min = arc.a
    end
    isfinite(nominal(arc.offset)) && nominal(arc.offset) > -nominal(ρ_min) ||
        throw(DomainError(arc.offset,
            "EllipticArc.offset must be finite and above minus its smallest radius of curvature, $(ρ_min) m"))
    return arc
end

"""
$(TYPEDSIGNATURES)

Return ``at((a\\cos θ,b\\sin θ)+d(b\\cos θ,a\\sin θ)/ρ)`` \\[m\\], with ``θ=φ0+s\\,span`` \\[rad\\]
for ``s∈[0,1]`` \\[dimensionless\\], ``ρ=hypot(b cos θ,a sin θ)`` \\[m\\]
and ``d=offset`` \\[m\\].
"""
function (arc::EllipticArc)(s::Real)
    θ = arc.φ0 + s * arc.span
    ρ = hypot(arc.b * cos(θ), arc.a * sin(θ))
    return arc.at((arc.a * cos(θ) + arc.offset * arc.b * cos(θ) / ρ,
        arc.b * sin(θ) + arc.offset * arc.a * sin(θ) / ρ))
end

"""
$(TYPEDSIGNATURES)

Return the tangent direction of `curve` at `t` in \\[rad\\].

``θ=atan(b_2-a_2,b_1-a_1)``.
"""
function direction(curve::Segment, t::Real)
    return atan(curve.b[2] - curve.a[2], curve.b[1] - curve.a[1])
end

"""
$(TYPEDSIGNATURES)

Return the tangent direction of `curve` at `t` in \\[rad\\].

``θ=at.φ+t\\,span+sign(span)π/2``.
"""
function direction(curve::Arc, t::Real)
    return curve.at.φ + t * curve.span + sign(curve.span) * π / 2
end

"""
$(TYPEDSIGNATURES)

Return the tangent direction of `curve` at `t` in \\[rad\\].

``α=at.φ+atan(sign(span)b cos θ,-sign(span)a sin θ)``, with ``θ=φ0+t\\,span``.
"""
function direction(curve::EllipticArc, t::Real)
    θ = curve.φ0 + t * curve.span
    return curve.at.φ + atan(sign(curve.span) * curve.b * cos(θ),
        -sign(curve.span) * curve.a * sin(θ))
end

"""
$(TYPEDSIGNATURES)

Return the fraction of a revolution \\[dimensionless\\] and length \\[m\\] of the straight
curve `segment`. The fraction is zero, and the length is the distance between its endpoints
\\[m\\].
"""
function turn(segment::Segment{T}) where {T}
    return (fraction = zero(T), length = hypot((segment.b .- segment.a)...))
end

"""
$(TYPEDSIGNATURES)

Return the fraction of a revolution \\[dimensionless\\] and length \\[m\\] of the circular
curve `arc`. Its radius is in \\[m\\] and its signed span in \\[rad\\]. The fraction is
`abs(span) / (2π)` and the length is `abs(span) * r`.
"""
function turn(arc::Arc)
    span = abs(arc.span)
    return (fraction = span / (2π), length = span * arc.r)
end

"""
$(TYPEDSIGNATURES)

Return the normal's fraction of a revolution \\[dimensionless\\] and length \\[m\\] of the
elliptic curve `arc`, ``|Δα|/(2π)`` and ``|L+d Δα|``. Its semi-axes ``a,b`` and normal
offset ``d`` are in \\[m\\], and its parametric ends ``θ₁,θ₂`` and signed span are in
\\[rad\\]. The normal turns by ``Δα=span+δ(θ₂)-δ(θ₁)``, with
``δ(θ)=atan((a-b)sin(θ)cos(θ),b cos²(θ)+a sin²(θ))``. The base arc's signed length is
``L=∫_{θ₁}^{θ₂}hypot(a sin(θ),b cos(θ))dθ``.
"""
function turn(arc::EllipticArc{T}) where {T}
    a, b = arc.a, arc.b
    θ₁, θ₂ = arc.φ0, arc.φ0 + arc.span
    Δα = arc.span + atan((a - b) * sin(θ₂) * cos(θ₂),
        b * cos(θ₂)^2 + a * sin(θ₂)^2) -
         atan((a - b) * sin(θ₁) * cos(θ₁),
        b * cos(θ₁)^2 + a * sin(θ₁)^2)
    L = first(quadgk(θ -> hypot(a * sin(θ), b * cos(θ)), θ₁, θ₂;
        rtol = eps(float(nominal(one(T))))^(3 / 4),
        norm = x -> abs(nominal(x))))
    return (fraction = abs(Δα) / (2π), length = abs(L + arc.offset * Δα))
end

"""
$(TYPEDSIGNATURES)

Return the parallel of `curve` at outward distance `d` in \\[m\\].

``p_d=p+d(cos φ,sin φ)``, with ``φ=atan(-dx,dy)``.
"""
function parallel(curve::Segment, d::Real)
    dx, dy = curve.b .- curve.a
    φ = atan(-dx, dy)
    return Segment(curve.a .+ d .* (cos(φ), sin(φ)),
        curve.b .+ d .* (cos(φ), sin(φ)))
end

"""
$(TYPEDSIGNATURES)

Return the parallel of `curve` at outward distance `d` in \\[m\\].

``r_d=r+sign(span)d``.
"""
function parallel(curve::Arc, d::Real)
    return Arc(curve.at, curve.r + sign(curve.span) * d, curve.span)
end

"""
$(TYPEDSIGNATURES)

Return the parallel of `curve` at outward distance `d` in \\[m\\].

``offset_d=offset+d``.
"""
function parallel(curve::EllipticArc, d::Real)
    return EllipticArc(curve.at, curve.a, curve.b, curve.φ0, curve.span,
        curve.offset + d)
end

"""
$(TYPEDSIGNATURES)

Return the piece of the parallel of `curve` between `p` and `q` in \\[m\\].

``p_d(t)=p+t(q-p)``.
"""
function parallel(curve::Segment, d::Real, p, q)
    return Segment(p, q)
end

"""
$(TYPEDSIGNATURES)

Return the piece of the parallel of `curve` between `p` and `q` in \\[m\\].

``r_d=r+sign(span)d``, with span ``sign(span)2π\\,turn(curve,d,p,q).fraction``.
"""
function parallel(curve::Arc, d::Real, p, q)
    return Arc(
        Placement(curve.at.x, curve.at.y,
            atan(p[2] - curve.at.y, p[1] - curve.at.x)),
        curve.r + sign(curve.span) * d,
        sign(curve.span) * 2π * turn(curve, d, p, q).fraction)
end

"""
$(TYPEDSIGNATURES)

Return the piece of the parallel of `curve` between `p` and `q` in \\[m\\].

``offset_d=offset+d``.
"""
function parallel(curve::EllipticArc, d::Real, p, q)
    return parallel(curve, d)
end

"""
$(TYPEDSIGNATURES)

Return the end of the parallel of `a` and the start of the parallel of `b` at the distance
`d` \\[m\\]. The curve `a` ends at ``j = a(1)``, where `b` starts. At a tangent junction,
both points lie at the distance `d` from ``j`` along the outward normal. Where the two
parallels converge, both points are their intersection nearest ``j``. Where they diverge,
each point is the end of its own parallel.
"""
function parallel(a::Segment, b::Segment, d::Real)
    j = a(1)
    θa, θb = direction(a, 1), direction(b, 0)
    τ = π - mod(π - (θb - θa), 2π)
    u, v = (cos(θa), sin(θa)), (cos(θb), sin(θb))
    sine = nominal(u[1] * v[2] - u[2] * v[1])
    if nominal(u[1] * v[1] + u[2] * v[2]) > 0 &&
       abs(sine) <= sqrt(eps(float(one(sine))))
        p = j .+ d .* (cos(θa - π / 2), sin(θa - π / 2))
        return (p, p)
    elseif nominal(sin(τ) * d) < 0
        A, B = parallel(a, d), parallel(b, d)
        U, V = A.b .- A.a, B.b .- B.a
        w = B.a .- A.a
        t = (w[1] * V[2] - w[2] * V[1]) / (U[1] * V[2] - U[2] * V[1])
        p = A.a .+ t .* U
        return (p, p)
    else
        return (parallel(a, d)(1), parallel(b, d)(0))
    end
end

"""
$(TYPEDSIGNATURES)

Return the end of the parallel of `a` and the start of the parallel of `b` at the distance
`d` \\[m\\]. The curve `a` ends at ``j = a(1)``, where `b` starts. At a tangent junction,
both points lie at the distance `d` from ``j`` along the outward normal. Where the two
parallels converge, both points are their intersection nearest ``j``. Where they diverge,
each point is the end of its own parallel.
"""
function parallel(a::Segment, b::Arc, d::Real)
    j = a(1)
    θa, θb = direction(a, 1), direction(b, 0)
    τ = π - mod(π - (θb - θa), 2π)
    u, v = (cos(θa), sin(θa)), (cos(θb), sin(θb))
    sine = nominal(u[1] * v[2] - u[2] * v[1])
    if nominal(u[1] * v[1] + u[2] * v[2]) > 0 &&
       abs(sine) <= sqrt(eps(float(one(sine))))
        p = j .+ d .* (cos(θa - π / 2), sin(θa - π / 2))
        return (p, p)
    elseif nominal(sin(τ) * d) < 0
        A, B = parallel(a, d), parallel(b, d)
        U = A.b .- A.a
        w = A.a .- (B.at.x, B.at.y)
        L = U[1]^2 + U[2]^2
        h = -(w[1] * U[1] + w[2] * U[2]) / L
        Δ = h^2 - (w[1]^2 + w[2]^2 - B.r^2) / L
        δ = sqrt(max(Δ, zero(Δ)))
        points = (A.a .+ (h - δ) .* U, A.a .+ (h + δ) .* U)
        p = argmin(xy -> nominal(hypot((xy .- j)...)), points)
        return (p, p)
    else
        return (parallel(a, d)(1), parallel(b, d)(0))
    end
end

"""
$(TYPEDSIGNATURES)

Return the end of the parallel of `a` and the start of the parallel of `b` at the distance
`d` \\[m\\]. The curve `a` ends at ``j = a(1)``, where `b` starts. At a tangent junction,
both points lie at the distance `d` from ``j`` along the outward normal. Where the two
parallels converge, both points are their intersection nearest ``j``. Where they diverge,
each point is the end of its own parallel.
"""
function parallel(a::Arc, b::Segment, d::Real)
    j = a(1)
    θa, θb = direction(a, 1), direction(b, 0)
    τ = π - mod(π - (θb - θa), 2π)
    u, v = (cos(θa), sin(θa)), (cos(θb), sin(θb))
    sine = nominal(u[1] * v[2] - u[2] * v[1])
    if nominal(u[1] * v[1] + u[2] * v[2]) > 0 &&
       abs(sine) <= sqrt(eps(float(one(sine))))
        p = j .+ d .* (cos(θa - π / 2), sin(θa - π / 2))
        return (p, p)
    elseif nominal(sin(τ) * d) < 0
        A, B = parallel(b, d), parallel(a, d)
        U = A.b .- A.a
        w = A.a .- (B.at.x, B.at.y)
        L = U[1]^2 + U[2]^2
        h = -(w[1] * U[1] + w[2] * U[2]) / L
        Δ = h^2 - (w[1]^2 + w[2]^2 - B.r^2) / L
        δ = sqrt(max(Δ, zero(Δ)))
        points = (A.a .+ (h - δ) .* U, A.a .+ (h + δ) .* U)
        p = argmin(xy -> nominal(hypot((xy .- j)...)), points)
        return (p, p)
    else
        return (parallel(a, d)(1), parallel(b, d)(0))
    end
end

"""
$(TYPEDSIGNATURES)

Return the end of the parallel of `a` and the start of the parallel of `b` at the distance
`d` \\[m\\]. The curve `a` ends at ``j = a(1)``, where `b` starts. At a tangent junction,
both points lie at the distance `d` from ``j`` along the outward normal. Where the two
parallels converge, both points are their intersection nearest ``j``. Where they diverge,
each point is the end of its own parallel.
"""
function parallel(a::Arc, b::Arc, d::Real)
    j = a(1)
    θa, θb = direction(a, 1), direction(b, 0)
    τ = π - mod(π - (θb - θa), 2π)
    u, v = (cos(θa), sin(θa)), (cos(θb), sin(θb))
    sine = nominal(u[1] * v[2] - u[2] * v[1])
    if nominal(u[1] * v[1] + u[2] * v[2]) > 0 &&
       abs(sine) <= sqrt(eps(float(one(sine))))
        p = j .+ d .* (cos(θa - π / 2), sin(θa - π / 2))
        return (p, p)
    elseif nominal(sin(τ) * d) < 0
        A, B = parallel(a, d), parallel(b, d)
        w = (B.at.x - A.at.x, B.at.y - A.at.y)
        L = hypot(w...)
        U = w ./ L
        h = (A.r^2 - B.r^2 + L^2) / (2L)
        Δ = A.r^2 - h^2
        δ = sqrt(max(Δ, zero(Δ)))
        center = (A.at.x, A.at.y) .+ h .* U
        points = (center .+ δ .* (-U[2], U[1]), center .- δ .* (-U[2], U[1]))
        p = argmin(xy -> nominal(hypot((xy .- j)...)), points)
        return (p, p)
    else
        return (parallel(a, d)(1), parallel(b, d)(0))
    end
end

"""
$(TYPEDSIGNATURES)

Return the end of the parallel of `a` and the start of the parallel of `b` at the distance
`d` \\[m\\]. The curve `a` ends at ``j = a(1)``, where `b` starts. At a tangent junction,
both points lie at the distance `d` from ``j`` along the outward normal. Where the two
parallels diverge, each point is the end of its own parallel. Where they converge, the
method throws `DomainError`: the parallels of an elliptic arc and its neighbour meet at no
point of a conic.
"""
function parallel(a::EllipticArc, b::Segment, d::Real)
    j = a(1)
    θa, θb = direction(a, 1), direction(b, 0)
    τ = π - mod(π - (θb - θa), 2π)
    u, v = (cos(θa), sin(θa)), (cos(θb), sin(θb))
    sine = nominal(u[1] * v[2] - u[2] * v[1])
    if nominal(u[1] * v[1] + u[2] * v[2]) > 0 &&
       abs(sine) <= sqrt(eps(float(one(sine))))
        p = j .+ d .* (cos(θa - π / 2), sin(θa - π / 2))
        return (p, p)
    elseif nominal(sin(τ) * d) < 0
        throw(DomainError(j,
            "the parallels of an elliptic arc and its neighbour meet at no point of a conic at $(j) m"))
    else
        return (parallel(a, d)(1), parallel(b, d)(0))
    end
end

"""
$(TYPEDSIGNATURES)

Return the end of the parallel of `a` and the start of the parallel of `b` at the distance
`d` \\[m\\]. The curve `a` ends at ``j = a(1)``, where `b` starts. At a tangent junction,
both points lie at the distance `d` from ``j`` along the outward normal. Where the two
parallels diverge, each point is the end of its own parallel. Where they converge, the
method throws `DomainError`: the parallels of an elliptic arc and its neighbour meet at no
point of a conic.
"""
function parallel(a::EllipticArc, b::Arc, d::Real)
    j = a(1)
    θa, θb = direction(a, 1), direction(b, 0)
    τ = π - mod(π - (θb - θa), 2π)
    u, v = (cos(θa), sin(θa)), (cos(θb), sin(θb))
    sine = nominal(u[1] * v[2] - u[2] * v[1])
    if nominal(u[1] * v[1] + u[2] * v[2]) > 0 &&
       abs(sine) <= sqrt(eps(float(one(sine))))
        p = j .+ d .* (cos(θa - π / 2), sin(θa - π / 2))
        return (p, p)
    elseif nominal(sin(τ) * d) < 0
        throw(DomainError(j,
            "the parallels of an elliptic arc and its neighbour meet at no point of a conic at $(j) m"))
    else
        return (parallel(a, d)(1), parallel(b, d)(0))
    end
end

"""
$(TYPEDSIGNATURES)

Return the end of the parallel of `a` and the start of the parallel of `b` at the distance
`d` \\[m\\]. The curve `a` ends at ``j = a(1)``, where `b` starts. At a tangent junction,
both points lie at the distance `d` from ``j`` along the outward normal. Where the two
parallels diverge, each point is the end of its own parallel. Where they converge, the
method throws `DomainError`: the parallels of an elliptic arc and its neighbour meet at no
point of a conic.
"""
function parallel(a::EllipticArc, b::EllipticArc, d::Real)
    j = a(1)
    θa, θb = direction(a, 1), direction(b, 0)
    τ = π - mod(π - (θb - θa), 2π)
    u, v = (cos(θa), sin(θa)), (cos(θb), sin(θb))
    sine = nominal(u[1] * v[2] - u[2] * v[1])
    if nominal(u[1] * v[1] + u[2] * v[2]) > 0 &&
       abs(sine) <= sqrt(eps(float(one(sine))))
        p = j .+ d .* (cos(θa - π / 2), sin(θa - π / 2))
        return (p, p)
    elseif nominal(sin(τ) * d) < 0
        throw(DomainError(j,
            "the parallels of an elliptic arc and its neighbour meet at no point of a conic at $(j) m"))
    else
        return (parallel(a, d)(1), parallel(b, d)(0))
    end
end

"""
$(TYPEDSIGNATURES)

Return the end of the parallel of `a` and the start of the parallel of `b` at the distance
`d` \\[m\\]. The curve `a` ends at ``j = a(1)``, where `b` starts. At a tangent junction,
both points lie at the distance `d` from ``j`` along the outward normal. Where the two
parallels diverge, each point is the end of its own parallel. Where they converge, the
method throws `DomainError`: the parallels of an elliptic arc and its neighbour meet at no
point of a conic.
"""
function parallel(a::Segment, b::EllipticArc, d::Real)
    j = a(1)
    θa, θb = direction(a, 1), direction(b, 0)
    τ = π - mod(π - (θb - θa), 2π)
    u, v = (cos(θa), sin(θa)), (cos(θb), sin(θb))
    sine = nominal(u[1] * v[2] - u[2] * v[1])
    if nominal(u[1] * v[1] + u[2] * v[2]) > 0 &&
       abs(sine) <= sqrt(eps(float(one(sine))))
        p = j .+ d .* (cos(θa - π / 2), sin(θa - π / 2))
        return (p, p)
    elseif nominal(sin(τ) * d) < 0
        throw(DomainError(j,
            "the parallels of an elliptic arc and its neighbour meet at no point of a conic at $(j) m"))
    else
        return (parallel(a, d)(1), parallel(b, d)(0))
    end
end

"""
$(TYPEDSIGNATURES)

Return the end of the parallel of `a` and the start of the parallel of `b` at the distance
`d` \\[m\\]. The curve `a` ends at ``j = a(1)``, where `b` starts. At a tangent junction,
both points lie at the distance `d` from ``j`` along the outward normal. Where the two
parallels diverge, each point is the end of its own parallel. Where they converge, the
method throws `DomainError`: the parallels of an elliptic arc and its neighbour meet at no
point of a conic.
"""
function parallel(a::Arc, b::EllipticArc, d::Real)
    j = a(1)
    θa, θb = direction(a, 1), direction(b, 0)
    τ = π - mod(π - (θb - θa), 2π)
    u, v = (cos(θa), sin(θa)), (cos(θb), sin(θb))
    sine = nominal(u[1] * v[2] - u[2] * v[1])
    if nominal(u[1] * v[1] + u[2] * v[2]) > 0 &&
       abs(sine) <= sqrt(eps(float(one(sine))))
        p = j .+ d .* (cos(θa - π / 2), sin(θa - π / 2))
        return (p, p)
    elseif nominal(sin(τ) * d) < 0
        throw(DomainError(j,
            "the parallels of an elliptic arc and its neighbour meet at no point of a conic at $(j) m"))
    else
        return (parallel(a, d)(1), parallel(b, d)(0))
    end
end

"""
$(TYPEDSIGNATURES)

Return the turn of the parallel of `face` at the distance `d` \\[m\\] between the points `p`
and `q`. The fraction is zero, and the signed length \\[m\\] is ``(q-p)⋅u`` for the unit
direction ``u`` of `face`.
"""
function turn(face::Segment, d::Real, p, q)
    u = (face.b .- face.a) ./ hypot((face.b .- face.a)...)
    return (fraction = zero(face.a[1] + d),
        length = (q[1] - p[1]) * u[1] + (q[2] - p[2]) * u[2])
end

"""
$(TYPEDSIGNATURES)

Return the turn of the parallel of `face` at the distance `d` \\[m\\] between the points `p`
and `q`. The signed fraction of a revolution is ``θ/(2π)``, and the signed length is
``θ(r+\\mathrm{sign}(span)\\,d)`` \\[m\\]. The angle ``θ=|span|-m_0+m_1`` \\[rad\\] takes
the moves about the center from the start of `face` to `p` and from its end to `q`. Each
move counts in the sense in which `face` turns. Both values are negative once the face has
vanished.
"""
function turn(face::Arc, d::Real, p, q)
    P, Q = face(0), face(1)
    m₀ = sign(face.span) * (π - mod(
        π -
        (atan(p[2] - face.at.y, p[1] - face.at.x) -
         atan(P[2] - face.at.y, P[1] - face.at.x)), 2π))
    m₁ = sign(face.span) * (π - mod(
        π -
        (atan(q[2] - face.at.y, q[1] - face.at.x) -
         atan(Q[2] - face.at.y, Q[1] - face.at.x)), 2π))
    θ = abs(face.span) - m₀ + m₁
    return (fraction = θ / (2π), length = θ * (face.r + sign(face.span) * d))
end

"""
$(TYPEDSIGNATURES)

Return `(p, q, d, s, t)` for the nearest points `p` of `curve` and `q` of `other`. Their
distance is `d` \\[m\\], and their parameters give ``p = curve(s)`` and ``q = other(t)``.
Each parameter lies in ``[0, 1]`` and is exactly 0 or 1 at an end. Parallel segments take
the middle of their common range, at their perpendicular distance, or else their nearest
ends.
"""
function clearance(curve::Segment, other::Segment)
    u, v = curve.b .- curve.a, other.b .- other.a
    w = other.a .- curve.a
    z = zero(u[1] + v[1])
    candidates = [(s, t, hypot((curve(s) .- other(t))...))
                  for s in (z, one(z)) for t in (z, one(z))]
    cross = u[1] * v[2] - u[2] * v[1]
    if iszero(nominal(cross))
        L = hypot(u...)
        n = u ./ L
        α,
        β = w[1] * n[1] + w[2] * n[2],
        (other.b[1] - curve.a[1]) * n[1] + (other.b[2] - curve.a[2]) * n[2]
        lo = max(z, min(α, β))
        hi = min(L, max(α, β))
        if nominal(lo) <= nominal(hi)
            h = (lo + hi) / 2
            s = h / L
            t = (h - α) / (β - α)
            p, q = curve(s), other(t)
            return (p, q, abs(w[1] * n[2] - w[2] * n[1]), s, t)
        end
    else
        s = (w[1] * v[2] - w[2] * v[1]) / cross
        t = (w[1] * u[2] - w[2] * u[1]) / cross
        if 0 <= nominal(s) <= 1 && 0 <= nominal(t) <= 1
            push!(candidates, (s, t, z))
        end
    end
    for s in (z, one(z))
        p = curve(s)
        t = ((p[1] - other.a[1]) * v[1] + (p[2] - other.a[2]) * v[2]) /
            (v[1]^2 + v[2]^2)
        t = nominal(t) <= 0 ? z : nominal(t) >= 1 ? one(z) : t
        push!(candidates, (s, t, hypot((p .- other(t))...)))
    end
    for t in (z, one(z))
        q = other(t)
        s = ((q[1] - curve.a[1]) * u[1] + (q[2] - curve.a[2]) * u[2]) /
            (u[1]^2 + u[2]^2)
        s = nominal(s) <= 0 ? z : nominal(s) >= 1 ? one(z) : s
        push!(candidates, (s, t, hypot((curve(s) .- q)...)))
    end
    s, t, d = argmin(candidate -> nominal(candidate[3]), candidates)
    return (curve(s), other(t), d, s, t)
end

"""
$(TYPEDSIGNATURES)

Return `(p, q, d, s, t)` for the nearest points `p` of `curve` and `q` of `other`. Their
distance is `d` \\[m\\], and their parameters give ``p = curve(s)`` and ``q = other(t)``.
Each parameter lies in ``[0, 1]`` and is exactly 0 or 1 at an end. The pair is the nearest
among the ends, the feet of the normals and the crossings.
"""
function clearance(curve::Segment, other::Arc)
    u = curve.b .- curve.a
    z = zero(u[1] + other.r)
    candidates = [(s, t, hypot((curve(s) .- other(t))...))
                  for s in (z, one(z)) for t in (z, one(z))]
    L = u[1]^2 + u[2]^2
    for t in (z, one(z))
        q = other(t)
        s = ((q[1] - curve.a[1]) * u[1] + (q[2] - curve.a[2]) * u[2]) / L
        s = nominal(s) <= 0 ? z : nominal(s) >= 1 ? one(z) : s
        push!(candidates, (s, t, hypot((curve(s) .- q)...)))
    end
    for s in (z, one(z))
        p = curve(s)
        φ = atan(p[2] - other.at.y, p[1] - other.at.x)
        travel = mod(sign(other.span) * (φ - other.at.φ), 2π)
        if nominal(travel) <= abs(nominal(other.span))
            t = travel / abs(other.span)
            push!(candidates, (s, t, hypot((p .- other(t))...)))
        end
    end
    h = ((other.at.x - curve.a[1]) * u[1] +
         (other.at.y - curve.a[2]) * u[2]) / L
    if 0 <= nominal(h) <= 1
        s = h
        p = curve(s)
        for φ in (atan(-u[1], u[2]), atan(u[1], -u[2]))
            travel = mod(sign(other.span) * (φ - other.at.φ), 2π)
            if nominal(travel) <= abs(nominal(other.span))
                t = travel / abs(other.span)
                push!(candidates, (s, t, hypot((p .- other(t))...)))
            end
        end
    end
    w = curve.a .- (other.at.x, other.at.y)
    Δ = h^2 - (w[1]^2 + w[2]^2 - other.r^2) / L
    if nominal(Δ) >= 0
        for s in (h - sqrt(max(Δ, zero(Δ))), h + sqrt(max(Δ, zero(Δ))))
            if 0 <= nominal(s) <= 1
                p = curve(s)
                φ = atan(p[2] - other.at.y, p[1] - other.at.x)
                travel = mod(sign(other.span) * (φ - other.at.φ), 2π)
                if nominal(travel) <= abs(nominal(other.span))
                    t = travel / abs(other.span)
                    push!(candidates, (s, t, z))
                end
            end
        end
    end
    s, t, d = argmin(candidate -> nominal(candidate[3]), candidates)
    return (curve(s), other(t), d, s, t)
end

"""
$(TYPEDSIGNATURES)

Return `clearance(other, curve)` with its two points and its two parameters swapped,
`(q, p, d, t, s)`.
"""
function clearance(curve::Arc, other::Segment)
    q, p, d, t, s = clearance(other, curve)
    return (p, q, d, s, t)
end

"""
$(TYPEDSIGNATURES)

Return `(p, q, d, s, t)` for the nearest points `p` of `curve` and `q` of `other`. Their
distance is `d` \\[m\\], and their parameters give ``p = curve(s)`` and ``q = other(t)``.
Each parameter lies in ``[0, 1]`` and is exactly 0 or 1 at an end. Concentric arcs take the
middle of their common range, at ``|r_1-r_2|`` less the distance of their centers, or else
their nearest ends.
"""
function clearance(curve::Arc, other::Arc)
    z = zero(curve.r + other.r)
    candidates = [(s, t, hypot((curve(s) .- other(t))...))
                  for s in (z, one(z)) for t in (z, one(z))]
    w = (other.at.x - curve.at.x, other.at.y - curve.at.y)
    L = hypot(w...)
    if nominal(curve.at.x) == nominal(other.at.x) &&
       nominal(curve.at.y) == nominal(other.at.y)
        α = mod(curve.at.φ + min(curve.span, zero(curve.span)), 2π)
        β = mod(other.at.φ + min(other.span, zero(other.span)), 2π)
        ranges = [(max(α, β + k * 2π),
                      min(α + abs(curve.span), β + k * 2π + abs(other.span)), k)
                  for k in (-1, 0, 1)
                  if nominal(max(α, β + k * 2π)) <=
                     nominal(min(α + abs(curve.span), β + k * 2π + abs(other.span)))]
        if !isempty(ranges)
            lo, hi, k = argmax(interval -> nominal(interval[2] - interval[1]), ranges)
            φ = (lo + hi) / 2
            s = nominal(curve.span) > 0 ? (φ - α) / abs(curve.span) :
                (α + abs(curve.span) - φ) / abs(curve.span)
            t = nominal(other.span) > 0 ? (φ - β - k * 2π) / abs(other.span) :
                (β + k * 2π + abs(other.span) - φ) / abs(other.span)
            s = nominal(φ) == nominal(α) ?
                (nominal(curve.span) > 0 ? z : one(z)) :
                nominal(φ) == nominal(α + abs(curve.span)) ?
                (nominal(curve.span) > 0 ? one(z) : z) : s
            t = nominal(φ) == nominal(β + k * 2π) ?
                (nominal(other.span) > 0 ? z : one(z)) :
                nominal(φ) == nominal(β + k * 2π + abs(other.span)) ?
                (nominal(other.span) > 0 ? one(z) : z) : t
            return (curve(s), other(t), abs(curve.r - other.r) - L, s, t)
        end
    else
        for s in (z, one(z))
            p = curve(s)
            φ = atan(p[2] - other.at.y, p[1] - other.at.x)
            travel = mod(sign(other.span) * (φ - other.at.φ), 2π)
            if nominal(travel) <= abs(nominal(other.span))
                t = travel / abs(other.span)
                push!(candidates, (s, t, hypot((p .- other(t))...)))
            end
        end
        for t in (z, one(z))
            q = other(t)
            φ = atan(q[2] - curve.at.y, q[1] - curve.at.x)
            travel = mod(sign(curve.span) * (φ - curve.at.φ), 2π)
            if nominal(travel) <= abs(nominal(curve.span))
                s = travel / abs(curve.span)
                push!(candidates, (s, t, hypot((curve(s) .- q)...)))
            end
        end
        for φ in (atan(w[2], w[1]), atan(w[2], w[1]) + π),
            ψ in (atan(w[2], w[1]), atan(w[2], w[1]) + π)

            travel₁ = mod(sign(curve.span) * (φ - curve.at.φ), 2π)
            travel₂ = mod(sign(other.span) * (ψ - other.at.φ), 2π)
            if nominal(travel₁) <= abs(nominal(curve.span)) &&
               nominal(travel₂) <= abs(nominal(other.span))
                s, t = travel₁ / abs(curve.span), travel₂ / abs(other.span)
                push!(candidates, (s, t, hypot((curve(s) .- other(t))...)))
            end
        end
        h = (curve.r^2 - other.r^2 + L^2) / (2L)
        Δ = curve.r^2 - h^2
        if nominal(Δ) >= 0
            n = w ./ L
            center = (curve.at.x, curve.at.y) .+ h .* n
            for point in (center .+ sqrt(max(Δ, zero(Δ))) .* (-n[2], n[1]),
                center .- sqrt(max(Δ, zero(Δ))) .* (-n[2], n[1]))
                φ = atan(point[2] - curve.at.y, point[1] - curve.at.x)
                ψ = atan(point[2] - other.at.y, point[1] - other.at.x)
                travel₁ = mod(sign(curve.span) * (φ - curve.at.φ), 2π)
                travel₂ = mod(sign(other.span) * (ψ - other.at.φ), 2π)
                if nominal(travel₁) <= abs(nominal(curve.span)) &&
                   nominal(travel₂) <= abs(nominal(other.span))
                    s, t = travel₁ / abs(curve.span), travel₂ / abs(other.span)
                    push!(candidates, (s, t, z))
                end
            end
        end
    end
    s, t, d = argmin(candidate -> nominal(candidate[3]), candidates)
    return (curve(s), other(t), d, s, t)
end

"""
$(TYPEDSIGNATURES)

Return `(p, q, d, s, t)` for the nearest points `p` of `curve` and `q` of `other`. Their
distance is `d` \\[m\\], and their parameters give ``p = curve(s)`` and ``q = other(t)``.
Each parameter lies in ``[0, 1]`` and is exactly 0 or 1 at an end. Arcs of one ellipse take
the middle of their common range, at ``|o_1-o_2|`` for their offsets ``o_1`` and ``o_2``, or
else their nearest ends. Arcs of different ellipses throw `DomainError`.
"""
function clearance(curve::EllipticArc, other::EllipticArc)
    (nominal(curve.at.x), nominal(curve.at.y), nominal(curve.at.φ),
        nominal(curve.a), nominal(curve.b)) ==
    (nominal(other.at.x), nominal(other.at.y), nominal(other.at.φ),
        nominal(other.a), nominal(other.b)) ||
        throw(DomainError((curve, other), "the clearance of two elliptic arcs needs one ellipse"))
    z = zero(curve.a + other.a)
    candidates = [(s, t, hypot((curve(s) .- other(t))...))
                  for s in (z, one(z)) for t in (z, one(z))]
    α = mod(curve.φ0 + min(curve.span, zero(curve.span)), 2π)
    β = mod(other.φ0 + min(other.span, zero(other.span)), 2π)
    ranges = [(max(α, β + k * 2π),
                  min(α + abs(curve.span), β + k * 2π + abs(other.span)), k)
              for k in (-1, 0, 1)
              if nominal(max(α, β + k * 2π)) <=
                 nominal(min(α + abs(curve.span), β + k * 2π + abs(other.span)))]
    if !isempty(ranges)
        lo, hi, k = argmax(interval -> nominal(interval[2] - interval[1]), ranges)
        φ = (lo + hi) / 2
        s = nominal(curve.span) > 0 ? (φ - α) / abs(curve.span) :
            (α + abs(curve.span) - φ) / abs(curve.span)
        t = nominal(other.span) > 0 ? (φ - β - k * 2π) / abs(other.span) :
            (β + k * 2π + abs(other.span) - φ) / abs(other.span)
        s = nominal(φ) == nominal(α) ?
            (nominal(curve.span) > 0 ? z : one(z)) :
            nominal(φ) == nominal(α + abs(curve.span)) ?
            (nominal(curve.span) > 0 ? one(z) : z) : s
        t = nominal(φ) == nominal(β + k * 2π) ?
            (nominal(other.span) > 0 ? z : one(z)) :
            nominal(φ) == nominal(β + k * 2π + abs(other.span)) ?
            (nominal(other.span) > 0 ? one(z) : z) : t
        return (curve(s), other(t), abs(curve.offset - other.offset), s, t)
    end
    s, t, d = argmin(candidate -> nominal(candidate[3]), candidates)
    return (curve(s), other(t), d, s, t)
end

"""
A curve of a contour is a [`Segment`](@ref) an [`Arc`](@ref) or an [`EllipticArc`](@ref) in the floating type `T`.
"""
const Curve{T} = Union{Segment{T}, Arc{T}, EllipticArc{T}}

"""
$(TYPEDEF)

Represent the closed paths that bound a region of the plane. Each path runs
counterclockwise about the region it bounds. Each hole follows the same orientation.

$(TYPEDFIELDS)
"""
struct Contour{T <: Real}
    "Outer path of curves with coordinates in \\[m\\] and angles in \\[rad\\]."
    outer::Vector{Curve{T}}
    "Paths of holes with coordinates in \\[m\\] and angles in \\[rad\\]."
    holes::Vector{Vector{Curve{T}}}
    @doc """
    $(TYPEDSIGNATURES)

    Construct a contour from its `outer` path and the paths of its `holes`, whose
    coordinates are in \\[m\\] and angles in \\[rad\\]. Each path is nonempty and closed.
    The tolerance is `sqrt(eps(T))` times the largest absolute coordinate at a path's ends.
    That scale is at least one meter.

    # Errors

    - Throws `DomainError` if a path is empty or its ends do not meet.
    """
    function Contour(outer::Vector{Curve{T}},
            holes::Vector{Vector{Curve{T}}}) where {T <: Real}
        return validate(new{T}(outer, holes))
    end
end

function validate(contour::Contour{T}) where {T}
    for path in (contour.outer, contour.holes...)
        isempty(path) && throw(DomainError(path, "Contour paths must be nonempty"))
        scale = max(
            one(T), maximum(abs(coordinate)
            for curve in path for t in (0, 1) for coordinate in curve(t)))
        tolerance = sqrt(eps(T)) * scale
        for k in eachindex(path)
            endpoint = path[k](1)
            nextpoint = path[mod1(k + 1, length(path))](0)
            all(abs.(endpoint .- nextpoint) .<= tolerance) ||
                throw(DomainError(path, "Contour path must be closed at curve $k"))
        end
    end
    return contour
end

"""
$(TYPEDSIGNATURES)

Return the parallel of `contour` at the distance `d` \\[m\\], signed along the outward
normal of the region: positive outward, negative into it. Each curve of the outer path takes
`parallel(curve, d)`. Each curve of a hole's path takes `parallel(curve, -d)`, since that
path runs counterclockwise about the hole.

Each curve ends at its junction with the next, given by `parallel(a, b, d)`. Where the two
parallels meet, both end at that point. Where they diverge, an arc of radius ``|d|`` about
the corner joins the two ends. Outward on a convex outline, the result is the offset
outline.
"""
function parallel(contour::Contour, d::Real)
    T = promote_type(typeof(first(contour.outer)(0)[1]), typeof(d))
    paths = Vector{Curve{T}}[]
    for (i, path) in enumerate((contour.outer, contour.holes...))
        d′ = i == 1 ? d : -d
        if length(path) == 1
            push!(paths, Curve{T}[parallel(only(path), d′)])
            continue
        end
        junctions = [parallel(path[k], path[mod1(k + 1, length(path))], d′)
                     for k in eachindex(path)]
        outline = Curve{T}[]
        for k in eachindex(path)
            p, q = junctions[k]
            push!(outline, parallel(path[k], d′,
                junctions[mod1(k - 1, length(path))][2], p))
            if nominal.(p) != nominal.(q)
                j = path[k](1)
                θa = direction(path[k], 1)
                θb = direction(path[mod1(k + 1, length(path))], 0)
                τ = π - mod(π - (θb - θa), 2π)
                push!(outline, Arc(Placement(j..., atan(p[2] - j[2], p[1] - j[1])),
                    abs(d′), τ))
            end
        end
        push!(paths, outline)
    end
    return Contour(first(paths), paths[2:end])
end

"""
$(TYPEDSIGNATURES)

Return the distance across the metal of `shape` \\[m\\]: twice the depth at which the inward
parallel of its outline first degenerates. The distance is the least of three bounds:

```math
d = \\min\\left(\\min d_{ab},\\; 2 \\min \\rho,\\; 2 \\min e_c\\right).
```

The width ``d_{ab}`` is the `clearance` of two curves that are not adjacent, or that lie on
different paths. It counts where their inward parallels meet. There the chord between its
two points lies in the inward normal cone of the outline at both ends. Its midpoint lies in
the metal. The radius ``\\rho`` runs over the arcs that turn with the metal. For an elliptic
arc, it is the smallest radius of curvature of its parallel.

The depth ``e_c`` is where the `parallel` piece of a face ``c`` with a corner at an end
vanishes between the parallels of its neighbours. Bisection of the signed length of the
piece finds it. A segment of length ``L`` between two segments at the interior angles
``\\theta_1`` and ``\\theta_2`` \\[rad\\] has
``e_c = L / (\\cot(\\theta_1/2) + \\cot(\\theta_2/2))``.
"""
function dimension(shape::AbstractShape)
    contour = Contour(shape)
    T = eltype(shape)
    bounds = T[]
    paths = (contour.outer, contour.holes...)
    for (i, path) in enumerate(paths), (j, otherpath) in enumerate(paths)

        j < i && continue
        for (k, a) in enumerate(path), (l, b) in enumerate(otherpath)

            if i == j && (l <= k || l == k + 1 || (k == 1 && l == length(path)))
                continue
            end
            p, q, w, s, t = clearance(a, b)
            u = (q .- p) ./ w
            counts = true
            for (curve, parameter, chord, side, index, edges) in
                ((a, s, u, i == 1 ? 1 : -1, k, path),
                (b, t, .-u, j == 1 ? 1 : -1, l, otherpath))
                if 0 < nominal(parameter) < 1
                    θ = direction(curve, parameter) + side * π / 2
                    n = (cos(θ), sin(θ))
                    sine = nominal(chord[1] * n[2] - chord[2] * n[1])
                    counts &= nominal(chord[1] * n[1] + chord[2] * n[2]) > 0 &&
                              abs(sine) <= sqrt(eps(float(one(sine))))
                else
                    before = nominal(parameter) == 0 ?
                             edges[mod1(index - 1, length(edges))] : curve
                    after = nominal(parameter) == 0 ?
                            curve : edges[mod1(index + 1, length(edges))]
                    θa, θb = direction(before, 1), direction(after, 0)
                    τ = π - mod(π - (θb - θa), 2π)
                    na = (cos(θa + side * π / 2), sin(θa + side * π / 2))
                    nb = (cos(θb + side * π / 2), sin(θb + side * π / 2))
                    sine = nominal(na[1] * nb[2] - na[2] * nb[1])
                    if nominal(na[1] * nb[1] + na[2] * nb[2]) > 0 &&
                       abs(sine) <= sqrt(eps(float(one(sine))))
                        cross = nominal(chord[1] * na[2] - chord[2] * na[1])
                        counts &= nominal(chord[1] * na[1] + chord[2] * na[2]) > 0 &&
                                  abs(cross) <= sqrt(eps(float(one(cross))))
                    elseif nominal(side * τ) < 0
                        cross₁ = nominal(na[1] * chord[2] - na[2] * chord[1])
                        cross₂ = nominal(chord[1] * nb[2] - chord[2] * nb[1])
                        counts &= side * cross₁ <= sqrt(eps(float(one(cross₁)))) &&
                                  side * cross₂ <= sqrt(eps(float(one(cross₂))))
                    else
                        counts = false
                    end
                end
            end
            if counts
                m = (p .+ q) ./ 2
                if iszero(nominal(distance(Contour(resolve(Placement(-m[1], -m[2]), shape)))))
                    push!(bounds, w)
                end
            end
        end
    end
    for (i, path) in enumerate(paths)
        side = i == 1 ? 1 : -1
        for curve in [curve for curve in path if curve isa Arc]
            if nominal(side * sign(curve.span)) > 0
                push!(bounds, 2 * curve.r)
            end
        end
        for curve in [curve for curve in path if curve isa EllipticArc]
            if nominal(side * sign(curve.span)) > 0
                θ₁, θ₂ = minmax(curve.φ0, curve.φ0 + curve.span)
                ρ_min = min(hypot(curve.a * sin(θ₁), curve.b * cos(θ₁))^3,
                    hypot(curve.a * sin(θ₂), curve.b * cos(θ₂))^3) / (curve.a * curve.b)
                if nominal(curve.a) > nominal(curve.b)
                    if ceil(nominal(θ₁ / π)) <= floor(nominal(θ₂ / π))
                        ρ_min = curve.b^2 / curve.a
                    end
                elseif nominal(curve.a) < nominal(curve.b)
                    if ceil(nominal((θ₁ - π / 2) / π)) <= floor(nominal((θ₂ - π / 2) / π))
                        ρ_min = curve.a^2 / curve.b
                    end
                else
                    ρ_min = curve.a
                end

                push!(bounds, 2 * (ρ_min + curve.offset))
            end
        end
    end
    width = min(support(shape, 0) + support(shape, π),
        support(shape, π / 2) + support(shape, 3π / 2))
    limit = minimum([bounds; width]) / 2
    for (i, path) in enumerate(paths)
        side = i == 1 ? 1 : -1
        for k in eachindex(path)
            a, c,
            b = path[mod1(k - 1, length(path))], path[k],
            path[mod1(k + 1, length(path))]
            corners = false
            for (before, after) in ((a, c), (c, b))
                θa, θb = direction(before, 1), direction(after, 0)
                U, V = (cos(θa), sin(θa)), (cos(θb), sin(θb))
                sine = nominal(U[1] * V[2] - U[2] * V[1])
                corners |= !(nominal(U[1] * V[1] + U[2] * V[2]) > 0 &&
                             abs(sine) <= sqrt(eps(float(one(sine)))))
            end
            corners || continue
            lo, hi = zero(limit), limit
            while true
                mid = (lo + hi) / 2
                (nominal(mid) == nominal(lo) || nominal(mid) == nominal(hi)) && break
                if nominal(turn(c, -side * mid,
                    parallel(a, c, -side * mid)[2],
                    parallel(c, b, -side * mid)[1]).length) > 0
                    lo = mid
                else
                    hi = mid
                end
            end
            nominal(hi) == nominal(limit) || push!(bounds, 2 * hi)
        end
    end
    return minimum(bounds)
end

"""
$(TYPEDEF)

Represent a region bounded by one closed path of segments, circular arcs and elliptic arcs.

$(TYPEDFIELDS)
"""
struct ArcPolygon{T <: Real} <: AbstractShape{T}
    "Closed path of curves that bounds the region, counterclockwise, in the frame of the cable."
    outline::Vector{Curve{T}}
    @doc """
    $(TYPEDSIGNATURES)

    Construct a region from its counterclockwise `outline`, with coordinates in \\[m\\]
    and angles in \\[rad\\], through `validate`.
    """
    function ArcPolygon(outline::Vector{Curve{T}}) where {T <: Real}
        return validate(new{T}(outline))
    end
end

"""
$(TYPEDSIGNATURES)

Check that the outline of `shape` is one closed path, by the checks of `Contour`, and
return `shape`.
Coordinates are in \\[m\\] and angles in \\[rad\\].
"""
function validate(shape::ArcPolygon{T}) where {T}
    Contour(shape.outline, Vector{Curve{T}}[])
    return shape
end

"""
$(TYPEDSIGNATURES)

Return `shape.outline` as the outer path with no holes.
Coordinates are in \\[m\\] and angles in \\[rad\\].
"""
function Contour(shape::ArcPolygon{T}) where {T}
    return Contour(shape.outline, Vector{Curve{T}}[])
end

area(shape::ArcPolygon) = area(Contour(shape))
centroid(shape::ArcPolygon) = centroid(Contour(shape))
perimeter(shape::ArcPolygon) = perimeter(Contour(shape))
support(shape::ArcPolygon) = support(Contour(shape))
support(shape::ArcPolygon, φ::Real) = support(Contour(shape), φ)
distance(shape::ArcPolygon) = distance(Contour(shape))
boundary(shape::ArcPolygon) = shape

function Base.issubset(B::AbstractShape, A::Union{Rectangle, Polygon, Ellipse, ArcPolygon})
    T = promote_type(eltype(A), eltype(B))
    paths = (Contour(A).outer, Contour(B).outer)
    for edge in [curve for path in paths for curve in path if curve isa Segment]
        dx, dy = edge.b .- edge.a
        φ = atan(-dx, dy)
        radius = support(A, φ)
        tolerance = sqrt(eps(float(nominal(radius)))) *
                    max(one(nominal(radius)), nominal(radius))
        nominal(support(B, φ)) <= nominal(radius) + tolerance || return false
    end
    points = [[curve(t) for curve in path for t in (zero(T), one(T))] for path in paths]
    arcs = [[(at = curve.at, a = curve.r, b = curve.r,
                 φ0 = zero(T), span = curve.span, offset = zero(T))
             for curve in path if curve isa Arc] for path in paths]
    for (i, path) in enumerate(paths)
        append!(arcs[i],
            [(at = curve.at, a = curve.a, b = curve.b,
                 φ0 = curve.φ0, span = curve.span, offset = curve.offset)
             for curve in path if curve isa EllipticArc])
    end
    c = centroid(A)
    d = centroid(B) .- c
    direction = atan(d[2], d[1])
    for arc in Iterators.flatten(arcs)
        start = arc.at.φ + atan(arc.a * sin(arc.φ0), arc.b * cos(arc.φ0))
        finish = arc.at.φ + atan(arc.a * sin(arc.φ0 + arc.span),
            arc.b * cos(arc.φ0 + arc.span))
        span = sign(arc.span) * mod(sign(arc.span) * (finish - start), T(2π))
        abs(nominal(arc.span)) >= 2π && (span = sign(arc.span) * T(2π))
        nominal(arc.span) < 0 && (start += T(π))
        low = start + (nominal(span) < 0 ? span : zero(T))
        high = start + (nominal(span) > 0 ? span : zero(T))
        φ = direction + T(2π) * round(nominal(((low + high) / 2 - direction) / T(2π)))
        φ = nominal(φ) < nominal(low) ? low : (nominal(φ) > nominal(high) ? high : φ)
        for angle in (low, high)
            radius = support(A, angle)
            tolerance = sqrt(eps(float(nominal(radius)))) *
                        max(one(nominal(radius)), nominal(radius))
            nominal(support(B, angle)) <= nominal(radius) + tolerance || return false
        end
        for _ in 1:100
            values = [(zero(T), zero(T), zero(T)) for _ in paths]
            for i in eachindex(paths)
                point = first(points[i])
                value = point[1] * cos(φ) + point[2] * sin(φ)
                values[i] = (value, -point[1] * sin(φ) + point[2] * cos(φ), -value)
                for point in points[i]
                    value = point[1] * cos(φ) + point[2] * sin(φ)
                    if nominal(value) > nominal(values[i][1])
                        values[i] = (value, -point[1] * sin(φ) + point[2] * cos(φ), -value)
                    end
                end
                for curve in arcs[i]
                    ψ = φ - curve.at.φ
                    θ = atan(curve.b * sin(ψ), curve.a * cos(ψ))
                    nominal(mod(sign(curve.span) * (θ - curve.φ0), T(2π))) <=
                    abs(nominal(curve.span)) || continue
                    h = hypot(curve.a * cos(ψ), curve.b * sin(ψ))
                    projection = curve.at.x * cos(φ) + curve.at.y * sin(φ)
                    value = projection + h + curve.offset
                    if nominal(value) >= nominal(values[i][1])
                        derivative = (curve.b^2 - curve.a^2) * sin(ψ) * cos(ψ) / h
                        values[i] = (value,
                            -curve.at.x * sin(φ) + curve.at.y * cos(φ) + derivative,
                            -projection + (curve.b^2 - curve.a^2) * cos(2ψ) / h -
                            derivative^2 / h)
                    end
                end
            end
            derivative = values[1][2] - values[2][2]
            second = values[1][3] - values[2][3]
            iszero(nominal(derivative)) && break
            nominal(second) > 0 || break
            step = derivative / second
            angle = φ - step
            angle = nominal(angle) < nominal(low) ? low :
                    (nominal(angle) > nominal(high) ? high : angle)
            tolerance = sqrt(eps(float(nominal(φ)))) * max(one(nominal(φ)), abs(nominal(φ)))
            stop = abs(nominal(angle - φ)) <= tolerance
            φ = angle
            stop && break
        end
        radius = support(A, φ)
        tolerance = sqrt(eps(float(nominal(radius)))) *
                    max(one(nominal(radius)), nominal(radius))
        nominal(support(B, φ)) <= nominal(radius) + tolerance || return false
    end
    return true
end

"""
$(TYPEDSIGNATURES)

Return whether ``B⊆A.outer`` and `B` lies apart from every hole.
A circular hole of radius ``r`` needs ``distance(B-c)+tol≥r``.
Another hole needs ``support(hole,ψ)≤-support(B,ψ+π)+tol`` in the outward normal of a
segment of the hole, or in the reversed normal of a segment of `B`.
Lengths and tolerances are in \\[m\\], angles in \\[rad\\].
The test is exact for polygonal outlines. Only segments give directions, so the method
refuses a content beside an arc of a hole other than a circle.
"""
function Base.issubset(B::AbstractShape, A::PerforatedShape)
    B ⊆ A.outer || return false
    for hole in [hole for hole in A.holes if hole isa Circle]
        radius = hole.r
        tolerance = sqrt(eps(float(nominal(radius)))) *
                    max(one(nominal(radius)), nominal(radius))
        nominal(distance(resolve(Placement(-hole.at.x, -hole.at.y), B))) + tolerance >=
        nominal(radius) || return false
    end
    for hole in [hole for hole in A.holes if !(hole isa Circle)]
        directions = [atan(-(edge.b[1] - edge.a[1]), edge.b[2] - edge.a[2])
                      for edge in Contour(hole).outer if edge isa Segment]
        append!(directions,
            [atan(-(edge.b[1] - edge.a[1]), edge.b[2] - edge.a[2]) + π
             for edge in Contour(B).outer if edge isa Segment])
        apart = false
        for ψ in directions
            radius = support(B, ψ + π)
            tolerance = sqrt(eps(float(nominal(radius)))) *
                        max(one(nominal(radius)), nominal(radius))
            apart |= nominal(support(hole, ψ)) <= -nominal(radius) + tolerance
        end
        apart || return false
    end
    return true
end

"""
$(TYPEDSIGNATURES)

Place the outline by ``P\\mapsto\\mathrm{at}(P)`` for segment ends and
``c\\mapsto\\mathrm{at}\\,c`` for arc placements. Coordinates are in \\[m\\]
and rotations in \\[rad\\]. Keep each elliptic arc's normal offset \\[m\\].
"""
function resolve(at::Placement, shape::ArcPolygon)
    T = promote_type(eltype(at), eltype(shape))
    outline = Curve{T}[]
    for curve in [curve for curve in shape.outline if curve isa Segment]
        push!(outline, Segment(at(curve.a), at(curve.b)))
    end
    for (k, curve) in [(k, curve)
         for (k, curve) in enumerate(shape.outline)
         if curve isa Union{Arc, EllipticArc}]
        insert!(outline, k, curve)
    end
    for (k, curve) in
        [(k, curve) for (k, curve) in enumerate(shape.outline) if curve isa Arc]
        outline[k] = Arc(at * curve.at, curve.r, curve.span)
    end
    for (k, curve) in
        [(k, curve) for (k, curve) in enumerate(shape.outline) if curve isa EllipticArc]
        outline[k] = EllipticArc(
            at * curve.at, curve.a, curve.b, curve.φ0, curve.span, curve.offset)
    end
    return ArcPolygon(outline)
end

"""
$(TYPEDSIGNATURES)

Map the curves' coordinates \\[m\\], radii \\[m\\], angles \\[rad\\] and normal offsets \\[m\\] of `shape` by ``x\\mapsto f(x)``.
"""
function Base.map(f, shape::ArcPolygon)
    outline = similar(shape.outline)
    for (k, curve) in
        [(k, curve) for (k, curve) in enumerate(shape.outline) if curve isa Segment]
        outline[k] = Segment(map(f, curve.a), map(f, curve.b))
    end
    for (k, curve) in
        [(k, curve) for (k, curve) in enumerate(shape.outline) if curve isa Arc]
        outline[k] = Arc(map(f, curve.at), f(curve.r), f(curve.span))
    end
    for (k, curve) in
        [(k, curve) for (k, curve) in enumerate(shape.outline) if curve isa EllipticArc]
        outline[k] = EllipticArc(map(f, curve.at), f(curve.a), f(curve.b),
            f(curve.φ0), f(curve.span), f(curve.offset))
    end
    return ArcPolygon(outline)
end

function area(contour::Contour{T}) where {T}
    result = zero(T)
    for (k, path) in enumerate((contour.outer, contour.holes...))
        orientation = k == 1 ? 1 : -1
        for curve in [curve for curve in path if curve isa Segment]
            P, Q = curve.a, curve.b
            result += orientation * (P[1] * Q[2] - P[2] * Q[1]) / 2
        end
        for curve in [curve for curve in path if curve isa Arc]
            P, Q = curve(0), curve(1)
            result += orientation * (curve.r^2 * curve.span +
                       curve.at.x * (Q[2] - P[2]) - curve.at.y * (Q[1] - P[1])) / 2
        end
        for curve in [curve for curve in path if curve isa EllipticArc]
            θ₁, θ₂ = curve.φ0, curve.φ0 + curve.span
            P = curve.at((curve.a * cos(θ₁), curve.b * sin(θ₁)))
            Q = curve.at((curve.a * cos(θ₂), curve.b * sin(θ₂)))
            result += orientation * (curve.a * curve.b * curve.span +
                       curve.at.x * (Q[2] - P[2]) - curve.at.y * (Q[1] - P[1])) / 2
            a, b, d = curve.a, curve.b, curve.offset
            Δα = curve.span + atan((a - b) * sin(θ₂) * cos(θ₂),
                b * cos(θ₂)^2 + a * sin(θ₂)^2) -
                 atan((a - b) * sin(θ₁) * cos(θ₁),
                b * cos(θ₁)^2 + a * sin(θ₁)^2)
            L = first(quadgk(θ -> hypot(a * sin(θ), b * cos(θ)), θ₁, θ₂;
                rtol = eps(float(nominal(one(T))))^(3 / 4),
                norm = x -> abs(nominal(x))))
            c, s = cos(curve.at.φ), sin(curve.at.φ)
            τ₁ = ((-a * sin(θ₁) * c - b * cos(θ₁) * s),
                (-a * sin(θ₁) * s + b * cos(θ₁) * c)) ./ hypot(b * cos(θ₁), a * sin(θ₁))
            τ₂ = ((-a * sin(θ₂) * c - b * cos(θ₂) * s),
                (-a * sin(θ₂) * s + b * cos(θ₂) * c)) ./ hypot(b * cos(θ₂), a * sin(θ₂))
            result += orientation * (d * L -
                       d / 2 *
                       (Q[1] * τ₂[1] + Q[2] * τ₂[2] -
                        P[1] * τ₁[1] - P[2] * τ₁[2]) + d^2 / 2 * Δα)
        end
    end
    return result
end

function centroid(contour::Contour{T}) where {T}
    x, y = zero(T), zero(T)
    for (k, path) in enumerate((contour.outer, contour.holes...))
        orientation = k == 1 ? 1 : -1
        for curve in [curve for curve in path if curve isa Segment]
            P, Q = curve.a, curve.b
            cross = P[1] * Q[2] - P[2] * Q[1]
            x += orientation * (P[1] + Q[1]) * cross / 6
            y += orientation * (P[2] + Q[2]) * cross / 6
        end
        for curve in [curve for curve in path if curve isa Arc]
            P, Q = curve(0), curve(1)
            cross = P[1] * Q[2] - P[2] * Q[1]
            θ₁, θ₂ = curve.at.φ, curve.at.φ + curve.span
            ρ = curve.r
            A = ρ^2 * (curve.span - sin(curve.span)) / 2
            x += orientation * ((P[1] + Q[1]) * cross / 6 + curve.at.x * A +
                  ρ^3 * (sin(θ₂) - sin(θ₁)) / 3 -
                  ρ^2 * sin(curve.span) * (P[1] + Q[1] - 2curve.at.x) / 6)
            y += orientation * ((P[2] + Q[2]) * cross / 6 + curve.at.y * A +
                  ρ^3 * (cos(θ₁) - cos(θ₂)) / 3 -
                  ρ^2 * sin(curve.span) * (P[2] + Q[2] - 2curve.at.y) / 6)
        end
        for curve in [curve for curve in path if curve isa EllipticArc]
            θ₁, θ₂ = curve.φ0, curve.φ0 + curve.span
            P = curve.at((curve.a * cos(θ₁), curve.b * sin(θ₁)))
            Q = curve.at((curve.a * cos(θ₂), curve.b * sin(θ₂)))
            cross = P[1] * Q[2] - P[2] * Q[1]
            A = curve.a * curve.b * (curve.span - sin(curve.span)) / 2
            u = (sin(θ₂) - sin(θ₁)) / 3 -
                sin(curve.span) * (cos(θ₁) + cos(θ₂)) / 6
            v = (cos(θ₁) - cos(θ₂)) / 3 -
                sin(curve.span) * (sin(θ₁) + sin(θ₂)) / 6
            c, s = cos(curve.at.φ), sin(curve.at.φ)
            x += orientation * ((P[1] + Q[1]) * cross / 6 + curve.at.x * A +
                  curve.a * curve.b * (curve.a * c * u - curve.b * s * v))
            y += orientation * ((P[2] + Q[2]) * cross / 6 + curve.at.y * A +
                  curve.a * curve.b * (curve.a * s * u + curve.b * c * v))
            a, b, d = curve.a, curve.b, curve.offset
            Δα = curve.span + atan((a - b) * sin(θ₂) * cos(θ₂),
                b * cos(θ₂)^2 + a * sin(θ₂)^2) -
                 atan((a - b) * sin(θ₁) * cos(θ₁),
                b * cos(θ₁)^2 + a * sin(θ₁)^2)
            L = first(quadgk(θ -> hypot(a * sin(θ), b * cos(θ)), θ₁, θ₂;
                rtol = eps(float(nominal(one(T))))^(3 / 4),
                norm = x -> abs(nominal(x))))
            integrals = Tuple{T, T, T, T}[]
            for θ in (θ₁, θ₂)
                ξ, η = sin(θ), cos(θ)
                if nominal(a) > nominal(b)
                    e = sqrt(a^2 - b^2)
                    G = (ξ * hypot(b, e * ξ) + b^2 / e * asinh(e * ξ / b)) / 2
                    H = (η * sqrt(a^2 - e^2 * η^2) + a^2 / e * asin(e * η / a)) / 2
                    J = atan(e * ξ / b) / (b * e)
                    K = atanh(e * η / a) / (a * e)
                elseif nominal(a) < nominal(b)
                    e = sqrt(b^2 - a^2)
                    G = (ξ * sqrt(b^2 - e^2 * ξ^2) + b^2 / e * asin(e * ξ / b)) / 2
                    H = (η * hypot(a, e * η) + a^2 / e * asinh(e * η / a)) / 2
                    J = atanh(e * ξ / b) / (b * e)
                    K = atan(e * η / a) / (a * e)
                else
                    G, H, J, K = b * ξ, a * η, ξ / b^2, η / a^2
                end
                push!(integrals, (a * G, -b * H, a^2 * b * J, -a * b^2 * K))
            end
            I = last(integrals) .- first(integrals)
            x += orientation * (d * (curve.at.x * L + c * I[1] - s * I[2]) +
                  d^2 / 2 * (curve.at.x * Δα + c * I[3] - s * I[4]) +
                  d^2 / 2 * (Q[2] - P[2]))
            y += orientation * (d * (curve.at.y * L + s * I[1] + c * I[2]) +
                  d^2 / 2 * (curve.at.y * Δα + s * I[3] + c * I[4]) -
                  d^2 / 2 * (Q[1] - P[1]))
            n₁ = ((b * cos(θ₁) * c - a * sin(θ₁) * s),
                (b * cos(θ₁) * s + a * sin(θ₁) * c)) ./
                 hypot(b * cos(θ₁), a * sin(θ₁))
            n₂ = ((b * cos(θ₂) * c - a * sin(θ₂) * s),
                (b * cos(θ₂) * s + a * sin(θ₂) * c)) ./
                 hypot(b * cos(θ₂), a * sin(θ₂))
            x += orientation * d^3 / 3 * (n₂[2] - n₁[2])
            y -= orientation * d^3 / 3 * (n₂[1] - n₁[1])
            for (U, V) in ((Q, Q .+ d .* n₂), (P .+ d .* n₁, P))
                cross = U[1] * V[2] - U[2] * V[1]
                x += orientation * (U[1] + V[1]) * cross / 6
                y += orientation * (U[2] + V[2]) * cross / 6
            end
        end
    end
    return (x, y) ./ area(contour)
end

function perimeter(contour::Contour{T}) where {T}
    result = zero(T)
    for path in (contour.outer, contour.holes...)
        for curve in [curve for curve in path if curve isa Segment]
            result += hypot((curve.b .- curve.a)...)
        end
        for curve in [curve for curve in path if curve isa Arc]
            result += curve.r * abs(curve.span)
        end
        for curve in [curve for curve in path if curve isa EllipticArc]
            θ₁, θ₂ = curve.φ0, curve.φ0 + curve.span
            a, b = curve.a, curve.b
            Δα = curve.span + atan((a - b) * sin(θ₂) * cos(θ₂),
                b * cos(θ₂)^2 + a * sin(θ₂)^2) -
                 atan((a - b) * sin(θ₁) * cos(θ₁),
                b * cos(θ₁)^2 + a * sin(θ₁)^2)
            result += abs(first(quadgk(θ -> hypot(curve.a * sin(θ), curve.b * cos(θ)),
                curve.φ0, curve.φ0 + curve.span;
                rtol = eps(float(nominal(one(T))))^(3 / 4),
                norm = x -> abs(nominal(x)))) + curve.offset * Δα)
        end
    end
    return result
end

function support(contour::Contour{T}) where {T}
    result = zero(T)
    for curve in [curve for curve in contour.outer if curve isa Segment]
        result = max(result, hypot(curve.a...), hypot(curve.b...))
    end
    for curve in [curve for curve in contour.outer if curve isa Arc]
        φ = atan(curve.at.y, curve.at.x)
        result = max(result, hypot(curve(0)...), hypot(curve(1)...))
        if mod(sign(curve.span) * (φ - curve.at.φ), 2π) <= abs(curve.span)
            result = max(result, hypot(curve.at.x, curve.at.y) + curve.r)
        end
    end
    for curve in [curve for curve in contour.outer if curve isa EllipticArc]
        result = max(result, hypot(curve.at.x, curve.at.y) + max(curve.a, curve.b) +
                             curve.offset)
    end
    return result
end

function support(contour::Contour, φ::Real)
    point = first(contour.outer)(0)
    result = point[1] * cos(φ) + point[2] * sin(φ)
    for curve in [curve for curve in contour.outer if curve isa Segment]
        for point in (curve.a, curve.b)
            result = max(result, point[1] * cos(φ) + point[2] * sin(φ))
        end
    end
    for curve in [curve for curve in contour.outer if curve isa Arc]
        for point in (curve(0), curve(1))
            result = max(result, point[1] * cos(φ) + point[2] * sin(φ))
        end
        if mod(sign(curve.span) * (φ - curve.at.φ), 2π) <= abs(curve.span)
            result = max(result, curve.at.x * cos(φ) + curve.at.y * sin(φ) + curve.r)
        end
    end
    for curve in [curve for curve in contour.outer if curve isa EllipticArc]
        for point in (curve(0), curve(1))
            result = max(result, point[1] * cos(φ) + point[2] * sin(φ))
        end
        ψ = φ - curve.at.φ
        θ = atan(curve.b * sin(ψ), curve.a * cos(ψ))
        if nominal(mod(sign(curve.span) * (θ - curve.φ0), 2π)) <= abs(nominal(curve.span))
            result = max(result,
                curve.at.x * cos(φ) + curve.at.y * sin(φ) +
                hypot(curve.a * cos(ψ), curve.b * sin(ψ)) + curve.offset)
        end
    end
    return result
end

function distance(contour::Contour{T}) where {T}
    result = hypot(first(contour.outer)(0)...)
    inside = true
    for (k, path) in enumerate((contour.outer, contour.holes...))
        winding = 0
        for curve in [curve for curve in path if curve isa Segment]
            P, Q = curve.a, curve.b
            dx, dy = Q .- P
            t = clamp(-(P[1] * dx + P[2] * dy) / (dx^2 + dy^2), zero(T), one(T))
            result = min(result, hypot(curve(t)...))
            cross = P[1] * Q[2] - P[2] * Q[1]
            if P[2] <= zero(T) < Q[2] && cross > zero(T)
                winding += 1
            elseif Q[2] <= zero(T) < P[2] && cross < zero(T)
                winding -= 1
            end
        end
        for curve in [curve for curve in path if curve isa Arc]
            φ = atan(-curve.at.y, -curve.at.x)
            result = min(result, hypot(curve(0)...), hypot(curve(1)...))
            if mod(sign(curve.span) * (φ - curve.at.φ), 2π) <= abs(curve.span)
                result = min(result, abs(hypot(curve.at.x, curve.at.y) - curve.r))
            end
            if abs(curve.at.y) < curve.r
                θ = asin(-curve.at.y / curve.r)
                for angle in (θ, π - θ)
                    travel = mod(sign(curve.span) * (angle - curve.at.φ), 2π)
                    direction = sign(curve.span) * cos(angle)
                    if curve.at.x + curve.r * cos(angle) > zero(T) &&
                       ((direction > zero(T) && travel < abs(curve.span)) ||
                        (direction < zero(T) && zero(T) < travel <= abs(curve.span)))
                        winding += direction > zero(T) ? 1 : -1
                    end
                end
            end
        end
        for curve in [curve for curve in path if curve isa EllipticArc]
            radius = max(curve.a, curve.b) + curve.offset
            result = min(result, max(zero(T), hypot(curve.at.x, curve.at.y) - radius))
        end
        inside &= k == 1 ? !iszero(winding) : iszero(winding)
    end
    return inside ? zero(result) : result
end

"""
$(TYPEDSIGNATURES)

Construct four counterclockwise segments from
``\\mathrm{at}((-w/2,-h/2))`` around `rectangle`.
Width, height and coordinates are in \\[m\\].
"""
function Contour(rectangle::Rectangle{T}) where {T}
    w, h, at = rectangle.w, rectangle.h, rectangle.at
    points = [at((-w / 2, -h / 2)), at((w / 2, -h / 2)),
        at((w / 2, h / 2)), at((-w / 2, h / 2))]
    return Contour(Curve{T}[Segment(points[k], points[mod1(k + 1, 4)]) for k in 1:4],
        Vector{Curve{T}}[])
end

"""
$(TYPEDSIGNATURES)

Evaluate the segment at `t` \\[dimensionless\\] in `[0, 1]`:

```math
\\boldsymbol{p}(t) = \\boldsymbol{a} + t(\\boldsymbol{b} - \\boldsymbol{a}),
```

where ``\\boldsymbol{a}`` and ``\\boldsymbol{b}`` are its endpoint coordinates \\[m\\].
Return coordinates in \\[m\\].
"""
function (segment::Segment)(t::Real)
    return segment.a .+ t .* (segment.b .- segment.a)
end

"""
$(TYPEDSIGNATURES)

Evaluate the arc at `t` \\[dimensionless\\] in `[0, 1]`:

```math
\\boldsymbol{p}(t) = \\mathrm{at}\\bigl(r\\cos(t\\,\\mathrm{span}),
    r\\sin(t\\,\\mathrm{span})\\bigr),
```

where ``r`` is its radius \\[m\\], ``\\mathrm{span}`` its signed angular span \\[rad\\], and
``\\mathrm{at}`` rotates by the initial direction \\[rad\\] and translates by the center
coordinates \\[m\\]. Return coordinates in \\[m\\].
"""
function (arc::Arc)(t::Real)
    return arc.at((arc.r * cos(t * arc.span), arc.r * sin(t * arc.span)))
end

"""
$(TYPEDSIGNATURES)

Construct the four counterclockwise quadrant arcs of `circle`, whose radius and center
coordinates are in \\[m\\]. Arc ``k = 0, \\ldots, 3`` starts at the angle ``k \\pi / 2`` and
spans ``\\pi / 2`` \\[rad\\], independent of the placement's rotation. The contour has no
holes.
"""
function Contour(circle::Circle{T}) where {T}
    c = circle.at
    outer = Curve{T}[Arc(Placement(c.x, c.y, T(k * π / 2)), circle.r, T(π / 2))
                     for k in 0:3]
    return Contour(outer, Vector{Curve{T}}[])
end

"""
$(TYPEDSIGNATURES)

Construct the counterclockwise quadrants of `annulus`: its outer radius bounds the outer
path and its inner radius bounds one hole. Radii and center coordinates are in \\[m\\].
Both paths start at angle zero, independent of the placement's rotation \\[rad\\].
"""
function Contour(annulus::Annulus)
    outer = Contour(Circle(annulus.ro, annulus.at)).outer
    hole = Contour(Circle(annulus.ri, annulus.at)).outer
    return Contour(outer, [hole])
end

"""
$(TYPEDSIGNATURES)

Construct the closed path of `sector`, whose radii and center coordinates are in \\[m\\]
and whose angular span and placement rotation are in \\[rad\\]. A partial turn of span
``\\Delta\\phi`` about the rotation ``\\phi`` follows the outer arc from
``\\phi - \\Delta\\phi / 2`` with span ``\\Delta\\phi``, the radial edge, the inner arc
from ``\\phi + \\Delta\\phi / 2`` with span ``-\\Delta\\phi`` and the closing radial
edge. With a
zero inner radius, two segments join the outer arc to the center. A full turn gives the
contour of a circle or annulus.

With displaced sides or fillets, let ``α=span/2``, ``g=apex\\sin α`` and ``ρ=fillet``.
The upper fillet centers ``c`` satisfy ``(\\sin α,-\\cos α)\\cdot c=g+ρ``, with
``|c|=ro-ρ`` at the back and ``ri+ρ`` at the base. The contacts are their projections on
the side and on its circle. The tip center is ``(apex+ρ/\\sin α,0)``. Reflection gives
the lower contacts, and the placement moves every curve to the frame of the cable.
"""
function Contour(sector::Sector{T}) where {T}
    c = sector.at
    if !(iszero(sector.apex) && iszero(sector.fillet))
        α = sector.span / 2
        g = sector.apex * sin(α)
        ρ = sector.fillet
        β = α - asin((g + ρ) / (sector.ro - ρ))
        outer = Arc(Placement(c.x, c.y, c.φ - β), sector.ro, 2β)
        if iszero(ρ)
            if iszero(sector.ri)
                tip = c((sector.apex, zero(sector.apex)))
                return Contour(
                    Curve{T}[outer, Segment(outer(1), tip),
                        Segment(tip, outer(0))], Vector{Curve{T}}[])
            end
            γ = α - asin(g / sector.ri)
            inner = Arc(Placement(c.x, c.y, c.φ + γ), sector.ri, -2γ)
            return Contour(
                Curve{T}[outer, Segment(outer(1), inner(0)), inner,
                    Segment(inner(1), outer(0))],
                Vector{Curve{T}}[])
        end
        # The upper side has inward normal (sin α, -cos α).
        # Its inner parallel meets the back circle at the back fillet's center.
        d = sqrt((sector.ro - ρ)^2 - (g + ρ)^2)
        x = d * cos(α) + (g + ρ) * sin(α)
        y = d * sin(α) - (g + ρ) * cos(α)
        upper = c((x, y))
        lower = c((x, -y))
        upper_back = Arc(Placement(upper..., c.φ + β), ρ, α + π / 2 - β)
        lower_back = Arc(Placement(lower..., c.φ - α - π / 2), ρ, α + π / 2 - β)
        if sector.ri > zero(sector.ri)
            d = sqrt((sector.ri + ρ)^2 - (g + ρ)^2)
            x = d * cos(α) + (g + ρ) * sin(α)
            y = d * sin(α) - (g + ρ) * cos(α)
            upper = c((x, y))
            lower = c((x, -y))
            γ = α - asin((g + ρ) / (sector.ri + ρ))
            upper_base = Arc(Placement(upper..., c.φ + α + π / 2),
                ρ, γ + π / 2 - α)
            inner = Arc(Placement(c.x, c.y, c.φ + γ), sector.ri, -2γ)
            lower_base = Arc(Placement(lower..., c.φ + π - γ),
                ρ, γ + π / 2 - α)
            return Contour(
                Curve{T}[outer, upper_back,
                    Segment(upper_back(1), upper_base(0)), upper_base, inner, lower_base,
                    Segment(lower_base(1), lower_back(0)), lower_back],
                Vector{Curve{T}}[])
        end
        tip = c((sector.apex + ρ / sin(α), zero(ρ)))
        tip_arc = Arc(Placement(tip..., c.φ + α + π / 2), ρ, π - sector.span)
        return Contour(
            Curve{T}[outer, upper_back,
                Segment(upper_back(1), tip_arc(0)), tip_arc,
                Segment(tip_arc(1), lower_back(0)), lower_back],
            Vector{Curve{T}}[])
    end
    if sector.span == 2π
        if iszero(sector.ri)
            return Contour(Circle(sector.ro, c))
        end
        return Contour(Annulus(sector.ri, sector.ro, c))
    end
    outer = Arc(Placement(c.x, c.y, c.φ - sector.span / 2), sector.ro, sector.span)
    if iszero(sector.ri)
        center = (c.x, c.y)
        inward = Segment(outer(1), center)
        outward = Segment(center, outer(0))
        return Contour(Curve{T}[outer, inward, outward], Vector{Curve{T}}[])
    end
    inner = Arc(Placement(c.x, c.y, c.φ + sector.span / 2), sector.ri, -sector.span)
    inward = Segment(outer(1), inner(0))
    outward = Segment(inner(1), outer(0))
    return Contour(Curve{T}[outer, inward, inner, outward], Vector{Curve{T}}[])
end

"""
$(TYPEDSIGNATURES)

Construct the contour of `shape`, whose dimensions and positions are in \\[m\\] and angles
in \\[rad\\]. Keep the contour of `shape.outer` and append `Contour(hole).outer` for each
member of `shape.holes` to its holes. All paths run counterclockwise about the regions they
bound.
"""
function Contour(shape::PerforatedShape)
    contour = Contour(shape.outer)
    holes = copy(contour.holes)
    append!(holes, [Contour(hole).outer for hole in shape.holes])
    return Contour(contour.outer, holes)
end

"""
$(TYPEDSIGNATURES)

Construct the four counterclockwise elliptic quadrants of `ellipse`,
``φ0=kπ/2``, ``span=π/2`` \\[rad\\], for ``k=0,…,3``.
Semi-axes and center coordinates are in \\[m\\].
"""
function Contour(ellipse::Ellipse{T}) where {T}
    return Contour(
        Curve{T}[EllipticArc(ellipse.at, ellipse.a, ellipse.b, T(k * π / 2),
                     T(π / 2)) for k in 0:3],
        Vector{Curve{T}}[])
end

"""
$(TYPEDSIGNATURES)

Construct the placed edges ``Segment(at(v_i),at(v_{i+1}))`` \\[m\\] of `polygon`,
Close the last edge to the first vertex. The contour has no holes.
"""
function Contour(polygon::Polygon{T}) where {T}
    points = [polygon.at(vertex) for vertex in polygon.vertices]
    return Contour(
        Curve{T}[Segment(points[i], points[mod1(i + 1, length(points))])
                 for i in eachindex(points)],
        Vector{Curve{T}}[])
end
