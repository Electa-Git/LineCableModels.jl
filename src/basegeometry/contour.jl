"""
$(TYPEDEF)

Represent a straight curve between two points in the cross-section plane.

$(TYPEDFIELDS)
"""
struct Segment{T <: Real}
    "Coordinates at `t = 0` in \\[m\\]."
    p::NTuple{2, T}
    "Coordinates at `t = 1` in \\[m\\]."
    q::NTuple{2, T}
    function Segment{T}(p::NTuple{2, T}, q::NTuple{2, T}) where {T <: Real}
        return validate(new{T}(p, q))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a segment from the points `p` and `q` \\[m\\], in the type that their
coordinates promote to with `Float64`.

# Errors

- Throws `DomainError` when a coordinate is not finite or the points coincide.
"""
function Segment(p::Tuple{Real, Real}, q::Tuple{Real, Real})
    T = promote_type(Float64, typeof.(p)..., typeof.(q)...)
    return Segment{T}(convert.(T, p), convert.(T, q))
end

function validate(segment::Segment)
    all(isfinite, (segment.p..., segment.q...)) || throw(DomainError(
        (segment.p, segment.q), "Segment coordinates must be finite"))
    segment.p != segment.q || throw(DomainError(
        (segment.p, segment.q), "Segment endpoints must differ"))
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
    isfinite(nominal(arc.offset)) &&
    nominal(nominal(arc.span) > 0 ? dimension(arc) : dimension(reverse(arc))) > 0 ||
        throw(DomainError(arc.offset,
            "EllipticArc.offset must be finite and above minus its smallest radius of curvature"))
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
    return atan(curve.q[2] - curve.p[2], curve.q[1] - curve.p[1])
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
    return (fraction = zero(T), length = hypot((segment.q .- segment.p)...))
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
    dx, dy = curve.q .- curve.p
    φ = atan(-dx, dy)
    return Segment(curve.p .+ d .* (cos(φ), sin(φ)),
        curve.q .+ d .* (cos(φ), sin(φ)))
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

``offset_d=offset+sign(span)d``.
"""
function parallel(curve::EllipticArc, d::Real)
    return EllipticArc(curve.at, curve.a, curve.b, curve.φ0, curve.span,
        curve.offset + sign(curve.span) * d)
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

``offset_d=offset+sign(span)d``.
"""
function parallel(curve::EllipticArc, d::Real, p, q)
    return parallel(curve, d)
end

"""
$(TYPEDSIGNATURES)

Return the turn of the parallel of `face` at the distance `d` \\[m\\] between the points `p`
and `q`. The fraction is zero, and the signed length \\[m\\] is ``(q-p)⋅u`` for the unit
direction ``u`` of `face`.
"""
function turn(face::Segment, d::Real, p, q)
    u = (face.q .- face.p) ./ hypot((face.q .- face.p)...)
    return (fraction = zero(face.p[1] + d),
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

Return the turn of the parallel of `face` at the distance `d` \\[m\\] between the points `p`
and `q`: the turn of its whole parallel. The junctions of an elliptic arc are tangent or
diverge. Its piece runs from the start of its parallel to its end.
"""
turn(face::EllipticArc, d::Real, p, q) = turn(parallel(face, d))

"""
$(TYPEDSIGNATURES)

Return the widest common angular range \\[rad\\] of two arcs of one center, or of two
elliptic arcs of one ellipse. The range is an interval of the angle of `a`, read
counterclockwise. With ``α`` the start of `a` in that sense and ``β`` the start of `b`, it
is the widest of ``[α, α+|span_a|] ∩ [β+2πk, β+2πk+|span_b|]`` for ``k ∈ \\{-1, 0, 1\\}``. An
elliptic arc reads its parametric angle. For two disjoint arcs the interval is empty, and
its endpoints are the two ends that bound the gap. The caller checks that the arcs have one
center or one ellipse.
"""
function overlap(a::Arc, b::Arc)
    α = mod(a.at.φ + min(a.span, zero(a.span)), 2π)
    β = mod(b.at.φ + min(b.span, zero(b.span)), 2π)
    A = ClosedInterval(α, α + abs(a.span))
    return argmax(I -> nominal(rightendpoint(I) - leftendpoint(I)),
        (A ∩ ClosedInterval(β + k * 2π, β + k * 2π + abs(b.span)) for k in (-1, 0, 1)))
end

function overlap(a::EllipticArc, b::EllipticArc)
    α = mod(a.φ0 + min(a.span, zero(a.span)), 2π)
    β = mod(b.φ0 + min(b.span, zero(b.span)), 2π)
    A = ClosedInterval(α, α + abs(a.span))
    return argmax(I -> nominal(rightendpoint(I) - leftendpoint(I)),
        (A ∩ ClosedInterval(β + k * 2π, β + k * 2π + abs(b.span)) for k in (-1, 0, 1)))
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
    u, v = curve.q .- curve.p, other.q .- other.p
    w = other.p .- curve.p
    z = zero(u[1] + v[1])
    best = argmin(candidate -> nominal(candidate[3]),
        map(((s, t),) -> (s, t, hypot((curve(s) .- other(t))...)),
            ((z, z), (z, one(z)), (one(z), z), (one(z), one(z)))))
    cross = u[1] * v[2] - u[2] * v[1]
    if iszero(nominal(cross))
        L = hypot(u...)
        n = u ./ L
        α,
        β = w[1] * n[1] + w[2] * n[2],
        (other.q[1] - curve.p[1]) * n[1] + (other.q[2] - curve.p[2]) * n[2]
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
            d = z
            best = argmin(candidate -> nominal(candidate[3]), (best, (s, t, d)))
        end
    end
    for s in (z, one(z))
        p = curve(s)
        t = ((p[1] - other.p[1]) * v[1] + (p[2] - other.p[2]) * v[2]) /
            (v[1]^2 + v[2]^2)
        t = nominal(t) <= 0 ? z : nominal(t) >= 1 ? one(z) : t
        d = hypot((p .- other(t))...)
        best = argmin(candidate -> nominal(candidate[3]), (best, (s, t, d)))
    end
    for t in (z, one(z))
        q = other(t)
        s = ((q[1] - curve.p[1]) * u[1] + (q[2] - curve.p[2]) * u[2]) /
            (u[1]^2 + u[2]^2)
        s = nominal(s) <= 0 ? z : nominal(s) >= 1 ? one(z) : s
        d = hypot((curve(s) .- q)...)
        best = argmin(candidate -> nominal(candidate[3]), (best, (s, t, d)))
    end
    s, t, d = best
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
    u = curve.q .- curve.p
    z = zero(u[1] + other.r)
    best = argmin(candidate -> nominal(candidate[3]),
        map(((s, t),) -> (s, t, hypot((curve(s) .- other(t))...)),
            ((z, z), (z, one(z)), (one(z), z), (one(z), one(z)))))
    L = u[1]^2 + u[2]^2
    for t in (z, one(z))
        q = other(t)
        s = ((q[1] - curve.p[1]) * u[1] + (q[2] - curve.p[2]) * u[2]) / L
        s = nominal(s) <= 0 ? z : nominal(s) >= 1 ? one(z) : s
        d = hypot((curve(s) .- q)...)
        best = argmin(candidate -> nominal(candidate[3]), (best, (s, t, d)))
    end
    for s in (z, one(z))
        p = curve(s)
        φ = atan(p[2] - other.at.y, p[1] - other.at.x)
        travel = mod(sign(other.span) * (φ - other.at.φ), 2π)
        if nominal(travel) <= abs(nominal(other.span))
            t = travel / abs(other.span)
            d = hypot((p .- other(t))...)
            best = argmin(candidate -> nominal(candidate[3]), (best, (s, t, d)))
        end
    end
    h = ((other.at.x - curve.p[1]) * u[1] +
         (other.at.y - curve.p[2]) * u[2]) / L
    if 0 <= nominal(h) <= 1
        s = h
        p = curve(s)
        for φ in (atan(-u[1], u[2]), atan(u[1], -u[2]))
            travel = mod(sign(other.span) * (φ - other.at.φ), 2π)
            if nominal(travel) <= abs(nominal(other.span))
                t = travel / abs(other.span)
                d = hypot((p .- other(t))...)
                best = argmin(candidate -> nominal(candidate[3]), (best, (s, t, d)))
            end
        end
    end
    w = curve.p .- (other.at.x, other.at.y)
    Δ = h^2 - (w[1]^2 + w[2]^2 - other.r^2) / L
    if nominal(Δ) >= 0
        for s in (h - sqrt(max(Δ, zero(Δ))), h + sqrt(max(Δ, zero(Δ))))
            if 0 <= nominal(s) <= 1
                p = curve(s)
                φ = atan(p[2] - other.at.y, p[1] - other.at.x)
                travel = mod(sign(other.span) * (φ - other.at.φ), 2π)
                if nominal(travel) <= abs(nominal(other.span))
                    t = travel / abs(other.span)
                    d = z
                    best = argmin(candidate -> nominal(candidate[3]), (best, (s, t, d)))
                end
            end
        end
    end
    s, t, d = best
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
middle of their widest common interval, at ``|r_1-r_2|`` less the distance of their
centers, or else their nearest ends.
"""
function clearance(curve::Arc, other::Arc)
    z = zero(curve.r + other.r)
    best = argmin(candidate -> nominal(candidate[3]),
        map(((s, t),) -> (s, t, hypot((curve(s) .- other(t))...)),
            ((z, z), (z, one(z)), (one(z), z), (one(z), one(z)))))
    w = (other.at.x - curve.at.x, other.at.y - curve.at.y)
    L = hypot(w...)
    if nominal(curve.at.x) == nominal(other.at.x) &&
       nominal(curve.at.y) == nominal(other.at.y)
        α = mod(curve.at.φ + min(curve.span, zero(curve.span)), 2π)
        β = mod(other.at.φ + min(other.span, zero(other.span)), 2π)
        interval = overlap(curve, other)
        if !isempty(interval)
            φ = (leftendpoint(interval) + rightendpoint(interval)) / 2
            s = nominal(curve.span) > 0 ? (φ - α) / abs(curve.span) :
                (α + abs(curve.span) - φ) / abs(curve.span)
            ψ = mod(φ - β, 2π)
            t = nominal(other.span) > 0 ? ψ / abs(other.span) :
                (abs(other.span) - ψ) / abs(other.span)
            s = nominal(φ) == nominal(α) ?
                (nominal(curve.span) > 0 ? z : one(z)) :
                nominal(φ) == nominal(α + abs(curve.span)) ?
                (nominal(curve.span) > 0 ? one(z) : z) : s
            t = nominal(ψ) == 0 ?
                (nominal(other.span) > 0 ? z : one(z)) :
                nominal(ψ) == nominal(abs(other.span)) ?
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
                d = hypot((p .- other(t))...)
                best = argmin(candidate -> nominal(candidate[3]), (best, (s, t, d)))
            end
        end
        for t in (z, one(z))
            q = other(t)
            φ = atan(q[2] - curve.at.y, q[1] - curve.at.x)
            travel = mod(sign(curve.span) * (φ - curve.at.φ), 2π)
            if nominal(travel) <= abs(nominal(curve.span))
                s = travel / abs(curve.span)
                d = hypot((curve(s) .- q)...)
                best = argmin(candidate -> nominal(candidate[3]), (best, (s, t, d)))
            end
        end
        for φ in (atan(w[2], w[1]), atan(w[2], w[1]) + π),
            ψ in (atan(w[2], w[1]), atan(w[2], w[1]) + π)

            travel₁ = mod(sign(curve.span) * (φ - curve.at.φ), 2π)
            travel₂ = mod(sign(other.span) * (ψ - other.at.φ), 2π)
            if nominal(travel₁) <= abs(nominal(curve.span)) &&
               nominal(travel₂) <= abs(nominal(other.span))
                s, t = travel₁ / abs(curve.span), travel₂ / abs(other.span)
                d = hypot((curve(s) .- other(t))...)
                best = argmin(candidate -> nominal(candidate[3]), (best, (s, t, d)))
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
                    d = z
                    best = argmin(candidate -> nominal(candidate[3]), (best, (s, t, d)))
                end
            end
        end
    end
    s, t, d = best
    return (curve(s), other(t), d, s, t)
end

"""
$(TYPEDSIGNATURES)

Return `(p, q, d, s, t)` for the nearest points `p` of `curve` and `q` of `other`, as
`clearance(::Segment, ::Segment)`. The candidates are the four end pairs and each end with its
nearest point on the other curve. They also include the apsides of `other` along the
segment's normal, with their feet on the segment. Between those parameters and
the ends, the signed distance of `other` from the segment's line is monotone, and `bisect`
finds each change of its sign. A change whose point projects onto the segment is a crossing,
at ``d=0``.
"""
function clearance(curve::Segment, other::EllipticArc)
    u = curve.q .- curve.p
    L = u[1]^2 + u[2]^2
    z = zero(u[1] + other.a)
    best = argmin(candidate -> nominal(candidate[3]),
        map(((s, t),) -> (s, t, hypot((curve(s) .- other(t))...)),
            ((z, z), (z, one(z)), (one(z), z), (one(z), one(z)))))
    for s in (z, one(z))
        t = apsides(other, curve(s)).near
        best = argmin(candidate -> nominal(candidate[3]),
            (best, (s, t, hypot((curve(s) .- other(t))...))))
    end
    near, far = apsides(other, atan(u[2], u[1]) + π / 2)
    for t in (z, one(z), near, far)
        s = apsides(curve, other(t)).near
        best = argmin(candidate -> nominal(candidate[3]),
            (best, (s, t, hypot((curve(s) .- other(t))...))))
    end
    knots = nominal(near) <= nominal(far) ? (z, near, far, one(z)) : (z, far, near, one(z))
    for k in 1:3
        q₁, q₂ = other(knots[k]), other(knots[k + 1])
        g₁ = u[1] * (q₁[2] - curve.p[2]) - u[2] * (q₁[1] - curve.p[1])
        g₂ = u[1] * (q₂[2] - curve.p[2]) - u[2] * (q₂[1] - curve.p[1])
        nominal(g₁) * nominal(g₂) > 0 && continue
        t = iszero(nominal(g₁)) ? knots[k] :
            iszero(nominal(g₂)) ? knots[k + 1] :
            argmin(
            t -> (xy = other(t);
                abs(nominal(u[1] * (xy[2] - curve.p[2]) - u[2] * (xy[1] - curve.p[1])))),
            bisect(
                t -> (xy = other(t);
                    nominal(sign(u[1] * (xy[2] - curve.p[2]) -
                                 u[2] * (xy[1] - curve.p[1]))) ==
                    nominal(sign(g₁))),
                knots[k], knots[k + 1]))
        q = other(t)
        s = ((q[1] - curve.p[1]) * u[1] + (q[2] - curve.p[2]) * u[2]) / L
        if 0 <= nominal(s) <= 1
            best = argmin(candidate -> nominal(candidate[3]), (best, (s, t, z)))
        end
    end
    s, t, d = best
    return (curve(s), other(t), d, s, t)
end

"""
$(TYPEDSIGNATURES)

Return `(p, q, d, s, t)` for the nearest points `p` of `curve` and `q` of `other`, as
`clearance(::Segment, ::Segment)`. The candidates start with the end pairs and each end's
nearest point on the other curve. The method then cuts `other` where it crosses the rays from
the arc's center ``O`` through the arc's ends.

On a piece inside the arc's angle, the distance to the arc is
``\\lvert\\lvert q-O\\rvert-r\\rvert``. There, the piece's nearest apsis about ``O`` gives it
outside the circle and its farthest apsis inside. Otherwise `bisect` finds a crossing
between them, at ``d=0``. Outside the angle, the arc's ends are nearest.
"""
function clearance(curve::Arc, other::EllipticArc)
    z = zero(curve.r + other.a)
    O = (curve.at.x, curve.at.y)
    best = argmin(candidate -> nominal(candidate[3]),
        map(((s, t),) -> (s, t, hypot((curve(s) .- other(t))...)),
            ((z, z), (z, one(z)), (one(z), z), (one(z), one(z)))))
    for s in (z, one(z))
        t = apsides(other, curve(s)).near
        best = argmin(candidate -> nominal(candidate[3]),
            (best, (s, t, hypot((curve(s) .- other(t))...))))
    end
    for t in (z, one(z))
        s = apsides(curve, other(t)).near
        best = argmin(candidate -> nominal(candidate[3]),
            (best, (s, t, hypot((curve(s) .- other(t))...))))
    end
    cuts = [z, one(z)]
    for α in (curve.at.φ, curve.at.φ + curve.span)
        near, far = apsides(other, α + π / 2)
        knots = nominal(near) <= nominal(far) ? (z, near, far, one(z)) :
                (z, far, near, one(z))
        for k in 1:3
            q₁, q₂ = other(knots[k]), other(knots[k + 1])
            g₁ = cos(α) * (q₁[2] - O[2]) - sin(α) * (q₁[1] - O[1])
            g₂ = cos(α) * (q₂[2] - O[2]) - sin(α) * (q₂[1] - O[1])
            nominal(g₁) * nominal(g₂) > 0 && continue
            t = iszero(nominal(g₁)) ? knots[k] :
                iszero(nominal(g₂)) ? knots[k + 1] :
                argmin(
                t -> (xy = other(t);
                    abs(nominal(cos(α) * (xy[2] - O[2]) - sin(α) * (xy[1] - O[1])))),
                bisect(
                    t -> (xy = other(t);
                        nominal(sign(cos(α) * (xy[2] - O[2]) - sin(α) * (xy[1] - O[1]))) ==
                        nominal(sign(g₁))),
                    knots[k], knots[k + 1]))
            q = other(t)
            nominal(cos(α) * (q[1] - O[1]) + sin(α) * (q[2] - O[2])) > 0 && push!(cuts, t)
        end
    end
    sort!(cuts; by = nominal)
    for k in 1:(length(cuts) - 1)
        t₁, t₂ = cuts[k], cuts[k + 1]
        nominal(t₂) > nominal(t₁) || continue
        q = other((t₁ + t₂) / 2)
        travel = mod(sign(curve.span) * (atan(q[2] - O[2], q[1] - O[1]) - curve.at.φ), 2π)
        nominal(travel) <= abs(nominal(curve.span)) || continue
        piece = EllipticArc(other.at, other.a, other.b, other.φ0 + t₁ * other.span,
            (t₂ - t₁) * other.span, other.offset)
        near, far = apsides(piece, O)
        near, far = t₁ + near * (t₂ - t₁), t₁ + far * (t₂ - t₁)
        if nominal(hypot((other(near) .- O)...)) > nominal(curve.r)
            s = apsides(curve, other(near)).near
            best = argmin(candidate -> nominal(candidate[3]),
                (best, (s, near, hypot((curve(s) .- other(near))...))))
        elseif nominal(hypot((other(far) .- O)...)) < nominal(curve.r)
            s = apsides(curve, other(far)).near
            best = argmin(candidate -> nominal(candidate[3]),
                (best, (s, far, hypot((curve(s) .- other(far))...))))
        else
            t = nominal(hypot((other(near) .- O)...)) == nominal(curve.r) ? near :
                nominal(hypot((other(far) .- O)...)) == nominal(curve.r) ? far :
                argmin(t -> abs(nominal(hypot((other(t) .- O)...) - curve.r)),
                bisect(t -> nominal(hypot((other(t) .- O)...)) < nominal(curve.r), near, far))
            best = argmin(candidate -> nominal(candidate[3]),
                (best, (apsides(curve, other(t)).near, t, z)))
        end
    end
    s, t, d = best
    return (curve(s), other(t), d, s, t)
end

"""
$(TYPEDSIGNATURES)

Return `clearance(other, curve)` with its two points and its two parameters swapped,
`(q, p, d, t, s)`.
"""
function clearance(curve::EllipticArc, other::Union{Segment, Arc})
    q, p, d, t, s = clearance(other, curve)
    return (p, q, d, s, t)
end

"""
$(TYPEDSIGNATURES)

Return `(p, q, d, s, t)` for the nearest points `p` of `curve` and `q` of `other`. Their
distance is `d` \\[m\\], and their parameters give ``p = curve(s)`` and ``q = other(t)``.
Each parameter lies in ``[0, 1]`` and is exactly 0 or 1 at an end.

Arcs of one ellipse take the middle of their widest common interval, at ``|o_1-o_2|`` for
their offsets ``o_1`` and ``o_2``. Without one, the candidates are the end pairs, each end
with its nearest point on the other arc, and the pairs along the two axes.

Arcs of different ellipses throw `DomainError`. So do arcs of one ellipse with an offset at or below minus its least radius
of curvature, ``\\min(b^2/a,a^2/b)``. There its parallels are no longer level sets of the
distance from it.
"""
function clearance(curve::EllipticArc, other::EllipticArc)
    (nominal(curve.at.x), nominal(curve.at.y), nominal(curve.at.φ),
        nominal(curve.a), nominal(curve.b)) ==
    (nominal(other.at.x), nominal(other.at.y), nominal(other.at.φ),
        nominal(other.a), nominal(other.b)) ||
        throw(DomainError((curve, other), "the clearance of two elliptic arcs needs one ellipse"))
    reach = min(curve.b^2 / curve.a, curve.a^2 / curve.b)
    nominal(min(curve.offset, other.offset)) > -nominal(reach) ||
        throw(DomainError((curve, other),
            "the clearance of two elliptic arcs of one ellipse needs their offsets above minus its least radius of curvature"))
    z = zero(curve.a + other.a)
    best = argmin(candidate -> nominal(candidate[3]),
        map(((s, t),) -> (s, t, hypot((curve(s) .- other(t))...)),
            ((z, z), (z, one(z)), (one(z), z), (one(z), one(z)))))
    α = mod(curve.φ0 + min(curve.span, zero(curve.span)), 2π)
    β = mod(other.φ0 + min(other.span, zero(other.span)), 2π)
    interval = overlap(curve, other)
    if !isempty(interval)
        φ = (leftendpoint(interval) + rightendpoint(interval)) / 2
        s = nominal(curve.span) > 0 ? (φ - α) / abs(curve.span) :
            (α + abs(curve.span) - φ) / abs(curve.span)
        ψ = mod(φ - β, 2π)
        t = nominal(other.span) > 0 ? ψ / abs(other.span) :
            (abs(other.span) - ψ) / abs(other.span)
        s = nominal(φ) == nominal(α) ?
            (nominal(curve.span) > 0 ? z : one(z)) :
            nominal(φ) == nominal(α + abs(curve.span)) ?
            (nominal(curve.span) > 0 ? one(z) : z) : s
        t = nominal(ψ) == 0 ?
            (nominal(other.span) > 0 ? z : one(z)) :
            nominal(ψ) == nominal(abs(other.span)) ?
            (nominal(other.span) > 0 ? one(z) : z) : t
        return (curve(s), other(t), abs(curve.offset - other.offset), s, t)
    end
    for s in (z, one(z))
        t = apsides(other, curve(s)).near
        best = argmin(candidate -> nominal(candidate[3]),
            (best, (s, t, hypot((curve(s) .- other(t))...))))
    end
    for t in (z, one(z))
        s = apsides(curve, other(t)).near
        best = argmin(candidate -> nominal(candidate[3]),
            (best, (s, t, hypot((curve(s) .- other(t))...))))
    end
    for θ in (z, z + π / 2, z + π, z + 3π / 2)
        s = mod(sign(curve.span) * (θ - curve.φ0), 2π)
        t = mod(sign(other.span) * (θ + π - other.φ0), 2π)
        if nominal(s) <= abs(nominal(curve.span)) && nominal(t) <= abs(nominal(other.span))
            s, t = s / abs(curve.span), t / abs(other.span)
            best = argmin(candidate -> nominal(candidate[3]),
                (best, (s, t, hypot((curve(s) .- other(t))...))))
        end
    end
    s, t, d = best
    return (curve(s), other(t), d, s, t)
end

"""
$(TYPEDSIGNATURES)

Place the endpoints of `curve` \\[m\\] by `at`, whose translation is in \\[m\\]
and rotation in \\[rad\\].
"""
function resolve(at::Placement, curve::Segment)
    return Segment(at(curve.p), at(curve.q))
end

"""
$(TYPEDSIGNATURES)

Place the center and initial direction of `curve` by `at`. Coordinates and radius
are in \\[m\\], and angles in \\[rad\\].
"""
function resolve(at::Placement, curve::Arc)
    return Arc(at * curve.at, curve.r, curve.span)
end

"""
$(TYPEDSIGNATURES)

Place the center and axes of `curve` by `at`. Coordinates, semi-axes and normal
offset are in \\[m\\], and angles in \\[rad\\].
"""
function resolve(at::Placement, curve::EllipticArc)
    return EllipticArc(at * curve.at, curve.a, curve.b, curve.φ0, curve.span, curve.offset)
end

"""
$(TYPEDSIGNATURES)

Map the endpoint coordinates \\[m\\] of `curve` by `f`.
"""
function Base.map(f, curve::Segment)
    return Segment(map(f, curve.p), map(f, curve.q))
end

"""
$(TYPEDSIGNATURES)

Map the center coordinates and radius \\[m\\], and angles \\[rad\\], of `curve` by `f`.
"""
function Base.map(f, curve::Arc)
    return Arc(map(f, curve.at), f(curve.r), f(curve.span))
end

"""
$(TYPEDSIGNATURES)

Map the coordinates, semi-axes and normal offset \\[m\\], and angles \\[rad\\],
of `curve` by `f`.
"""
function Base.map(f, curve::EllipticArc)
    return EllipticArc(map(f, curve.at), f(curve.a), f(curve.b),
        f(curve.φ0), f(curve.span), f(curve.offset))
end

"""
$(TYPEDSIGNATURES)

Return `curve` traversed in the opposite direction. Coordinates and lengths are in \\[m\\],
and angles are in \\[rad\\].
"""
function Base.reverse(curve::Segment)
    return Segment(curve.q, curve.p)
end

"""
$(TYPEDSIGNATURES)

Return the distance across the metal of `curve`, the metal on its left \\[m\\]: infinite,
since the inward parallel of a straight curve never degenerates.
"""
function dimension(curve::Segment)
    return oftype(curve.p[1], Inf)
end

"""
$(TYPEDSIGNATURES)

Return `curve` traversed in the opposite direction. Coordinates and lengths are in \\[m\\],
and angles are in \\[rad\\].
"""
function Base.reverse(curve::Arc)
    return Arc(Placement(curve.at.x, curve.at.y, curve.at.φ + curve.span), curve.r, -curve.span)
end

"""
$(TYPEDSIGNATURES)

Return the distance across the metal of `curve`, the metal on its left \\[m\\]. It is twice
the radius when the arc turns toward its left, where its inward parallel shrinks to its
center. It is infinite when the arc turns away.
"""
function dimension(curve::Arc)
    return nominal(curve.span) > 0 ? 2 * curve.r : oftype(curve.r, Inf)
end

"""
$(TYPEDSIGNATURES)

Return `curve` traversed in the opposite direction. Coordinates and lengths are in \\[m\\],
and angles are in \\[rad\\].
"""
function Base.reverse(curve::EllipticArc)
    return EllipticArc(
        curve.at, curve.a, curve.b, curve.φ0 + curve.span, -curve.span, curve.offset)
end

"""
$(TYPEDSIGNATURES)

Return the distance across the metal of `curve`, the metal on its left \\[m\\]. It is twice
the smallest radius of curvature of its parallel over its interval when the arc turns
toward its left. It is infinite when the arc turns away. For its normal offset ``d``, the
radius of the parallel equals ``ρ(θ)+d``. The arc's own radius of curvature equals
``ρ(θ)=(a^2\\sin^2θ+b^2\\cos^2θ)^{3/2}/(ab)``.
"""
function dimension(curve::EllipticArc)
    nominal(curve.span) > 0 || return oftype(curve.a, Inf)
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
    return 2 * (ρ_min + curve.offset)
end

"""
$(TYPEDSIGNATURES)

Return the parameters `(lo, hi)`, adjacent in nominal value, between which `predicate`
changes from true to false. The bracket starts where `predicate` holds at `lo` and fails at
`hi`. Each step halves it at its midpoint, and `predicate` decides which half keeps the
change.
"""
function bisect(predicate, lo::Real, hi::Real)
    lo, hi = promote(lo, hi)
    while true
        mid = (lo + hi) / 2
        (nominal(mid) == nominal(lo) || nominal(mid) == nominal(hi)) && break
        if predicate(mid)
            lo = mid
        else
            hi = mid
        end
    end
    return (lo, hi)
end

"""
$(TYPEDSIGNATURES)

Return the parameters `(near, far)` of the points of `curve` nearest to and farthest from
the point `m` \\[m\\], its ends included. A segment's nearest point is the foot of `m`,
clamped to the segment. An arc's lie in the directions from its center toward `m` and away
from it, each clamped to the span, or at an end.

The points of an elliptic arc lie where its distance from `m` is stationary. In the
ellipse's frame, with ``(x, y)`` the center less `m`, the condition
``(b^2-a^2)\\sin θ\\cos θ-a x\\sin θ+b y\\cos θ=0`` holds there. Bracket its roots between
the quadrant ends and the angles of ``\\tan^3 θ=-b y/(a x)``, bisect each bracket, and
compare the roots in the arc with its ends. A regular normal offset keeps the stationary
parameters.
"""
function apsides(curve::Segment{T}, m) where {T}
    P = curve.p .- m
    dx, dy = curve.q .- curve.p
    t = clamp(-(P[1] * dx + P[2] * dy) / (dx^2 + dy^2), zero(T), one(T))
    return (
        near = t, far = argmax(t -> nominal(hypot((curve(t) .- m)...)), (zero(T), one(T))))
end

function apsides(curve::Arc{T}, m) where {T}
    φ = atan(m[2] - curve.at.y, m[1] - curve.at.x)
    toward = clamp(mod(sign(curve.span) * (φ - curve.at.φ), 2π) / abs(curve.span), zero(T),
        one(T))
    away = clamp(
        mod(sign(curve.span) * (φ + π - curve.at.φ), 2π) / abs(curve.span), zero(T),
        one(T))
    return (
        near = argmin(t -> nominal(hypot((curve(t) .- m)...)), (zero(T), one(T), toward)),
        far = argmax(t -> nominal(hypot((curve(t) .- m)...)), (zero(T), one(T), away)))
end

function apsides(curve::EllipticArc{T}, m) where {T}
    c, s = cos(curve.at.φ), sin(curve.at.φ)
    X, Y = curve.at.x - m[1], curve.at.y - m[2]
    x = X * c + Y * s
    y = -X * s + Y * c
    a, b = curve.a, curve.b
    θ = !iszero(nominal(x)) && !iszero(nominal(y)) ? mod(atan(cbrt(-b * y / (a * x))), π) :
        T(π / 2)
    angles = nominal(θ) < π / 2 ?
             (zero(T), θ, T(π / 2), T(π), θ + π, T(3π / 2), T(2π)) :
             (zero(T), T(π / 2), θ, T(π), T(3π / 2), θ + π, T(2π))
    near = argmin(t -> nominal(hypot((curve(t) .- m)...)), (zero(T), one(T)))
    far = argmax(t -> nominal(hypot((curve(t) .- m)...)), (zero(T), one(T)))
    for k in 1:6
        left = (b^2 - a^2) * sin(angles[k]) * cos(angles[k]) - a * x * sin(angles[k]) +
               b * y * cos(angles[k])
        right = (b^2 - a^2) * sin(angles[k + 1]) * cos(angles[k + 1]) -
                a * x * sin(angles[k + 1]) + b * y * cos(angles[k + 1])
        nominal(left) * nominal(right) > 0 && continue
        for ϑ in bisect(
            ϑ -> nominal(sign((b^2 - a^2) * sin(ϑ) * cos(ϑ) - a * x * sin(ϑ) +
                              b * y * cos(ϑ))) == nominal(sign(left)),
            angles[k], angles[k + 1])
            travel = mod(sign(curve.span) * (ϑ - curve.φ0), 2π)
            if nominal(travel) <= abs(nominal(curve.span))
                t = travel / abs(curve.span)
                near = argmin(t -> nominal(hypot((curve(t) .- m)...)), (near, t))
                far = argmax(t -> nominal(hypot((curve(t) .- m)...)), (far, t))
            end
        end
    end
    return (near = near, far = far)
end

"""
$(TYPEDSIGNATURES)

Return the parameters `(near, far)` of the points of `curve` with the greatest and the least
projection on the direction `φ` \\[rad\\]. Its ends are candidates. They are its apsides about the
point at infinity along `φ`.

At the distance ``R`` along `φ`, a point's distance from the curve is ``R`` less the
projection, to order ``1/R``. A segment's are its ends. An arc's lie
in the directions `φ` and ``φ+π`` from its center, each clamped to the span, or at an end. An
elliptic arc's lie at the parametric angles ``θ=\\operatorname{atan}(b\\sin ψ,a\\cos ψ)`` and
``θ+π``, with ``ψ=φ-at.φ``, where its outward normal points along `φ` and against it. A
regular normal offset keeps them.
"""
function apsides(curve::Segment{T}, φ::Real) where {T}
    c, s = cos(φ), sin(φ)
    return (
        near = argmax(t -> (xy = curve(t); nominal(xy[1] * c + xy[2] * s)), (
            zero(T), one(T))),
        far = argmin(t -> (xy = curve(t); nominal(xy[1] * c + xy[2] * s)), (
            zero(T), one(T))))
end

function apsides(curve::Arc{T}, φ::Real) where {T}
    c, s = cos(φ), sin(φ)
    toward = clamp(mod(sign(curve.span) * (φ - curve.at.φ), 2π) / abs(curve.span), zero(T),
        one(T))
    away = clamp(mod(sign(curve.span) * (φ + π - curve.at.φ), 2π) / abs(curve.span),
        zero(T), one(T))
    return (
        near = argmax(t -> (xy = curve(t); nominal(xy[1] * c + xy[2] * s)),
            (zero(T), one(T), toward)),
        far = argmin(t -> (xy = curve(t); nominal(xy[1] * c + xy[2] * s)),
            (zero(T), one(T), away)))
end

function apsides(curve::EllipticArc{T}, φ::Real) where {T}
    c, s = cos(φ), sin(φ)
    ψ = φ - curve.at.φ
    θ = atan(curve.b * sin(ψ), curve.a * cos(ψ))
    toward = clamp(mod(sign(curve.span) * (θ - curve.φ0), 2π) / abs(curve.span), zero(T),
        one(T))
    away = clamp(mod(sign(curve.span) * (θ + π - curve.φ0), 2π) / abs(curve.span),
        zero(T), one(T))
    return (
        near = argmax(t -> (xy = curve(t); nominal(xy[1] * c + xy[2] * s)),
            (zero(T), one(T), toward)),
        far = argmin(t -> (xy = curve(t); nominal(xy[1] * c + xy[2] * s)),
            (zero(T), one(T), away)))
end

"""
$(TYPEDSIGNATURES)

Return the distance \\[m\\] of the farther endpoint of `curve` from the origin.
"""
function support(curve::Segment)
    return max(distance(curve.p), distance(curve.q))
end

"""
$(TYPEDSIGNATURES)

Return the farthest distance \\[m\\] from the origin on `curve`: the center's
distance plus its radius when that direction lies in its span, otherwise an endpoint.
For a center at the origin, every point lies at the radius.
"""
function support(curve::Arc)
    iszero(distance(curve.at)) && return curve.r
    φ = atan(curve.at.y, curve.at.x)
    if nominal(mod(sign(curve.span) * (φ - curve.at.φ), 2π)) <= abs(nominal(curve.span))
        return distance(curve.at) + curve.r
    end
    return max(distance(curve(0)), distance(curve(1)))
end

"""
$(TYPEDSIGNATURES)

Return the distance of the farthest point of `curve` from the origin, at its normal offset,
in \\[m\\]: the far point of its `apsides` about the origin.
"""
function support(curve::EllipticArc{T}) where {T}
    return distance(curve(apsides(curve, (zero(T), zero(T))).far))
end

"""
A curve of a contour is a [`Segment`](@ref) an [`Arc`](@ref) or an [`EllipticArc`](@ref) in the floating type `T`.
"""
const Curve{T} = Union{Segment{T}, Arc{T}, EllipticArc{T}}

"""
$(TYPEDSIGNATURES)

Return the largest projection of `curve` on the direction `φ` \\[rad\\], in \\[m\\]: the
projection of its point at the near parameter of `apsides(curve, φ)`.
"""
function support(curve::Curve, φ::Real)
    x, y = curve(apsides(curve, φ).near)
    return x * cos(φ) + y * sin(φ)
end

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

Place `contour` by `at`: each curve of its outer path and of its holes, in order, by
``P\\mapsto\\mathrm{at}(P)`` for segment ends and ``c\\mapsto\\mathrm{at}\\,c`` for arc
placements. Coordinates are in \\[m\\] and rotations in \\[rad\\]. Keep each elliptic arc's
normal offset \\[m\\].
"""
function resolve(at::Placement, contour::Contour{T}) where {T}
    S = promote_type(eltype(at), T)
    return Contour(Curve{S}[resolve(at, curve) for curve in contour.outer],
        Vector{Curve{S}}[Curve{S}[resolve(at, curve) for curve in path]
                         for path in contour.holes])
end

"""
$(TYPEDSIGNATURES)

Return the end of the parallel of `a` and the start of the parallel of `b` at the distance
`d` \\[m\\]. The curve `a` ends at ``j = a(1)``, where `b` starts. At a tangent junction,
both points lie at the distance `d` from ``j`` along the outward normal. Where the two
parallels converge, both points are their `vertex` nearest ``j``, which an elliptic arc
refuses. Where they diverge, each point is the end of its own parallel.
"""
function parallel(a::Curve, b::Curve, d::Real)
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
        p = vertex(parallel(a, d), parallel(b, d), j)
        return (p, p)
    else
        return (parallel(a, d)(1), parallel(b, d)(0))
    end
end

"""
$(TYPEDSIGNATURES)

Return the crossing of the curves `A` and `B` nearest the point `j`, in \\[m\\]. Segments
cross where their lines meet. A segment and an arc cross where the line meets the circle,
and two arcs cross where their circles meet. Of the two crossings of a circle, the result is
the one nearest `j`.

A line and a circle that have ceased to meet give the foot of the normal from the center to
the line. Circles that have ceased to meet give the point where their radical axis crosses
the line of their centers. An elliptic arc's parallel is no conic, and the method throws
`DomainError` for it.
"""
function vertex(A::Segment, B::Segment, j)
    U, V = A.q .- A.p, B.q .- B.p
    w = B.p .- A.p
    t = (w[1] * V[2] - w[2] * V[1]) / (U[1] * V[2] - U[2] * V[1])
    return A.p .+ t .* U
end

function vertex(A::Segment, B::Arc, j)
    U = A.q .- A.p
    w = A.p .- (B.at.x, B.at.y)
    L = U[1]^2 + U[2]^2
    h = -(w[1] * U[1] + w[2] * U[2]) / L
    Δ = h^2 - (w[1]^2 + w[2]^2 - B.r^2) / L
    δ = sqrt(max(Δ, zero(Δ)))
    points = (A.p .+ (h - δ) .* U, A.p .+ (h + δ) .* U)
    return argmin(xy -> nominal(hypot((xy .- j)...)), points)
end

vertex(A::Arc, B::Segment, j) = vertex(B, A, j)

function vertex(A::Arc, B::Arc, j)
    w = (B.at.x - A.at.x, B.at.y - A.at.y)
    L = hypot(w...)
    U = w ./ L
    h = (A.r^2 - B.r^2 + L^2) / (2L)
    Δ = A.r^2 - h^2
    δ = sqrt(max(Δ, zero(Δ)))
    center = (A.at.x, A.at.y) .+ h .* U
    points = (center .+ δ .* (-U[2], U[1]), center .- δ .* (-U[2], U[1]))
    return argmin(xy -> nominal(hypot((xy .- j)...)), points)
end

function vertex(A::Curve, B::Curve, j)
    throw(DomainError(j,
        "the parallels of an elliptic arc and its neighbour meet at no point of a conic at $(j) m"))
end

"""
$(TYPEDSIGNATURES)

Return the parallel of `contour` at the distance `d` \\[m\\], signed along the outward
normal of the region: positive outward, negative into it. Each curve of the outer path takes
`parallel(curve, d)`. Each curve of a hole's path takes `parallel(curve, -d)`, since that
path runs counterclockwise about the hole.

Each curve ends at its junction with the next, given by `parallel(a, b, d)`. Where the two
parallels meet, both end at that point. Where they diverge, an arc of radius ``|d|`` about
the corner joins the two ends. The result is the offset outline.

# Errors

Throws `DomainError` where the parallel fails, and names the curve and its path. Inward, the
message gives the `dimension` of `contour`. The tests read nominal values against the
tolerance ``\\sqrt{ε}\\,s``, for the `support` ``s`` of `contour`:

- A curve whose own parallel degenerates: ``|d|`` reaches half its `dimension` on the side
  the parallel moves to.
- A piece of a curve's parallel whose signed length is not above the tolerance.
- Curves of the parallel closer to each other than the tolerance, except neighbours and the
  two pieces at one corner.
"""
function parallel(contour::Contour, d::Real)
    T = promote_type(typeof(first(contour.outer)(0)[1]), typeof(d))
    scale = support(contour)
    tolerance = sqrt(eps(float(typeof(nominal(scale))))) * nominal(scale)
    paths = Vector{Curve{T}}[]
    corners = NTuple{3, Int}[]
    for (i, path) in enumerate((contour.outer, contour.holes...))
        d′ = i == 1 ? d : -d
        for (k, curve) in enumerate(path)
            nominal(abs(d′)) < nominal(dimension(d′ < 0 ? curve : reverse(curve))) / 2 ||
                throw(DomainError(d,
                    nominal(d) < 0 ?
                    "the inward parallel at depth $(abs(d)) m has degenerated at curve $k of path $i: the distance across is $(dimension(contour)) m" :
                    "the outward parallel at distance $d m self-intersects at curve $k of path $i"))
        end
        if length(path) == 1
            push!(paths, Curve{T}[parallel(only(path), d′)])
            push!(corners, (i, 1, 1))
            continue
        end
        junctions = [parallel(path[k], path[mod1(k + 1, length(path))], d′)
                     for k in eachindex(path)]
        for k in eachindex(path)
            start, p = junctions[mod1(k - 1, length(path))][2], junctions[k][1]
            nominal(turn(path[k], d′, start, p).length) > tolerance ||
                throw(DomainError(d,
                    nominal(d) < 0 ?
                    "the inward parallel at depth $(abs(d)) m has degenerated at curve $k of path $i: the distance across is $(dimension(contour)) m" :
                    "the outward parallel at distance $d m self-intersects at curve $k of path $i"))
        end
        outline = Curve{T}[]
        for k in eachindex(path)
            p, q = junctions[k]
            push!(outline, parallel(path[k], d′,
                junctions[mod1(k - 1, length(path))][2], p))
            push!(corners, (i, mod1(k - 1, length(path)), k))
            if nominal.(p) != nominal.(q)
                j = path[k](1)
                θa = direction(path[k], 1)
                θb = direction(path[mod1(k + 1, length(path))], 0)
                τ = π - mod(π - (θb - θa), 2π)
                push!(outline, Arc(Placement(j..., atan(p[2] - j[2], p[1] - j[1])),
                    abs(d′), τ))
                push!(corners, (i, k, k))
            end
        end
        push!(paths, outline)
    end
    curves = reduce(vcat, paths)
    for (k, A) in enumerate(curves), (k′, B) in enumerate(curves)

        k < k′ || continue
        first(corners[k]) == first(corners[k′]) &&
            !isdisjoint(corners[k][2:3], corners[k′][2:3]) && continue
        nominal(clearance(A, B)[3]) > tolerance ||
            throw(DomainError(d,
                nominal(d) < 0 ?
                "the inward parallel at depth $(abs(d)) m has degenerated at curve $(corners[k][3]) of path $(corners[k][1]) and curve $(corners[k′][3]) of path $(corners[k′][1]): the distance across is $(dimension(contour)) m" :
                "the outward parallel at distance $d m self-intersects at curve $(corners[k][3]) of path $(corners[k][1]) and curve $(corners[k′][3]) of path $(corners[k′][1])"))
    end
    return Contour(first(paths), paths[2:end])
end

"""
$(TYPEDSIGNATURES)

Return the winding of `curve` about the point `m` \\[m\\]: the signed fraction of a
revolution that the direction from `m` to the curve sweeps along it. For a segment, that is
the angle between its two ends as seen from `m`.

An arc or an elliptic arc sweeps the angle of its chord. When `m` lies in the bulge between
the curve and its chord, it sweeps one revolution more, in the sense of its span. There `m`
lies strictly on the inner side of the tangents at both ends and at the point of the curve
nearest it. It also lies strictly on the side of the chord where the middle of the curve
lies, unless the curve makes a full turn. The points of the curve itself count as outside.

For a closed path, return the sum of the `winding` of its curves about `m`, rounded to an
integer: the number of turns the path makes about `m`.
"""
function winding(curve::Segment, m)
    u, v = curve.p .- m, curve.q .- m
    return atan(u[1] * v[2] - u[2] * v[1], u[1] * v[1] + u[2] * v[2]) / 2π
end

function winding(curve::Union{Arc, EllipticArc}, m)
    P, Q = curve(0), curve(1)
    u, v = P .- m, Q .- m
    σ = sign(curve.span)
    inside = all((0, 1, apsides(curve, m).near)) do t
        ν = direction(curve, t) - σ * π / 2
        point = curve(t)
        nominal((m[1] - point[1]) * cos(ν) + (m[2] - point[2]) * sin(ν)) < 0
    end
    if abs(nominal(curve.span)) < 2π
        chord, middle = Q .- P, curve(1 / 2)
        inside &= nominal((chord[1] * (middle[2] - P[2]) - chord[2] * (middle[1] - P[1])) *
                          (chord[1] * (m[2] - P[2]) - chord[2] * (m[1] - P[1]))) > 0
    end
    return (atan(u[1] * v[2] - u[2] * v[1], u[1] * v[1] + u[2] * v[2]) + 2π * σ * inside) /
           2π
end

function winding(path::AbstractVector{<:Curve}, m)
    return round(Int, nominal(sum(curve -> winding(curve, m), path)))
end

"""
$(TYPEDSIGNATURES)

Return `true` when the point `m` \\[m\\] lies in the region that `contour` bounds. The
`winding` of its outer path about `m` is then nonzero, and that of each hole is zero. The
points of its curves count as outside.
"""
function Base.in(m::Tuple{Real, Real}, contour::Contour)
    return winding(contour.outer, m) != 0 &&
           all(path -> winding(path, m) == 0, contour.holes)
end

"""
$(TYPEDSIGNATURES)

Return the distance across the region that `contour` bounds \\[m\\]: twice the depth at
which its inward parallel first degenerates. The distance is the least of three bounds:

```math
d = \\min\\left(\\min d_{ab},\\; 2 \\min \\rho,\\; 2 \\min e_c\\right).
```

The width ``d_{ab}`` is the `clearance` of two curves that are not adjacent, or that lie on
different paths. It counts where their inward parallels meet. There the chord between its
two points lies in the inward normal cone of the contour at both ends. Its midpoint lies in
the region.

A clearance of zero is a contact of two such curves, where the region has no thickness,
and the distance is then 0. The bound ``2 \\min \\rho`` is the least `dimension` of a curve
with the region on its left. The outer path's curves count as they run, and a hole's
curves reversed.

The depth ``e_c`` is where the `parallel` piece of a face ``c`` with a corner at an end
vanishes between the parallels of its neighbours. Bisection of the signed length of the
piece finds it. A segment of length ``L`` between two segments at the interior angles
``\\theta_1`` and ``\\theta_2`` \\[rad\\] has
``e_c = L / (\\cot(\\theta_1/2) + \\cot(\\theta_2/2))``.
"""
function dimension(contour::Contour{T}) where {T}
    bounds = T[]
    paths = (contour.outer, contour.holes...)
    for (i, path) in enumerate(paths), (j, otherpath) in enumerate(paths)

        j < i && continue
        for (k, a) in enumerate(path), (l, b) in enumerate(otherpath)

            if i == j && (l <= k || l == k + 1 || (k == 1 && l == length(path)))
                continue
            end
            p, q, w, s, t = clearance(a, b)
            iszero(nominal(w)) && return zero(T)
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
                if m ∈ contour
                    push!(bounds, w)
                end
            end
        end
    end
    for curve in contour.outer
        push!(bounds, dimension(curve))
    end
    for path in contour.holes, curve in path

        push!(bounds, dimension(reverse(curve)))
    end
    width = min(support(contour, 0) + support(contour, π),
        support(contour, π / 2) + support(contour, 3π / 2))
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
            lo,
            hi = bisect(
                e -> nominal(turn(c, -side * e, parallel(a, c, -side * e)[2],
                    parallel(c, b, -side * e)[1]).length) > 0,
                zero(limit), limit)
            nominal(hi) == nominal(limit) || push!(bounds, 2 * hi)
        end
    end
    return minimum(bounds)
end

"""
    area(contour)
    area(path)
    area(curve, o)

Return the signed area of the region that `contour` bounds \\[m²\\], positive
counterclockwise: the area of its outer path less the area of each hole. The area of a
closed `path` is Green's integral about its first point ``\\mathbf o``,
``A = \\frac12 \\oint (\\mathbf x-\\mathbf o)\\times d\\mathbf x``, the sum of
`area(curve, o)` over its curves.

`area(curve, o)` is the signed area that the segments from the point `o` \\[m\\] to `curve`
sweep. The formulas below measure every point and center from ``\\mathbf o``. A segment
from ``P`` to ``Q`` sweeps ``(P \\times Q)/2``. An arc of radius ``ρ`` about ``c`` from
``θ_1`` to ``θ_2`` sweeps ``ρ^2(θ_2-θ_1)/2 + c\\times(Q-P)/2``.

An elliptic arc sweeps ``ab\\,span/2+c\\times(Q-P)/2``.
Its normal offset ``d`` adds ``dL-d[p⋅τ̂]/2+d²Δα/2``, where
``L=∫ρ\\,dθ``, ``ρ=hypot(b cos θ,a sin θ)``, and
``Δα=span+[atan((a-b)sin θ cos θ,b cos²θ+a sin²θ)]``.
Here ``p`` is the base point and ``τ̂`` its unit tangent.
"""
function area end

function area(curve::Segment, o)
    P, Q = curve.p .- o, curve.q .- o
    return (P[1] * Q[2] - P[2] * Q[1]) / 2
end

function area(curve::Arc, o)
    P, Q = curve(0) .- o, curve(1) .- o
    X, Y = curve.at.x - o[1], curve.at.y - o[2]
    return (curve.r^2 * curve.span +
            X * (Q[2] - P[2]) - Y * (Q[1] - P[1])) / 2
end

function area(curve::EllipticArc, o)
    θ₁, θ₂ = curve.φ0, curve.φ0 + curve.span
    P = curve.at((curve.a * cos(θ₁), curve.b * sin(θ₁))) .- o
    Q = curve.at((curve.a * cos(θ₂), curve.b * sin(θ₂))) .- o
    X, Y = curve.at.x - o[1], curve.at.y - o[2]
    a, b, d = curve.a, curve.b, curve.offset
    measure = turn(curve)
    Δα = sign(curve.span) * 2π * measure.fraction
    c, s = cos(curve.at.φ), sin(curve.at.φ)
    τ₁ = ((-a * sin(θ₁) * c - b * cos(θ₁) * s),
        (-a * sin(θ₁) * s + b * cos(θ₁) * c)) ./ hypot(b * cos(θ₁), a * sin(θ₁))
    τ₂ = ((-a * sin(θ₂) * c - b * cos(θ₂) * s),
        (-a * sin(θ₂) * s + b * cos(θ₂) * c)) ./ hypot(b * cos(θ₂), a * sin(θ₂))
    return (curve.a * curve.b * curve.span +
            X * (Q[2] - P[2]) - Y * (Q[1] - P[1])) / 2 +
           (d * sign(curve.span) * measure.length -
            d / 2 *
            (Q[1] * τ₂[1] + Q[2] * τ₂[2] -
             P[1] * τ₁[1] - P[2] * τ₁[2]) - d^2 / 2 * Δα)
end

function area(path::AbstractVector{<:Curve})
    o = first(path)(0)
    return sum(curve -> area(curve, o), path)
end

function area(contour::Contour{T}) where {T}
    area(contour.outer) - sum(area, contour.holes; init = zero(T))
end

"""
    moment(path)
    moment(curve, o)

Return the first moment of area \\[m³\\] of the region that the closed `path` bounds, about
the origin, as a coordinate pair: ``\\mathbf M = \\mathbf M_o + A\\,\\mathbf o``. Here
``\\mathbf M_o`` is the sum of `moment(curve, o)` over the curves of `path` about its first
point ``\\mathbf o``, and ``A`` is the area of `path`.

`moment(curve, o)` is the first moment about the point `o` \\[m\\] of the region that the
segments from `o` to `curve` sweep. The formulas below measure every point and center
from ``\\mathbf o``. A segment from ``P`` to ``Q`` contributes ``(P+Q)(P\\times Q)/6``. An
arc adds that chord moment and the circular segment's moment:
``c ρ^2(Δ-\\sin Δ)/2 + ρ^3(\\sin θ_2-\\sin θ_1,\\cos θ_1-\\cos θ_2)/3
-ρ^2\\sin Δ(P+Q-2c)/6``, where ``Δ=θ_2-θ_1``.

For an elliptic arc, let ``L=R(at.φ)diag(a,b)``, ``u_i=(cos θ_i,sin θ_i)``,
and ``B=ab(Δ-sin Δ)/2``. Its Green moment is the chord moment plus
``cB+abL((sin θ_2-sin θ_1,cos θ_1-cos θ_2)/3-sin Δ(u_1+u_2)/6)``.
An offset ``d`` adds the band's moment
``d∫p ds+(d²/2)∫p dα+(d²/2)[(p_y,-p_x)]+(d³/3)[(sin α,-cos α)]``
and the segment moments from ``p_2`` to ``p_2+dn_2`` and from ``p_1+dn_1`` to ``p_1``.
The integrals use ``∫p ds=cL+R(a[G(sin θ)],-b[H(cos θ)])`` and
``∫p dα=cΔα+R(a²b[J(sin θ)],-ab²[K(cos θ)])``, with
``G′(u)=√(b²+(a²-b²)u²)``, ``H′(v)=√(a²+(b²-a²)v²)``,
``J′(u)=1/(b²+(a²-b²)u²)``, and ``K′(v)=1/(a²+(b²-a²)v²)``.
These primitives use inverse hyperbolic and trigonometric functions that depend on
the order of the semi-axes, and polynomials for equal axes.
"""
function moment end

function moment(curve::Segment, o)
    P, Q = curve.p .- o, curve.q .- o
    cross = P[1] * Q[2] - P[2] * Q[1]
    x = (P[1] + Q[1]) * cross / 6
    y = (P[2] + Q[2]) * cross / 6
    return (x, y)
end

function moment(curve::Arc, o)
    P, Q = curve(0) .- o, curve(1) .- o
    X, Y = curve.at.x - o[1], curve.at.y - o[2]
    cross = P[1] * Q[2] - P[2] * Q[1]
    θ₁, θ₂ = curve.at.φ, curve.at.φ + curve.span
    ρ = curve.r
    A = ρ^2 * (curve.span - sin(curve.span)) / 2
    x = ((P[1] + Q[1]) * cross / 6 + X * A +
         ρ^3 * (sin(θ₂) - sin(θ₁)) / 3 -
         ρ^2 * sin(curve.span) * (P[1] + Q[1] - 2X) / 6)
    y = ((P[2] + Q[2]) * cross / 6 + Y * A +
         ρ^3 * (cos(θ₁) - cos(θ₂)) / 3 -
         ρ^2 * sin(curve.span) * (P[2] + Q[2] - 2Y) / 6)
    return (x, y)
end

function moment(curve::EllipticArc, o)
    θ₁, θ₂ = curve.φ0, curve.φ0 + curve.span
    P = curve.at((curve.a * cos(θ₁), curve.b * sin(θ₁))) .- o
    Q = curve.at((curve.a * cos(θ₂), curve.b * sin(θ₂))) .- o
    X, Y = curve.at.x - o[1], curve.at.y - o[2]
    cross = P[1] * Q[2] - P[2] * Q[1]
    A = curve.a * curve.b * (curve.span - sin(curve.span)) / 2
    u = (sin(θ₂) - sin(θ₁)) / 3 -
        sin(curve.span) * (cos(θ₁) + cos(θ₂)) / 6
    v = (cos(θ₁) - cos(θ₂)) / 3 -
        sin(curve.span) * (sin(θ₁) + sin(θ₂)) / 6
    c, s = cos(curve.at.φ), sin(curve.at.φ)
    x = ((P[1] + Q[1]) * cross / 6 + X * A +
         curve.a * curve.b * (curve.a * c * u - curve.b * s * v))
    y = ((P[2] + Q[2]) * cross / 6 + Y * A +
         curve.a * curve.b * (curve.a * s * u + curve.b * c * v))
    a, b, d = curve.a, curve.b, curve.offset
    measure = turn(curve)
    Δα = sign(curve.span) * 2π * measure.fraction
    integrals = map((θ₁, θ₂)) do θ
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
        (a * G, -b * H, a^2 * b * J, -a * b^2 * K)
    end
    I = last(integrals) .- first(integrals)
    x += (d * (X * sign(curve.span) * measure.length + c * I[1] -
               s * I[2]) +
          d^2 / 2 * (-X * Δα + c * I[3] - s * I[4]) +
          d^2 / 2 * (Q[2] - P[2]))
    y += (d * (Y * sign(curve.span) * measure.length + s * I[1] +
               c * I[2]) +
          d^2 / 2 * (-Y * Δα + s * I[3] + c * I[4]) -
          d^2 / 2 * (Q[1] - P[1]))
    n₁ = ((b * cos(θ₁) * c - a * sin(θ₁) * s),
        (b * cos(θ₁) * s + a * sin(θ₁) * c)) ./
         hypot(b * cos(θ₁), a * sin(θ₁))
    n₂ = ((b * cos(θ₂) * c - a * sin(θ₂) * s),
        (b * cos(θ₂) * s + a * sin(θ₂) * c)) ./
         hypot(b * cos(θ₂), a * sin(θ₂))
    x += d^3 / 3 * (n₂[2] - n₁[2])
    y -= d^3 / 3 * (n₂[1] - n₁[1])
    for (U, V) in ((Q, Q .+ d .* n₂), (P .+ d .* n₁, P))
        cross = U[1] * V[2] - U[2] * V[1]
        x += (U[1] + V[1]) * cross / 6
        y += (U[2] + V[2]) * cross / 6
    end
    return (x, y)
end

function moment(path::AbstractVector{<:Curve})
    o = first(path)(0)
    return mapreduce(curve -> moment(curve, o), .+, path) .+ o .* area(path)
end

"""
    centroid(contour)

Return the centroid `(x, y)` of the region that `contour` bounds \\[m\\],
``\\mathbf c = \\mathbf M/A``. The first moment of area ``\\mathbf M`` is that of its outer
path less that of each hole, by `moment`, and ``A`` is its area.
"""
function centroid end

function centroid(contour::Contour{T}) where {T}
    return (moment(contour.outer) .-
            mapreduce(moment, .+, contour.holes; init = (zero(T), zero(T)))) ./
           area(contour)
end

"""
    perimeter(contour)

Return the length of every outer and hole curve of `contour` \\[m\\], the sum of ``|Q-P|``
over its segments and of ``ρ|θ_2-θ_1|`` over its arcs. For an elliptic arc, QuadGK
integrates ``√(a²sin²θ+b²cos²θ)`` over its parametric interval, the incomplete elliptic
integral of the second kind.
A normal offset ``d`` changes the signed length ``L`` to ``L+dΔα``,
where ``Δα=span+[atan((a-b)sin θ cos θ,b cos²θ+a sin²θ)]``.
"""
function perimeter end

function perimeter(contour::Contour{T}) where {T}
    result = zero(T)
    for path in (contour.outer, contour.holes...), curve in path

        result += turn(curve).length
    end
    return result
end

"""
    support(contour)

Return the radius of the smallest circle about the origin that contains the region
`contour` bounds \\[m\\], the distance of the farthest point of its outer curves from the
origin. Holes leave the support unchanged.

On an arc, the farthest point lies at the radius beyond the center when the direction of
the center falls within the arc. Otherwise it lies at an end. Every point of an arc about
the origin lies at its radius. An elliptic arc takes the point of its curve farthest from
the origin, at its normal offset.

`support(contour, φ)` returns the support function ``\\max(x\\cos φ+y\\sin φ)`` \\[m\\]
of the region in the direction `φ` \\[rad\\]. It is the largest projection of the outer
curves, each at the near parameter of its `apsides` along `φ`.

A normal offset ``d`` adds ``d`` to the elliptic arc's interior directional maximum. Its
ends use ``p+dn``. Each curve has the supports of its points, and a contour the largest of
its outer curves' supports.
"""
function support end

function support(contour::Contour)
    return maximum(support, contour.outer)
end

function support(contour::Contour, φ::Real)
    return maximum(support(curve, φ) for curve in contour.outer)
end

"""
$(TYPEDSIGNATURES)

Return the lowest point of `contour`'s outer path, ordered by least `y`, then least `x`.
Coordinates and radii are in \\[m\\] and angles in \\[rad\\]. The returned coordinate
pair is in \\[m\\].
"""
function lowest(contour::Contour)
    return argmin(xy -> (xy[2], xy[1]), [lowest(curve) for curve in contour.outer])
end

"""
$(TYPEDSIGNATURES)

Return the lowest point of `curve` \\[m\\], ordered by least `y`, then least `x`, among its
ends and its point at the near parameter of `apsides(curve, 3π / 2)`.
"""
function lowest(curve::Curve)
    argmin(xy -> (xy[2], xy[1]),
        (curve(0), curve(1), curve(apsides(curve, 3π / 2).near)))
end

"""
    distance(contour)
    distance(curve)

Return the radius of a circle about the origin that the region `contour` bounds does not
enter \\[m\\]. It is zero when the origin lies in the region, by `in`, else the least
distance of the curves of its paths.

`distance(curve)` is the distance \\[m\\] of the nearest point of `curve` from the origin.
A segment or an elliptic arc takes the near point of its `apsides` about the origin. On an
arc about the origin, every point lies at the radius ``ρ``. For another arc, of center
``c``, the distance is ``\\lvert |c|-ρ \\rvert`` when the direction from ``c`` toward the
origin falls within the arc, else that of its nearer end.
"""
function distance end

function distance(curve::Arc)
    iszero(distance(curve.at)) && return curve.r
    φ = atan(-curve.at.y, -curve.at.x)
    if nominal(mod(sign(curve.span) * (φ - curve.at.φ), 2π)) <= abs(nominal(curve.span))
        return abs(distance(curve.at) - curve.r)
    end
    return min(distance(curve(0)), distance(curve(1)))
end

function distance(curve::Union{Segment{T}, EllipticArc{T}}) where {T}
    distance(curve(apsides(curve, (zero(T), zero(T))).near))
end

function distance(contour::Contour{T}) where {T}
    (zero(T), zero(T)) ∈ contour && return zero(T)
    return minimum((minimum(distance, path) for path in contour.holes); init = minimum(distance, contour.outer))
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
    return segment.p .+ t .* (segment.q .- segment.p)
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
