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
A positive span runs counterclockwise.

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
    @doc """
    $(TYPEDSIGNATURES)

    Construct an elliptic arc through `validate`, with semi-axes in \\[m\\]
    and parametric angles in \\[rad\\].
    """
    function EllipticArc{T}(at::Placement{T}, a::T, b::T, φ0::T, span::T) where {T <: Real}
        return validate(new{T}(at, a, b, φ0, span))
    end
end

"""
$(TYPEDSIGNATURES)

Construct an elliptic arc from its placement, semi-axes \\[m\\], initial parametric
angle and signed span \\[rad\\], promoted with `Float64`.
"""
function EllipticArc(at::Placement, a::Real, b::Real, φ0::Real, span::Real)
    T = promote_type(Float64, eltype(at), typeof(a), typeof(b), typeof(φ0), typeof(span))
    return EllipticArc{T}(convert(Placement{T}, at), convert(T, a), convert(T, b),
        convert(T, φ0), convert(T, span))
end

"""
$(TYPEDSIGNATURES)

Check that the semi-axes are positive and finite \\[m\\], a finite initial angle and
``0<|span|≤2π`` \\[rad\\]. Return `arc`, or throw `DomainError`.
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
    return arc
end

"""
$(TYPEDSIGNATURES)

Return ``at(a\\cos θ,b\\sin θ)`` \\[m\\], with ``θ=φ0+s\\,span`` \\[rad\\]
for ``s∈[0,1]`` \\[dimensionless\\].
"""
function (arc::EllipticArc)(s::Real)
    θ = arc.φ0 + s * arc.span
    return arc.at((arc.a * cos(θ), arc.b * sin(θ)))
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
    arcs = [[(at = curve.at, a = curve.r, b = curve.r, φ0 = zero(T), span = curve.span)
             for curve in path if curve isa Arc] for path in paths]
    for (i, path) in enumerate(paths)
        append!(arcs[i],
            [(at = curve.at, a = curve.a, b = curve.b,
                 φ0 = curve.φ0, span = curve.span)
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
                    value = projection + h
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

Place the outline by ``P\\mapsto\\mathrm{at}(P)`` for segment ends and
``c\\mapsto\\mathrm{at}\\,c`` for arc placements. Coordinates are in \\[m\\]
and rotations in \\[rad\\].
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
        outline[k] = EllipticArc(at * curve.at, curve.a, curve.b, curve.φ0, curve.span)
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
            P, Q = curve(0), curve(1)
            result += orientation * (curve.a * curve.b * curve.span +
                       curve.at.x * (Q[2] - P[2]) - curve.at.y * (Q[1] - P[1])) / 2
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
            P, Q = curve(0), curve(1)
            θ₁, θ₂ = curve.φ0, curve.φ0 + curve.span
            cross = P[1] * Q[2] - P[2] * Q[1]
            A = curve.a * curve.b * (curve.span - sin(curve.span)) / 2
            u = (sin(θ₂) - sin(θ₁)) / 3 -
                sin(curve.span) * (cos(θ₁) + cos(θ₂)) / 6
            v = (cos(θ₁) - cos(θ₂)) / 3 -
                sin(curve.span) * (sin(θ₁) + sin(θ₂)) / 6
            c, d = cos(curve.at.φ), sin(curve.at.φ)
            x += orientation * ((P[1] + Q[1]) * cross / 6 + curve.at.x * A +
                  curve.a * curve.b * (curve.a * c * u - curve.b * d * v))
            y += orientation * ((P[2] + Q[2]) * cross / 6 + curve.at.y * A +
                  curve.a * curve.b * (curve.a * d * u + curve.b * c * v))
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
            result += abs(first(quadgk(θ -> hypot(curve.a * sin(θ), curve.b * cos(θ)),
                curve.φ0, curve.φ0 + curve.span;
                rtol = eps(float(nominal(one(T))))^(3 / 4),
                norm = x -> abs(nominal(x)))))
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
        result = max(result, hypot(curve.at.x, curve.at.y) + max(curve.a, curve.b))
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
                hypot(curve.a * cos(ψ), curve.b * sin(ψ)))
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
            radius = max(curve.a, curve.b)
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
