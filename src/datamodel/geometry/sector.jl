"""
$(TYPEDEF)

Define one material-neutral cable-sector primitive with optional fillets.

The symmetry axis is local `+x`. Translation and rotation belong to
[`Pose2`](@ref), while material and terminal identity belong to `Region` and
its containing physical tree.

$(TYPEDFIELDS)
"""
struct Sector{T <: Real} <: AbstractPrimitive{T}
    "Angular opening \\[rad\\]."
    span::T
    "Radial coordinate of the innermost midpoint \\[m\\]."
    r_base::T
    "Radius of the circular outer back \\[m\\]."
    r_back::T
    "Common base and side fillet radius \\[m\\]."
    fillet::T

    function Sector{T}(
            span::T,
            r_base::T,
            r_back::T,
            fillet::T
    ) where {T <: Real}
        all(isfinite, (span, r_base, r_back, fillet)) || throw(ArgumentError(
            "sector dimensions must be finite"
        ))
        zero(span) < span < oftype(span, π) || throw(DomainError(
            span, "sector span must lie in (0, π)"
        ))
        zero(r_base) <= r_base < r_back || throw(DomainError(
            r_base, "sector base radius must lie in [0, r_back)"
        ))
        zero(fillet) <= fillet < r_back || throw(DomainError(
            fillet, "sector fillet must lie in [0, r_back)"
        ))
        value = new{T}(span, r_base, r_back, fillet)
        _sector_contacts(value)
        return value
    end
end

function Sector(
        span::Real,
        r_base::Real,
        r_back::Real,
        fillet::Real = 0
)
    values = map(float, promote(span, r_base, r_back, fillet))
    return Sector{typeof(first(values))}(values...)
end

function Base.convert(
        ::Type{<:AbstractPrimitive{T}},
        value::Sector
) where {T <: Real}
    return Sector{T}(
        convert(T, value.span),
        convert(T, value.r_base),
        convert(T, value.r_back),
        convert(T, value.fillet)
    )
end

"""
$(TYPEDEF)

Store the exact contact geometry of one placed [`Sector`](@ref).

`contacts` contains passive derived arc and segment data. It is rebuilt from
the authoritative primitive and is never serialized.
"""
struct SectorShape{
        T <: Real,
        P <: Sector{T},
        C <: NamedTuple,
        A <: Pose2{T}
} <: AbstractShape{T}
    "Authoritative intrinsic cable-sector primitive."
    primitive::P
    "Exact derived arc contacts and straight boundary segments."
    contacts::C
    "Pose in the completed cable coordinate system."
    at::A
end

"""
$(TYPEDEF)

Represent the exact material domain between two resolved boundaries.
"""
struct ShellShape{
        T <: Real,
        I <: AbstractShape{T},
        O <: AbstractShape{T}
} <: AbstractShape{T}
    "Resolved inner boundary."
    inner::I
    "Resolved parallel outer boundary."
    outer::O

    function ShellShape(inner::I, outer::O) where {
            T <: Real, I <: AbstractShape{T}, O <: AbstractShape{T}
    }
        shell_area = area(outer) - area(inner)
        shell_area > zero(shell_area) || throw(DomainError(
            shell_area, "a conformal shell must have positive area"
        ))
        return new{T, I, O}(inner, outer)
    end
end

_geometry_scalar(value) = float(nominal(value))

function _geometry_tolerance(value)
    scale = abs(_geometry_scalar(value))
    return 64 * eps(iszero(scale) ? one(scale) : scale)
end

function _positive_arc(start, stop)
    period = oftype(start + stop, 2π)
    while stop < start
        stop += period
    end
    return (start = start, stop = stop)
end

function _arc(center, radius, first_point, last_point)
    start = atan(first_point[2] - center[2], first_point[1] - center[1])
    stop = atan(last_point[2] - center[2], last_point[1] - center[1])
    interval = _positive_arc(start, stop)
    return (; center, radius, interval.start, interval.stop)
end

function _sector_contacts(primitive::Sector)
    α = primitive.span
    b = primitive.r_base
    radius = primitive.r_back
    fillet = primitive.fillet
    half_complement = π / 2 - α / 2
    cosine = cos(half_complement)
    sine = sin(half_complement)
    tangent = tan(half_complement)

    base_offset = b - fillet * (inv(cosine) - one(cosine))
    line_offset = b + fillet
    qa = one(tangent) + tangent^2
    qb = 2 * line_offset * tangent
    qc = line_offset^2 - (radius - fillet)^2
    discriminant = qb^2 - 4 * qa * qc
    tolerance = _geometry_tolerance(max(abs(qb^2), abs(4 * qa * qc)))
    discriminant_value = _geometry_scalar(discriminant)
    discriminant_value < -tolerance && throw(DomainError(
        discriminant,
        "sector side fillet has no valid tangent contact"
    ))
    discriminant_value < 0 && (discriminant = zero(discriminant))

    old_x = (-qb + sqrt(discriminant)) / (2 * qa)
    old_y = old_x * tangent + line_offset
    side_distance = hypot(old_x, old_y)
    distance_tolerance = _geometry_tolerance(radius)
    abs(_geometry_scalar(side_distance - (radius - fillet))) <=
        distance_tolerance || throw(DomainError(
        side_distance,
        "sector side contact is inconsistent with its outer back"
    ))

    base_center = (b + fillet, zero(b))
    base_x = base_offset + fillet * sine * tangent
    base_lower = (base_x, -fillet * sine)
    base_upper = (base_x, fillet * sine)
    lower_center = (old_y, -old_x)
    upper_center = (old_y, old_x)
    side_lower = (old_y - fillet * cosine, -(old_x + fillet * sine))
    side_upper = (side_lower[1], -side_lower[2])
    back_lower = (old_y * radius / side_distance, -old_x * radius / side_distance)
    back_upper = (back_lower[1], -back_lower[2])

    lower_direction = (
        side_lower[1] - base_lower[1],
        side_lower[2] - base_lower[2]
    )
    lower_direction[1] * cosine - lower_direction[2] * sine >=
        -distance_tolerance || throw(DomainError(
        lower_direction, "sector straight side has reversed contact order"
    ))

    arcs = (
        base = _arc(base_center, fillet, base_upper, base_lower),
        lower = _arc(lower_center, fillet, side_lower, back_lower),
        back = _arc((zero(radius), zero(radius)), radius, back_lower, back_upper),
        upper = _arc(upper_center, fillet, back_upper, side_upper)
    )
    all(arc -> begin
        iszero(arc.radius) && return true
        extent = _geometry_scalar(arc.stop - arc.start)
        -distance_tolerance <= extent <= π + distance_tolerance
    end, values(arcs)) || throw(DomainError(
        arcs, "sector arc contacts do not form a convex boundary"
    ))
    segments = (
        lower = (base_lower, side_lower),
        upper = (side_upper, base_upper)
    )
    points = (;
        base_upper,
        base_lower,
        side_lower,
        back_lower,
        back_upper,
        side_upper
    )
    return (; points, arcs, segments)
end

function SectorShape(primitive::Sector, at::Pose2 = _origin(eltype(primitive)))
    T = promote_type(eltype(primitive), eltype(at))
    converted = convert(AbstractPrimitive{T}, primitive)
    pose = convert(Pose2{T}, at)
    contacts = _sector_contacts(converted)
    return SectorShape{T, typeof(converted), typeof(contacts), typeof(pose)}(
        converted, contacts, pose
    )
end

function _line_moments(first_point, last_point)
    x0, y0 = first_point
    x1, y1 = last_point
    return (
        area = (x0 * y1 - x1 * y0) / 2,
        xmoment = (y1 - y0) * (x0^2 + x0 * x1 + x1^2) / 6,
        ymoment = -(x1 - x0) * (y0^2 + y0 * y1 + y1^2) / 6
    )
end

function _arc_moments(arc)
    cx, cy = arc.center
    radius = arc.radius
    start = arc.start
    stop = arc.stop
    span = stop - start
    cosine = sin(stop) - sin(start)
    sine = -cos(stop) + cos(start)
    cosine_squared = span / 2 + (sin(2 * stop) - sin(2 * start)) / 4
    sine_squared = span / 2 - (sin(2 * stop) - sin(2 * start)) / 4
    cosine_cubed =
        sin(stop) - sin(stop)^3 / 3 - sin(start) + sin(start)^3 / 3
    sine_cubed =
        -cos(stop) + cos(stop)^3 / 3 + cos(start) - cos(start)^3 / 3
    return (
        area = (
            radius * cx * cosine + radius * cy * sine + radius^2 * span
        ) / 2,
        xmoment = radius * (
            cx^2 * cosine + 2 * cx * radius * cosine_squared +
            radius^2 * cosine_cubed
        ) / 2,
        ymoment = radius * (
            cy^2 * sine + 2 * cy * radius * sine_squared +
            radius^2 * sine_cubed
        ) / 2
    )
end

function _sector_moments(shape::SectorShape)
    moments = (
        _arc_moments(shape.contacts.arcs.base),
        _line_moments(shape.contacts.segments.lower...),
        _arc_moments(shape.contacts.arcs.lower),
        _arc_moments(shape.contacts.arcs.back),
        _arc_moments(shape.contacts.arcs.upper),
        _line_moments(shape.contacts.segments.upper...)
    )
    return (
        area = sum(getproperty.(moments, :area)),
        xmoment = sum(getproperty.(moments, :xmoment)),
        ymoment = sum(getproperty.(moments, :ymoment))
    )
end

area(shape::SectorShape) = _sector_moments(shape).area
function perimeter(shape::SectorShape)
    arc_length = sum(values(shape.contacts.arcs)) do arc
        arc.radius * (arc.stop - arc.start)
    end
    line_length = sum(values(shape.contacts.segments)) do segment
        hypot(
            segment[2][1] - segment[1][1],
            segment[2][2] - segment[1][2]
        )
    end
    return arc_length + line_length
end

function _transform_point(point, at::Pose2)
    cosine = cos(at.φ)
    sine = sin(at.φ)
    return (
        at.x + cosine * point[1] - sine * point[2],
        at.y + sine * point[1] + cosine * point[2]
    )
end

function centroid(shape::SectorShape)
    moments = _sector_moments(shape)
    local_centroid = (
        moments.xmoment / moments.area,
        moments.ymoment / moments.area
    )
    return _transform_point(local_centroid, shape.at)
end

function _angle_in_arc(angle, arc)
    period = oftype(angle, 2π)
    candidate = angle + ceil((arc.start - angle) / period) * period
    tolerance = _geometry_tolerance(arc.stop - arc.start)
    return _geometry_scalar(candidate - arc.stop) <= tolerance
end

function _local_support(shape::SectorShape, angle::Real)
    cosine = cos(angle)
    sine = sin(angle)
    projection(point) = point[1] * cosine + point[2] * sine
    candidates = [projection(point) for point in values(shape.contacts.points)]
    for arc in values(shape.contacts.arcs)
        _angle_in_arc(angle, arc) || continue
        push!(candidates, projection(arc.center) + arc.radius)
    end
    return maximum(candidates)
end

function support(shape::SectorShape, angle::Real)
    return shape.at.x * cos(angle) + shape.at.y * sin(angle) +
           _local_support(shape, angle - shape.at.φ)
end

support(shape::SectorShape) = hypot(shape.at.x, shape.at.y) + shape.primitive.r_back
boundary(shape::SectorShape) = shape
r_in(shape::SectorShape) = shape.primitive.r_base
r_ex(shape::SectorShape) = shape.primitive.r_back
thickness(shape::SectorShape) = r_ex(shape) - r_in(shape)

resolve(::EmptyBoundary, primitive::Sector) = SectorShape(primitive)
resolve(at::Pose2, primitive::Sector) = SectorShape(primitive, at)
function resolve(at::Pose2, shape::SectorShape)
    return SectorShape(shape.primitive, at * shape.at)
end

function resolve(inner::SectorShape, layer::Shell)
    primitive = inner.primitive
    outer = Sector(
        primitive.span,
        primitive.r_base - layer.t,
        primitive.r_back + layer.t,
        primitive.fillet + layer.t
    )
    return ShellShape(inner, SectorShape(outer, inner.at))
end

boundary(shape::ShellShape) = shape.outer
area(shape::ShellShape) = area(shape.outer) - area(shape.inner)
perimeter(shape::ShellShape) = perimeter(shape.inner) + perimeter(shape.outer)
support(shape::ShellShape, angle::Real) = support(shape.outer, angle)
support(shape::ShellShape) = support(shape.outer)
function centroid(shape::ShellShape)
    inner_area = area(shape.inner)
    outer_area = area(shape.outer)
    shell_area = outer_area - inner_area
    inner_centroid = centroid(shape.inner)
    outer_centroid = centroid(shape.outer)
    return (
        (outer_area * outer_centroid[1] - inner_area * inner_centroid[1]) /
        shell_area,
        (outer_area * outer_centroid[2] - inner_area * inner_centroid[2]) /
        shell_area
    )
end

function resolve(at::Pose2, shape::ShellShape)
    return ShellShape(resolve(at, shape.inner), resolve(at, shape.outer))
end

function _arc_points(arc, count::Int)
    iszero(arc.radius) && return [arc.center]
    return [
        (
            arc.center[1] + arc.radius * cos(angle),
            arc.center[2] + arc.radius * sin(angle)
        )
        for angle in range(arc.start, arc.stop; length = count)
    ]
end

"""
    tessellate(shape::SectorShape; points_per_arc=32)

Approximate the exact boundary with ordinary coordinate tuples for rendering
or meshing. Exact geometric properties never use these points.
"""
function tessellate(shape::SectorShape; points_per_arc::Integer = 32)
    points_per_arc >= 2 || throw(ArgumentError(
        "points_per_arc must be at least two"
    ))
    contacts = shape.contacts
    points = collect(_arc_points(contacts.arcs.base, Int(points_per_arc)))
    push!(points, contacts.points.side_lower)
    append!(points, _arc_points(contacts.arcs.lower, Int(points_per_arc))[2:end])
    append!(points, _arc_points(contacts.arcs.back, Int(points_per_arc))[2:end])
    append!(points, _arc_points(contacts.arcs.upper, Int(points_per_arc))[2:end])
    return [_transform_point(point, shape.at) for point in points]
end

function tessellate(shape::ShellShape; points_per_arc::Integer = 32)
    return (
        outer = tessellate(shape.outer; points_per_arc),
        inner = tessellate(shape.inner; points_per_arc)
    )
end
