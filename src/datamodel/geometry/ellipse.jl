"""
$(TYPEDEF)

Represent the exact outward parallel boundary of an ellipse.

For ellipse support ``h(\\phi)`` and normal offset ``t``, the parallel boundary
has support

```math
h_t(\\phi) = h(\\phi) + t.
```

`EllipseOffset` is resolved geometry produced by applying [`Shell`](@ref); it
is not a modelling declaration.

$(TYPEDFIELDS)
"""
struct EllipseOffset{T <: Real, P <: Pose2{T}} <: AbstractShape{T}
    "Original ellipse semi-axis along the local x-axis \\[m\\]."
    a::T
    "Original ellipse semi-axis along the local y-axis \\[m\\]."
    b::T
    "Outward normal offset \\[m\\]."
    t::T
    "Pose in the completed cable coordinate system."
    at::P

    function EllipseOffset{T, P}(a::T, b::T, t::T, at::P) where {
            T <: Real, P <: Pose2{T}
    }
        isfinite(a) && a > zero(a) || throw(DomainError(
            a, "ellipse-offset semi-axis a must be positive and finite"
        ))
        isfinite(b) && b > zero(b) || throw(DomainError(
            b, "ellipse-offset semi-axis b must be positive and finite"
        ))
        isfinite(t) && t > zero(t) || throw(DomainError(
            t, "ellipse offset must be positive and finite"
        ))
        return new{T, P}(a, b, t, at)
    end
end

function EllipseOffset(ellipse::Ellipse, t::Real)
    T = promote_type(eltype(ellipse), typeof(t))
    return EllipseOffset{T, Pose2{T}}(
        convert(T, ellipse.a),
        convert(T, ellipse.b),
        convert(T, t),
        convert(Pose2{T}, ellipse.at)
    )
end

boundary(shape::EllipseOffset) = shape
area(shape::EllipseOffset) =
    π * shape.a * shape.b + _ellipse_perimeter(shape.a, shape.b) * shape.t +
    π * shape.t^2
perimeter(shape::EllipseOffset) =
    _ellipse_perimeter(shape.a, shape.b) + 2π * shape.t
centroid(shape::EllipseOffset) = (shape.at.x, shape.at.y)

function support(shape::EllipseOffset, angle::Real)
    local_angle = angle - shape.at.φ
    return shape.at.x * cos(angle) + shape.at.y * sin(angle) +
           hypot(shape.a * cos(local_angle), shape.b * sin(local_angle)) + shape.t
end

support(shape::EllipseOffset) =
    hypot(shape.at.x, shape.at.y) + max(shape.a, shape.b) + shape.t

function resolve(at::Pose2, shape::EllipseOffset)
    return EllipseOffset(
        Ellipse(shape.a, shape.b, at * shape.at),
        shape.t
    )
end

function resolve(inner::Ellipse, layer::Shell)
    return ShellShape(inner, EllipseOffset(inner, layer.t))
end

function resolve(inner::EllipseOffset, layer::Shell)
    outer = EllipseOffset(
        Ellipse(inner.a, inner.b, inner.at),
        inner.t + layer.t
    )
    return ShellShape(inner, outer)
end

function tessellate(shape::Ellipse; points_per_arc::Integer = 128)
    points_per_arc >= 8 || throw(ArgumentError(
        "an ellipse requires at least eight boundary points"
    ))
    angles = range(zero(shape.a), oftype(shape.a, 2π); length = Int(points_per_arc))
    local_points = [
        (shape.a * cos(angle), shape.b * sin(angle)) for angle in angles
    ]
    return [_transform_point(point, shape.at) for point in local_points]
end

function tessellate(shape::EllipseOffset; points_per_arc::Integer = 128)
    points_per_arc >= 8 || throw(ArgumentError(
        "an ellipse offset requires at least eight boundary points"
    ))
    angles = range(zero(shape.a), oftype(shape.a, 2π); length = Int(points_per_arc))
    local_points = map(angles) do angle
        cosine = cos(angle)
        sine = sin(angle)
        normalizer = hypot(shape.b * cosine, shape.a * sine)
        (
            shape.a * cosine + shape.t * shape.b * cosine / normalizer,
            shape.b * sine + shape.t * shape.a * sine / normalizer
        )
    end
    return [_transform_point(point, shape.at) for point in local_points]
end
