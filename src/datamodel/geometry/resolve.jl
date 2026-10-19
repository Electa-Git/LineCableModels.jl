resolve(::EmptyBoundary, primitive::AbstractPrimitive) = primitive

function resolve(context::Disk, definition::Shell)
    values = map(float, promote(context.r, definition.t))
    inner, layer = values
    return Annulus(inner, inner + layer, context.at)
end

function resolve(context::Annulus, definition::Shell)
    values = map(float, promote(context.ro, definition.t))
    inner, layer = values
    return Annulus(inner, inner + layer, context.at)
end

# A shell around a rectangle is the larger rectangle at the same pose, with sharp corners.
# The layer is t thick on the sides and t√2 at the corners.
function resolve(inner::Rectangle, layer::Shell)
    w, h, t = map(float, promote(inner.w, inner.h, layer.t))
    return ShellShape(Rectangle(w, h, inner.at), Rectangle(w + 2t, h + 2t, inner.at))
end

# A shell around a polygon moves every edge outward by t along its normal. Each new vertex
# is the intersection of the two adjacent moved edges, a sharp corner as around a
# rectangle. Where adjacent edges lie on one line, their vertex moves along their common
# normal.
function resolve(inner::Polygon, layer::Shell)
    T = float(promote_type(eltype(inner), typeof(layer.t)))
    t = convert(T, layer.t)
    points = [(convert(T, point[1]), convert(T, point[2])) for point in inner.points]
    vertex_count = length(points)
    orientation = sign(nominal(signed_polygon_area(points)))
    edges = map(eachindex(points)) do index
        first_point = points[index]
        last_point = points[mod1(index + 1, vertex_count)]
        dx, dy = last_point[1] - first_point[1], last_point[2] - first_point[2]
        normal = (orientation * dy / hypot(dx, dy), -orientation * dx / hypot(dx, dy))
        (normal, normal[1] * first_point[1] + normal[2] * first_point[2] + t)
    end
    outer = map(eachindex(points)) do index
        (a, c), (b, d) = edges[mod1(index - 1, vertex_count)], edges[index]
        determinant = a[1] * b[2] - a[2] * b[1]
        abs(nominal(determinant)) <= sqrt(eps(Float64)) &&
            return (points[index][1] + t * b[1], points[index][2] + t * b[2])
        return ((c * b[2] - d * a[2]) / determinant, (a[1] * d - b[1] * c) / determinant)
    end
    at = convert(Pose2{T}, inner.at)
    inside = Polygon{T, Pose2{T}}(points, at)
    return ShellShape(inside,
        validate(Polygon{T, Pose2{T}}(outer, at), ShellShape, inside, layer.t))
end

# The outer polygon of a shell must be simple: no moved edge reverses its direction, and
# no two edges cross. Otherwise the shell thickness is too large for the outline.
function validate(outer::Polygon, ::Type{ShellShape}, inner::Polygon, t)
    vertex_count = length(outer.points)
    reversed = any(1:vertex_count) do index
        next = mod1(index + 1, vertex_count)
        before = inner.points[next] .- inner.points[index]
        after = outer.points[next] .- outer.points[index]
        nominal(before[1] * after[1] + before[2] * after[2]) <= 0
    end
    crossed = any(Iterators.product(1:vertex_count, 1:vertex_count)) do (i, j)
        j > i + 1 && !(i == 1 && j == vertex_count) || return false
        a, b = outer.points[i], outer.points[mod1(i + 1, vertex_count)]
        c, d = outer.points[j], outer.points[mod1(j + 1, vertex_count)]
        side_c = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
        side_d = (b[1] - a[1]) * (d[2] - a[2]) - (b[2] - a[2]) * (d[1] - a[1])
        side_a = (d[1] - c[1]) * (a[2] - c[2]) - (d[2] - c[2]) * (a[1] - c[1])
        side_b = (d[1] - c[1]) * (b[2] - c[2]) - (d[2] - c[2]) * (b[1] - c[1])
        nominal(side_c * side_d) < 0 && nominal(side_a * side_b) < 0
    end
    (reversed || crossed) && throw(ArgumentError(
        "a shell thickness of $t m is too large for this Polygon outline: its outer " *
        "outline would not be simple"
    ))
    return outer
end

function resolve(
        context::Union{Disk, Annulus},
        definition::Annulus
)
    inner = r_ex(context)
    tolerance = sqrt(eps(float(inner))) * max(one(inner), inner)
    definition.ri + tolerance >= inner || throw(DomainError(
        definition.ri,
        "annulus inner radius must not overlap the current outer radius $inner"
    ))
    return Annulus(definition.ri, definition.ro, context.at)
end

"""
Compose `at` with the existing absolute primitive pose.
"""
resolve(at::Pose2, primitive::AbstractPrimitive) = _with_pose(primitive, at * primitive.at)

function resolve(at::Pose2, primitive::DifferenceShape)
    return DifferenceShape(
        resolve(at, primitive.outer),
        map(hole -> resolve(at, hole), primitive.holes)
    )
end
