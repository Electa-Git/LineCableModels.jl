"Resolve an intrinsic primitive against a geometric context."
function resolve end

resolve(::EmptyBoundary, primitive::Disk) = DiskShape(primitive.r)
resolve(::EmptyBoundary, primitive::Rectangle) =
    RectangleShape(primitive.w, primitive.h)
resolve(::EmptyBoundary, primitive::Ellipse) = EllipseShape(primitive.a, primitive.b)
resolve(::EmptyBoundary, primitive::Sector) =
    SectorShape(primitive.ri, primitive.ro, primitive.φ0, primitive.span)
resolve(::EmptyBoundary, primitive::Annulus) =
    AnnulusShape(primitive.ri, primitive.ro)
resolve(::EmptyBoundary, primitive::Polygon) = PolygonShape(primitive.points)

function resolve(context::DiskShape, primitive::Shell)
    values = map(float, promote(context.r, primitive.t))
    inner, layer = values
    return AnnulusShape(inner, inner + layer)
end

function resolve(context::AnnulusShape, primitive::Shell)
    values = map(float, promote(context.ro, primitive.t))
    inner, layer = values
    return AnnulusShape(inner, inner + layer)
end

function resolve(
        context::Union{DiskShape, AnnulusShape},
        primitive::Annulus
)
    inner = r_ex(context)
    isapprox(primitive.ri, inner) || throw(DomainError(
        primitive.ri,
        "annulus inner radius must equal the current outer radius $inner"
    ))
    return AnnulusShape(primitive.ri, primitive.ro)
end

function PlacedShape(shape::S, at::P) where {S <: AbstractShape, P <: Pose2}
    T = promote_type(eltype(shape), eltype(at))
    return PlacedShape{T, S, P}(shape, at)
end

PlacedShape(shape::PlacedShape, at::Pose2) = PlacedShape(shape.shape, at * shape.at)

resolve(at::Pose2, primitive::AbstractPrimitive) =
    PlacedShape(resolve(EmptyBoundary(), primitive), at)

area(shape::PlacedShape) = area(shape.shape)
function centroid(shape::PlacedShape)
    x, y = centroid(shape.shape)
    return (
        shape.at.x + cos(shape.at.φ) * x - sin(shape.at.φ) * y,
        shape.at.y + sin(shape.at.φ) * x + cos(shape.at.φ) * y
    )
end
support(shape::PlacedShape, φ::Real) =
    shape.at.x * cos(φ) + shape.at.y * sin(φ) +
    support(shape.shape, φ - shape.at.φ)
support(shape::PlacedShape) = hypot(shape.at.x, shape.at.y) + support(shape.shape)
boundary(shape::PlacedShape) = PlacedShape(boundary(shape.shape), shape.at)

r_in(shape::PlacedShape) = r_in(shape.shape)
r_ex(shape::PlacedShape) = r_ex(shape.shape)
thickness(shape::PlacedShape) = thickness(shape.shape)

function resolve(context::PlacedShape{<:Real, <:Union{DiskShape, AnnulusShape}}, primitive::Shell)
    return PlacedShape(resolve(context.shape, primitive), context.at)
end

function resolve(
        context::PlacedShape{<:Real, <:Union{DiskShape, AnnulusShape}},
        primitive::Annulus
)
    return PlacedShape(resolve(context.shape, primitive), context.at)
end
