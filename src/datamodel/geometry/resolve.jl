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

function PlacedShape(shape::S, at::P) where {S <: AbstractShape, P <: Pose2}
    T = promote_type(eltype(shape), eltype(at))
    return PlacedShape{T, S, P}(shape, at)
end

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
boundary(shape::PlacedShape) = PlacedShape(boundary(shape.shape), shape.at)

const _PlacedRadialShape = PlacedShape{
    <:Real,
    <:Union{DiskShape, AnnulusShape, SectorShape}
}
r_in(shape::_PlacedRadialShape) = r_in(shape.shape)
r_ex(shape::_PlacedRadialShape) = r_ex(shape.shape)
thickness(shape::_PlacedRadialShape) = thickness(shape.shape)

function resolve(context::PlacedShape{<:Real, <:Union{DiskShape, AnnulusShape}}, primitive::Shell)
    return PlacedShape(resolve(context.shape, primitive), context.at)
end
