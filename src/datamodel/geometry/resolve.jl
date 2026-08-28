"Resolve one primitive definition against a geometric context."
function resolve end

resolve(::EmptyBoundary, definition::DiskDefinition) = Disk(definition.r)
resolve(::EmptyBoundary, definition::RectangleDefinition) =
    Rectangle(definition.w, definition.h)
resolve(::EmptyBoundary, definition::EllipseDefinition) =
    Ellipse(definition.a, definition.b)
resolve(::EmptyBoundary, definition::SectorDefinition) =
    Sector(definition.ri, definition.ro, definition.φ0, definition.span)
resolve(::EmptyBoundary, definition::AnnulusDefinition) =
    Annulus(definition.ri, definition.ro)
resolve(::EmptyBoundary, definition::PolygonDefinition) = Polygon(definition.points)

function resolve(context::Disk, definition::ShellDefinition)
    values = map(float, promote(context.r, definition.t))
    inner, layer = values
    return Annulus(inner, inner + layer, context.at)
end

function resolve(context::Annulus, definition::ShellDefinition)
    values = map(float, promote(context.ro, definition.t))
    inner, layer = values
    return Annulus(inner, inner + layer, context.at)
end

function resolve(context::Sector, definition::ShellDefinition)
    values = map(float, promote(context.ro, definition.t))
    inner, layer = values
    return Sector(
        inner,
        inner + layer,
        context.φ0,
        context.span,
        context.at
    )
end

function resolve(
        context::Union{Disk, Annulus},
        definition::AnnulusDefinition
)
    inner = r_ex(context)
    tolerance = sqrt(eps(float(inner))) * max(one(inner), inner)
    definition.ri + tolerance >= inner || throw(DomainError(
        definition.ri,
        "annulus inner radius must not overlap the current outer radius $inner"
    ))
    return Annulus(definition.ri, definition.ro, context.at)
end

function resolve(
        context::Union{Disk, Annulus, Sector},
        definition::SectorDefinition
)
    inner = r_ex(context)
    tolerance = sqrt(eps(float(inner))) * max(one(inner), inner)
    definition.ri + tolerance >= inner || throw(DomainError(
        definition.ri,
        "sector inner radius must not overlap the current outer radius $inner"
    ))
    return Sector(
        definition.ri,
        definition.ro,
        definition.φ0,
        definition.span,
        context.at
    )
end

resolve(at::Pose2, definition::AbstractPrimitiveDefinition) =
    resolve(at, resolve(EmptyBoundary(), definition))

"Compose `at` with the existing absolute primitive pose."
resolve(at::Pose2, primitive::AbstractPrimitive) =
    _with_pose(primitive, at * primitive.at)

function resolve(at::Pose2, primitive::_DifferencePrimitive)
    return _DifferencePrimitive(
        resolve(at, primitive.outer),
        map(hole -> resolve(at, hole), primitive.holes)
    )
end
