"Resolve one primitive or contextual layer against a geometric context."
function resolve end

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

"Compose `at` with the existing absolute primitive pose."
resolve(at::Pose2, primitive::AbstractPrimitive) = _with_pose(primitive, at * primitive.at)

function resolve(at::Pose2, primitive::DifferenceShape)
    return DifferenceShape(
        resolve(at, primitive.outer),
        map(hole -> resolve(at, hole), primitive.holes)
    )
end
