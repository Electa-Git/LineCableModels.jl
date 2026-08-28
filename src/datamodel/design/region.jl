"""
$(TYPEDEF)

Represent one homogeneous physical domain with intrinsic geometry.

$(TYPEDFIELDS)
"""
struct Region{P <: AbstractPrimitive, M} <: AbstractCablePart
    "Physical identity within its containing cable object."
    tag::Symbol
    "Intrinsic cross-sectional geometry."
    primitive::P
    "Constitutive material data."
    material::M

    function Region(tag::Symbol, primitive::P, material::M) where {
            P <: AbstractPrimitive, M
    }
        isempty(String(tag)) && throw(ArgumentError("region tag cannot be empty"))
        material isa Material ||
            throw(ArgumentError("region material must resolve to Material"))
        return new{P, M}(tag, primitive, material)
    end
end

"Pair one physical region with resolved geometry and derived placement data."
struct ResolvedRegion{R <: Region, S <: AbstractShape, T <: Real}
    region::R
    shape::S
    terminal::Union{Nothing, Symbol}
    overlength::T
    turns::T
    pattern_depth::Int
    path_depth::Int
end

function ResolvedRegion(
        region::R,
        shape::S,
        terminal::Union{Nothing, Symbol},
        overlength::Real,
        turns::Real,
        pattern_depth::Int,
        path_depth::Int
) where {R <: Region, S <: AbstractShape}
    T = promote_type(
        eltype(shape),
        eltype(region.material),
        typeof(overlength),
        typeof(turns)
    )
    return ResolvedRegion{R, S, T}(
        region,
        shape,
        terminal,
        convert(T, overlength),
        convert(T, turns),
        pattern_depth,
        path_depth
    )
end

function ResolvedRegion(region::Region, shape::AbstractShape)
    T = promote_type(eltype(shape), eltype(region.material))
    return ResolvedRegion(
        region,
        shape,
        nothing,
        convert(T, one(T)),
        convert(T, zero(T)),
        0,
        0
    )
end

boundary(region::ResolvedRegion) = boundary(region.shape)
area(region::ResolvedRegion) = area(region.shape)
centroid(region::ResolvedRegion) = centroid(region.shape)
support(region::ResolvedRegion, φ::Real) = support(region.shape, φ)

function resolve(context::Union{EmptyBoundary, AbstractShape}, region::Region)
    shape = resolve(context, region.primitive)
    resolved = ResolvedRegion(region, shape)
    return ResolvedPart(ResolvedRegion[resolved], boundary(resolved), Symbol[])
end
