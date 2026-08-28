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

"""
$(TYPEDEF)

Pair one physical region with its resolved shape, retained terminal, and
backend-neutral placement and path declarations.

$(TYPEDFIELDS)
"""
struct PlacedRegion{
        R <: Region,
        S <: AbstractShape,
        P <: NamedTuple,
        H <: Tuple
}
    "Authoritative homogeneous physical region."
    source::R
    "Resolved absolute cross-sectional shape."
    shape::S
    "Retained electrical terminal, or `nothing` for a nonconductive region."
    terminal::Union{Nothing, Symbol}
    "Resolved pose and the ordered placement declarations that produced it."
    placement::P
    "Ordered `(path, radius)` declarations traversed by the region."
    paths::H

    function PlacedRegion(
            source::R,
            shape::S,
            terminal::Union{Nothing, Symbol},
            placement::P,
            paths::H
    ) where {
            R <: Region,
            S <: AbstractShape,
            P <: NamedTuple,
            H <: Tuple
    }
        terminal === nothing || !isempty(String(terminal)) || throw(
            ArgumentError("placed-region terminal cannot be empty")
        )
        keys(placement) == (:pose, :patterns) || throw(ArgumentError(
            "placed-region placement must contain pose and patterns"
        ))
        placement.pose isa Pose2 || throw(ArgumentError(
            "placed-region pose must resolve to Pose2"
        ))
        placement.patterns isa Tuple || throw(ArgumentError(
            "placed-region patterns must be a tuple"
        ))
        all(placement.patterns) do entry
            entry isa NamedTuple && keys(entry) == (:pattern, :member, :pose) &&
                entry.member isa Integer && entry.member > 0 &&
                entry.pose isa Pose2
        end || throw(ArgumentError(
            "placed-region patterns must contain positive member indices and poses"
        ))
        all(paths) do entry
            entry isa NamedTuple && keys(entry) == (:path, :radius) &&
                entry.radius isa Real && isfinite(entry.radius) && entry.radius >= 0
        end || throw(ArgumentError(
            "placed-region paths must contain finite nonnegative radii"
        ))
        return new{R, S, P, H}(source, shape, terminal, placement, paths)
    end
end

function PlacedRegion(source::Region, shape::AbstractShape)
    T = promote_type(eltype(shape), eltype(source.material))
    return PlacedRegion(
        source,
        shape,
        nothing,
        (pose = Pose2(zero(T), zero(T), zero(T)), patterns = ()),
        ()
    )
end

boundary(region::PlacedRegion) = boundary(region.shape)
area(region::PlacedRegion) = area(region.shape)
centroid(region::PlacedRegion) = centroid(region.shape)
support(region::PlacedRegion, φ::Real) = support(region.shape, φ)

function resolve(context::Union{EmptyBoundary, AbstractShape}, region::Region)
    shape = resolve(context, region.primitive)
    placed = PlacedRegion(region, shape)
    return CableGeometry(PlacedRegion[placed], boundary(placed))
end
