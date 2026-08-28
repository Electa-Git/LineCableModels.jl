"""
$(TYPEDEF)

Represent one homogeneous physical domain with intrinsic geometry.

$(TYPEDFIELDS)
"""
struct Region{P <: AbstractPrimitiveDefinition, M} <: AbstractCablePart
    "Physical identity within its containing cable object."
    tag::Symbol
    "Intrinsic cross-sectional geometry."
    primitive::P
    "Constitutive material data."
    material::M

    function Region(tag::Symbol, primitive::P, material::M) where {
            P <: AbstractPrimitiveDefinition, M
    }
        isempty(String(tag)) && throw(ArgumentError("region tag cannot be empty"))
        material isa Material ||
            throw(ArgumentError("region material must resolve to Material"))
        return new{P, M}(tag, primitive, material)
    end
end

"""
$(TYPEDEF)

Pair one physical region with its resolved primitive, retained terminal, and
backend-neutral placement and path declarations.

$(TYPEDFIELDS)
"""
struct PlacedRegion{
        R <: Region,
        S <: AbstractPrimitive,
        P <: NamedTuple,
        H <: Tuple
}
    "Authoritative homogeneous physical region."
    source::R
    "Resolved cross-sectional primitive with its absolute pose."
    primitive::S
    "Retained electrical terminal, or `nothing` for a nonconductive region."
    terminal::Union{Nothing, Symbol}
    "Ordered placement declarations that produced the resolved primitive pose."
    placement::P
    "Ordered `(path, radius)` declarations traversed by the region."
    paths::H

    function PlacedRegion(
            source::R,
            primitive::S,
            terminal::Union{Nothing, Symbol},
            placement::P,
            paths::H
    ) where {
            R <: Region,
            S <: AbstractPrimitive,
            P <: NamedTuple,
            H <: Tuple
    }
        terminal === nothing || !isempty(String(terminal)) || throw(
            ArgumentError("placed-region terminal cannot be empty")
        )
        keys(placement) == (:patterns,) || throw(ArgumentError(
            "placed-region placement must contain patterns"
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
        return new{R, S, P, H}(source, primitive, terminal, placement, paths)
    end
end

Base.:(==)(left::PlacedRegion, right::PlacedRegion) =
    left.source == right.source && left.primitive == right.primitive &&
    left.terminal == right.terminal && left.placement == right.placement &&
    left.paths == right.paths

function PlacedRegion(source::Region, primitive::AbstractPrimitive)
    return PlacedRegion(
        source,
        primitive,
        nothing,
        (patterns = (),),
        ()
    )
end

boundary(region::PlacedRegion) = boundary(region.primitive)
area(region::PlacedRegion) = area(region.primitive)
centroid(region::PlacedRegion) = centroid(region.primitive)
support(region::PlacedRegion, φ::Real) = support(region.primitive, φ)

function resolve(context::Union{EmptyBoundary, AbstractPrimitive}, region::Region)
    primitive = resolve(context, region.primitive)
    placed = PlacedRegion(region, primitive)
    return CableGeometry(PlacedRegion[placed], boundary(placed))
end
