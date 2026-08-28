"""
$(TYPEDEF)

Represent one retained electrical terminal and its physical members.

$(TYPEDFIELDS)
"""
struct Group{A, E <: AbstractCablePart, P, H, C} <: AbstractCablePart
    "Retained terminal name."
    name::Symbol
    "Pose relative to the containing frame."
    at::A
    "Repeated physical member."
    item::E
    "Cross-sectional placement definition."
    pattern::P
    "Longitudinal path definition or `nothing`."
    path::H
    "Compaction definition or `nothing`."
    compact::C

    function Group(
            name::Symbol,
            at::A,
            item::E,
            pattern::P,
            path::H,
            compact::C
    ) where {A, E <: AbstractCablePart, P, H, C}
        isempty(String(name)) && throw(ArgumentError("group name cannot be empty"))
        at isa Pose2 || throw(ArgumentError("group pose must resolve to Pose2"))
        return new{A, E, P, H, C}(name, at, item, pattern, path, compact)
    end
end

function _path_radius(
        pattern::Ring,
        pose::Pose2,
        shape::Union{
            AnnulusShape,
            SectorShape,
            PlacedShape{<:Real, <:Union{AnnulusShape, SectorShape}}
        }
)
    return iszero(pattern.r) ? (r_in(shape) + r_ex(shape)) / 2 : pattern.r
end
_path_radius(pattern::Ring, pose::Pose2, shape::AbstractShape) = pattern.r
_path_radius(
    ::Nothing,
    pose::Pose2,
    shape::Union{
        AnnulusShape,
        SectorShape,
        PlacedShape{<:Real, <:Union{AnnulusShape, SectorShape}}
    }
) = (r_in(shape) + r_ex(shape)) / 2
_path_radius(::Nothing, pose::Pose2, shape::AbstractShape) = hypot(pose.x, pose.y)
_path_radius(pattern, pose::Pose2, shape::AbstractShape) = hypot(pose.x, pose.y)

_minimum_radius(
    shape::PlacedShape{<:Real, <:Union{AnnulusShape, SectorShape}}
) = r_in(shape)
function _minimum_radius(shape::PlacedShape)
    center = centroid(shape)
    center_radius = hypot(center...)
    iszero(center_radius) && return zero(eltype(shape))
    φ = atan(center[2], center[1])
    return -support(shape, φ + pi)
end

function resolve(context::EmptyBoundary, group::Group)
    child = resolve(EmptyBoundary(), group.item)
    poses = placements(group.pattern, child, group.compact)
    isempty(poses) && throw(ArgumentError("group placement cannot be empty"))

    regions = PlacedRegion[]
    has_conductor = false
    for (member, pose) in enumerate(poses)
        placed_at = group.at * pose
        for source in child.regions
            terminal = source.source.material.kind === :conductor ? group.name :
                       source.terminal
            has_conductor |= terminal === group.name
            shape = PlacedShape(source.shape, placed_at)
            patterns = group.pattern === nothing ? source.placement.patterns :
                       (source.placement.patterns...,
                        (pattern = group.pattern, member = member, pose = pose))
            paths = group.path === nothing ? source.paths :
                    (source.paths..., (
                        path = group.path,
                        radius = _path_radius(group.pattern, pose, source.shape)
                    ))
            push!(regions, PlacedRegion(
                source.source,
                shape,
                terminal,
                (pose = shape.at, patterns = patterns),
                paths
            ))
        end
    end
    has_conductor || throw(ArgumentError(
        "group :$(group.name) has no conductive descendant"
    ))

    child_extent = support(boundary(child))
    outer_radius = maximum(hypot(pose.x, pose.y) + child_extent for pose in poses)
    local_boundary = DiskShape(outer_radius)
    return CableGeometry(regions, PlacedShape(local_boundary, group.at))
end

function resolve(context::AbstractShape, group::Group)
    result = resolve(EmptyBoundary(), group)
    current_radius = support(context)
    tolerance = sqrt(eps(typeof(float(current_radius)))) *
                max(one(current_radius), current_radius)
    for source in result.regions
        minimum_radius = _minimum_radius(source.shape)
        minimum_radius + tolerance >= current_radius || throw(DomainError(
            minimum_radius,
            "group geometry overlaps the current stack boundary at radius $current_radius"
        ))
    end
    return result
end
