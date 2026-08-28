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
        primitive::Union{Annulus, Sector}
)
    return iszero(pattern.r) ? (r_in(primitive) + r_ex(primitive)) / 2 : pattern.r
end
_path_radius(pattern::Ring, pose::Pose2, primitive::AbstractPrimitive) = pattern.r
_path_radius(
    ::Nothing,
    pose::Pose2,
    primitive::Union{Annulus, Sector}
) = (r_in(primitive) + r_ex(primitive)) / 2
_path_radius(::Nothing, pose::Pose2, primitive::AbstractPrimitive) = hypot(pose.x, pose.y)
_path_radius(pattern, pose::Pose2, primitive::AbstractPrimitive) = hypot(pose.x, pose.y)

function _minimum_radius(primitive::Union{Annulus, Sector})
    iszero(primitive.at.x) && iszero(primitive.at.y) && return r_in(primitive)
    return _minimum_radius_general(primitive)
end
_minimum_radius(primitive::AbstractPrimitive) = _minimum_radius_general(primitive)
function _minimum_radius_general(primitive::AbstractPrimitive)
    center = centroid(primitive)
    center_radius = hypot(center...)
    iszero(center_radius) && return zero(eltype(primitive))
    φ = atan(center[2], center[1])
    return -support(primitive, φ + pi)
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
            primitive = resolve(placed_at, source.primitive)
            patterns = group.pattern === nothing ? source.placement.patterns :
                       (source.placement.patterns...,
                        (pattern = group.pattern, member = member, pose = pose))
            paths = group.path === nothing ? source.paths :
                    (source.paths..., (
                        path = group.path,
                        radius = _path_radius(group.pattern, pose, source.primitive)
                    ))
            push!(regions, PlacedRegion(
                source.source,
                primitive,
                terminal,
                (patterns = patterns,),
                paths
            ))
        end
    end
    has_conductor || throw(ArgumentError(
        "group :$(group.name) has no conductive descendant"
    ))

    child_extent = support(boundary(child))
    outer_radius = maximum(hypot(pose.x, pose.y) + child_extent for pose in poses)
    local_boundary = Disk(outer_radius)
    return CableGeometry(regions, resolve(group.at, local_boundary))
end

function resolve(context::AbstractPrimitive, group::Group)
    result = resolve(EmptyBoundary(), group)
    current_radius = support(context)
    tolerance = sqrt(eps(typeof(float(current_radius)))) *
                max(one(current_radius), current_radius)
    for source in result.regions
        minimum_radius = _minimum_radius(source.primitive)
        minimum_radius + tolerance >= current_radius || throw(DomainError(
            minimum_radius,
            "group geometry overlaps the current stack boundary at radius $current_radius"
        ))
    end
    return result
end
