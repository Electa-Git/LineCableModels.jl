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

function resolve(context::EmptyBoundary, group::Group)
    child = resolve(EmptyBoundary(), group.item)
    poses = placements(group.pattern, child, group.compact)
    isempty(poses) && throw(ArgumentError("group placement cannot be empty"))

    regions = ResolvedRegion[]
    has_conductor = false
    for pose in poses
        placed_at = group.at * pose
        radius = if group.pattern isa Ring && iszero(group.pattern.r) &&
                    length(child.regions) == 1 &&
                    only(child.regions).region.primitive isa Union{Annulus, Sector}
            primitive = only(child.regions).region.primitive
            (primitive.ri + primitive.ro) / 2
        elseif group.pattern isa Ring
            group.pattern.r
        elseif group.pattern === nothing && length(child.regions) == 1 &&
                    only(child.regions).region.primitive isa Union{Annulus, Sector}
            primitive = only(child.regions).region.primitive
            (primitive.ri + primitive.ro) / 2
        else
            hypot(pose.x, pose.y)
        end
        path_factor = group.path === nothing ? one(eltype(placed_at)) :
                      overlength(group.path, radius)
        local_turns = group.path === nothing ? zero(eltype(placed_at)) :
                      inv(pitch(group.path, radius))
        pattern_depth = group.pattern === nothing ? 0 : 1
        path_depth = group.path === nothing ? 0 : 1
        for source in child.regions
            terminal = source.region.material.kind === :conductor ? group.name :
                       source.terminal
            has_conductor |= terminal === group.name
            push!(regions, ResolvedRegion(
                source.region,
                PlacedShape(source.shape, placed_at),
                terminal,
                source.overlength * path_factor,
                source.turns + local_turns,
                source.pattern_depth + pattern_depth,
                source.path_depth + path_depth
            ))
        end
    end
    has_conductor || throw(ArgumentError(
        "group :$(group.name) has no conductive descendant"
    ))

    child_extent = support(boundary(child))
    outer_radius = maximum(hypot(pose.x, pose.y) + child_extent for pose in poses)
    local_boundary = DiskShape(outer_radius)
    return ResolvedPart(regions, PlacedShape(local_boundary, group.at), [group.name])
end

function resolve(context::AbstractShape, group::Group)
    result = resolve(EmptyBoundary(), group)
    current_radius = support(context)
    tolerance = sqrt(eps(typeof(float(current_radius)))) *
                max(one(current_radius), current_radius)
    for source in result.regions
        center = centroid(source.shape)
        region_radius = sqrt(area(source.shape) / (one(current_radius) * pi))
        minimum_radius = source.region.primitive isa Disk ?
                         hypot(center...) - region_radius :
                         source.region.primitive isa Union{Annulus, Sector} ?
                         r_in(source.shape) : -oftype(current_radius, Inf)
        minimum_radius + tolerance >= current_radius || throw(DomainError(
            minimum_radius,
            "group geometry overlaps the current stack boundary at radius $current_radius"
        ))
    end
    if group.pattern isa Ring && group.item isa Region &&
            group.item.primitive isa Disk && iszero(group.at.x) &&
            iszero(group.at.y) && iszero(group.at.φ)
        member_radius = group.item.primitive.r
        declared_inner = group.pattern.r - member_radius
        if isapprox(declared_inner, current_radius)
            outer = PlacedShape(
                DiskShape(current_radius + 2member_radius),
                group.at
            )
            return ResolvedPart(result.regions, outer, result.terminals)
        end
    end
    return result
end
