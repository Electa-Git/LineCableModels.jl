"""
$(TYPEDEF)

Represent one repeated-member coalescing scope.

All conductive descendants resolve to terminal `name`. A group containing no
conductive descendant is a valid physical group and contributes no terminal.

$(TYPEDFIELDS)
"""
struct Group{A, E <: AbstractCablePart, P, H, C} <: AbstractCablePart
    "Coalesced terminal name when the group contains conductive descendants."
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

Base.:(==)(left::Group, right::Group) =
    left.name == right.name && left.at == right.at && left.item == right.item &&
    left.pattern == right.pattern && left.path == right.path &&
    left.compact == right.compact

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

_resolved_path_radius(compact, pattern, pose, primitive) =
    _path_radius(pattern, pose, primitive)
_resolved_path_radius(
    ::FillFactor, pattern, pose, primitive::Union{Annulus, Sector}
) = (r_in(primitive) + r_ex(primitive)) / 2

_member_definition(region::Region) = region.primitive
_member_definition(::AbstractCablePart) = nothing

_radial_half_extent(definition::DiskDefinition) = definition.r
_radial_half_extent(definition::RectangleDefinition) = definition.h / 2
_radial_half_extent(definition::AbstractPrimitiveDefinition) =
    support(resolve(EmptyBoundary(), definition))

function _contextual_ring(
        pattern::Ring,
        item::AbstractCablePart,
        child::CableGeometry,
        compact,
        context::Union{EmptyBoundary, AbstractPrimitive}
)
    definition = _member_definition(item)
    inner = context isa EmptyBoundary ? zero(support(boundary(child))) : support(context)
    radial = definition === nothing ? support(boundary(child)) :
             _radial_half_extent(definition)
    radius = something(pattern.r, inner + radial)
    count = if pattern.n isa Int
        pattern.n
    elseif definition === nothing
        capacity(Ring, radius, support(boundary(child)); gap_frac = pattern.gap_frac)
    else
        capacity(
            Ring(pattern.n; r = radius, φ0 = pattern.φ0,
                span = pattern.span, gap_frac = pattern.gap_frac),
            definition,
            compact
        )
    end
    count > 0 || throw(ArgumentError(
        "the group geometry cannot admit one member"
    ))
    return Ring(
        count;
        r = radius,
        φ0 = pattern.φ0,
        span = pattern.span,
        gap_frac = pattern.gap_frac
    )
end

_contextual_pattern(pattern, item, child, compact, context) = pattern
_contextual_pattern(
    pattern::Ring, item, child, compact, context
) = _contextual_ring(pattern, item, child, compact, context)

function _group_placements(pattern, item, child, compact, context)
    concrete = _contextual_pattern(pattern, item, child, compact, context)
    subject = _member_definition(item)
    subject === nothing && (subject = child)
    return concrete, placements(concrete, subject, compact)
end

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

function _resolve_group(
        context::Union{EmptyBoundary, AbstractPrimitive},
        group::Group
)
    if group.pattern === nothing
        child = resolve(context, group.item)
        regions = PlacedRegion[]
        for source in child.regions
            primitive = resolve(group.at, source.primitive)
            terminal = source.source.material.kind === :conductor ? group.name :
                       source.terminal
            centre = centroid(source.primitive)
            paths = group.path === nothing ? source.paths :
                    (source.paths..., (
                        path = group.path,
                        radius = _path_radius(
                            nothing,
                            Pose2(centre[1], centre[2], 0),
                            source.primitive
                        )
                    ))
            push!(regions, PlacedRegion(
                source.source,
                primitive,
                terminal,
                source.placement,
                paths
            ))
        end
        outer = group.path === nothing ? boundary(child) :
                Disk(support(boundary(child)))
        return CableGeometry(regions, resolve(group.at, outer))
    end

    child = resolve(EmptyBoundary(), group.item)
    pattern, members = _group_placements(
        group.pattern, group.item, child, group.compact, context
    )
    isempty(members) && throw(ArgumentError("group placement cannot be empty"))

    regions = PlacedRegion[]
    local_extent = nothing
    for (member, placement) in enumerate(members)
        pose = _placement_pose(placement)
        placed_at = group.at * pose
        for source in child.regions
            terminal = source.source.material.kind === :conductor ? group.name :
                       source.terminal
            definition = _placement_definition(
                placement,
                source.source.primitive
            )
            local_primitive = placement isa _ResolvedPlacement ?
                              resolve(pose, definition) :
                              resolve(pose, source.primitive)
            primitive = resolve(group.at, local_primitive)
            extent = support(local_primitive)
            local_extent = local_extent === nothing ? extent : max(local_extent, extent)
            patterns = pattern === nothing ? source.placement.patterns :
                       (source.placement.patterns...,
                        (pattern = pattern, member = member, pose = pose))
            paths = group.path === nothing ? source.paths :
                    (source.paths..., (
                        path = group.path,
                        radius = _resolved_path_radius(
                            group.compact,
                            pattern,
                            pose,
                            local_primitive
                        )
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
    local_extent === nothing && throw(ArgumentError("group placement cannot be empty"))
    local_boundary = Disk(local_extent)
    return CableGeometry(regions, resolve(group.at, local_boundary))
end

resolve(context::EmptyBoundary, group::Group) = _resolve_group(context, group)

function resolve(context::AbstractPrimitive, group::Group)
    result = _resolve_group(context, group)
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
