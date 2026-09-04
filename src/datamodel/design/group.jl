"""
$(TYPEDEF)

Represent one repeated-member coalescing scope.

All conductive descendants resolve to terminal `name`. A group containing no
conductive descendant is a valid physical group and contributes no terminal.

$(TYPEDFIELDS)
"""
struct Group{A, E <: AbstractCablePart, P, H, C, B} <: AbstractCablePart
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
    "Authoritative finished boundary for a bounded formation, or `nothing`."
    boundary::B

    function Group(
            name::Symbol,
            at::A,
            item::E,
            pattern::P,
            path::H,
            compact::C,
            boundary::B
    ) where {A, E <: AbstractCablePart, P, H, C, B}
        isempty(String(name)) && throw(ArgumentError("group name cannot be empty"))
        at isa Pose2 || throw(ArgumentError("group pose must resolve to Pose2"))
        boundary === nothing || boundary isa Union{Disk, Sector} || throw(ArgumentError(
            "group boundary must be a Disk, Sector, or nothing"
        ))
        boundary === nothing || pattern === nothing || throw(ArgumentError(
            "a bounded formation owns its member placements and cannot carry a group pattern"
        ))
        boundary === nothing || path === nothing || throw(ArgumentError(
            "a bounded formation cannot carry a group-level longitudinal path"
        ))
        return new{A, E, P, H, C, B}(
            name, at, item, pattern, path, compact, boundary
        )
    end
end

Group(name::Symbol, at, item, pattern, path, compact) =
    Group(name, at, item, pattern, path, compact, nothing)

function Base.:(==)(left::Group, right::Group)
    left.name == right.name && left.at == right.at && left.item == right.item &&
        left.pattern == right.pattern && left.path == right.path &&
        left.compact == right.compact && left.boundary == right.boundary
end

function _path_radius(pattern::Ring, pose::Pose2, primitive::Annulus)
    return iszero(pattern.r) ? (r_in(primitive) + r_ex(primitive)) / 2 : pattern.r
end
_path_radius(pattern::Ring, pose::Pose2, primitive::AbstractShape) = pattern.r
function _path_radius(::Nothing, pose::Pose2, primitive::Annulus)
    (r_in(primitive) + r_ex(primitive)) / 2
end
_path_radius(::Nothing, pose::Pose2, primitive::AbstractShape) = hypot(pose.x, pose.y)
_path_radius(pattern, pose::Pose2, primitive::AbstractShape) = hypot(pose.x, pose.y)

function _resolved_path_radius(compact, pattern, pose, primitive)
    _path_radius(pattern, pose, primitive)
end
function _resolved_path_radius(
        compact, pattern, pose, primitive::Union{Polygon, BentStrip}
)
    centre = centroid(primitive)
    return hypot(centre...)
end
function _resolved_path_radius(::FillFactor, pattern, pose, primitive::Annulus)
    (r_in(primitive) + r_ex(primitive)) / 2
end

_member_definition(region::Region) = region.primitive
_member_definition(::AbstractCablePart) = nothing

function natural_sector_members(group::Group)
    group.item isa Stack && length(group.item.items) == 1 || throw(ArgumentError(
        "an uncompacted sector stranded formation uses one circular-wire inventory"
    ))
    part = only(group.item.items)
    part isa Group && part.item isa Region || throw(ArgumentError(
        "an uncompacted sector stranded formation requires one circular wire declaration"
    ))
    part.at == Pose2(0, 0, 0) || throw(ArgumentError(
        "sector wire packing owns every member placement"
    ))
    part.item.primitive isa Disk || throw(ArgumentError(
        "an uncompacted sector stranded formation requires circular wires"
    ))
    pattern = part.pattern
    pattern isa Ring || throw(ArgumentError(
        "an uncompacted sector stranded formation requires one wire inventory"
    ))
    boundary_shape = resolve(group.at, group.boundary)
    sites = sector_wire_sites(boundary_shape, part.item.primitive)
    count = pattern.n isa Int ? pattern.n : length(sites)
    sites = select_sector_sites(sites, count, centroid(boundary_shape))
    concrete = Ring(
        count;
        r = nothing,
        φ0 = pattern.φ0,
        span = pattern.span,
        gap_frac = pattern.gap_frac
    )
    centre = centroid(boundary_shape)
    return [(
                source = part.item,
                course = 1,
                member,
                pattern = concrete,
                path = part.path,
                angle = atan(site[2] - centre[2], site[1] - centre[1]),
                site
            ) for (member, site) in enumerate(sites)]
end

function bounded_members(group::Group)
    group.boundary isa Sector && group.compact === nothing &&
        return natural_sector_members(group)
    group.item isa Stack || throw(ArgumentError(
        "a bounded formation must own an ordered Stack of member courses"
    ))
    members = NamedTuple[]
    course = 0
    for (part_index, part) in enumerate(group.item.items)
        part isa Group || throw(ArgumentError(
            "a bounded formation Stack may contain only member Groups"
        ))
        part.item isa Region || throw(ArgumentError(
            "each bounded-formation course must repeat one Region"
        ))
        part.at == Pose2(0, 0, 0) || throw(ArgumentError(
            "bounded-formation courses must share the formation origin"
        ))
        part.boundary === nothing || throw(ArgumentError(
            "nested bounded formations are not supported"
        ))
        if part.pattern === nothing
            part_index == 1 || throw(ArgumentError(
                "only the first bounded-formation member may omit its course pattern"
            ))
            push!(members,
                (
                    source = part.item,
                    course,
                    member = 1,
                    pattern = nothing,
                    path = part.path,
                    angle = zero(eltype(group.at))
            ))
            continue
        end
        course += 1
        pattern = part.pattern
        pattern isa Ring || throw(ArgumentError(
            "bounded-formation outer courses require a Ring inventory declaration"
        ))
        pattern.n isa Int || throw(ArgumentError(
            "bounded-formation member counts must be explicit integers"
        ))
        pattern.r === nothing || throw(ArgumentError(
            "bounded-formation course radii are resolved by compaction"
        ))
        step = pattern.span / pattern.n
        for member in 1:pattern.n
            push!(members,
                (
                    source = part.item,
                    course,
                    member,
                    pattern,
                    path = part.path,
                    angle = pattern.φ0 + (member - 1) * step
                ))
        end
    end
    isempty(members) && throw(ArgumentError(
        "a bounded formation requires at least one member"
    ))
    return members
end

function bounded_pose(absolute, origin::Pose2, angle)
    dx = absolute[1] - origin.x
    dy = absolute[2] - origin.y
    cosine = cos(origin.φ)
    sine = sin(origin.φ)
    return Pose2(
        cosine * dx + sine * dy,
        -sine * dx + cosine * dy,
        angle
    )
end

function resolve_bounded(group::Group)
    members = bounded_members(group)
    primitives, outer = compact_bounded_members(
        group.boundary, members, group.compact, group.at
    )
    length(primitives) == length(members) || throw(DimensionMismatch(
        "bounded compaction must preserve the declared member count"
    ))
    formation_centre = centroid(outer)
    regions = PlacedRegion[]
    for (formation_member, (member, primitive)) in enumerate(zip(members, primitives))
        absolute_centre = centroid(primitive)
        pose = bounded_pose(absolute_centre, group.at, member.angle)
        declared = member.pattern === nothing ? () :
                   ((pattern = member.pattern, member = member.member, pose = pose),)
        patterns = (
            declared...,
            (
                pattern = BoundedPlacement(outer),
                member = formation_member,
                pose = pose
            )
        )
        paths = member.path === nothing ? () :
                ((
                    path = member.path,
                    radius = hypot(
                        absolute_centre[1] - formation_centre[1],
                        absolute_centre[2] - formation_centre[2]
                    )
                ),)
        terminal = member.source.material.kind === :conductor ? group.name : nothing
        push!(regions,
            PlacedRegion(
                member.source,
                primitive,
                terminal,
                (patterns = patterns,),
                paths
            ))
    end
    return CableGeometry(regions, outer)
end

_radial_half_extent(definition::Disk) = definition.r
_radial_half_extent(definition::Rectangle) = definition.h / 2
function _radial_half_extent(definition::AbstractPrimitive)
    support(resolve(EmptyBoundary(), definition))
end

function _contextual_ring(
        pattern::Ring,
        item::AbstractCablePart,
        child::CableGeometry,
        compact,
        context::Union{EmptyBoundary, AbstractShape}
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
function _contextual_pattern(
        pattern::Ring, item, child, compact, context
)
    _contextual_ring(pattern, item, child, compact, context)
end

function _group_placements(pattern, item, child, compact, context)
    concrete = _contextual_pattern(pattern, item, child, compact, context)
    subject = _member_definition(item)
    subject === nothing && (subject = child)
    return concrete, placements(concrete, subject, compact)
end

function _minimum_radius(primitive::Annulus)
    iszero(primitive.at.x) && iszero(primitive.at.y) && return r_in(primitive)
    return _minimum_radius_general(primitive)
end
function _minimum_radius(primitive::BentStrip)
    iszero(primitive.at.x) && iszero(primitive.at.y) && return primitive.ri
    return _minimum_radius_general(primitive)
end
_minimum_radius(primitive::AbstractShape) = _minimum_radius_general(primitive)
function _minimum_radius_general(primitive::AbstractShape)
    center = centroid(primitive)
    center_radius = hypot(center...)
    iszero(center_radius) && return zero(eltype(primitive))
    φ = atan(center[2], center[1])
    return -support(primitive, φ + pi)
end

function _resolve_group(
        context::Union{EmptyBoundary, AbstractShape},
        group::Group
)
    if group.boundary !== nothing
        context isa EmptyBoundary || throw(ArgumentError(
            "a bounded stranded formation must be the innermost physical formation"
        ))
        return resolve_bounded(group)
    end
    if group.pattern === nothing
        child = resolve(context, group.item)
        regions = PlacedRegion[]
        for source in child.regions
            placed = resolve(group.at, source)
            terminal = source.source.material.kind === :conductor ? group.name :
                       source.terminal
            centre = centroid(source.primitive)
            paths = group.path === nothing ? source.paths :
                    (source.paths...,
                (
                    path = group.path,
                    radius = _path_radius(
                        nothing,
                        Pose2(centre[1], centre[2], 0),
                        source.primitive
                    )
                ))
            push!(regions, PlacedRegion(
                source.source,
                placed.primitive,
                terminal,
                placed.placement,
                paths
            ))
        end
        outer = group.path === nothing ? boundary(child) : Disk(support(boundary(child)))
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
            placed = resolve(placed_at, source)
            extent = support(local_primitive)
            local_extent = local_extent === nothing ? extent : max(local_extent, extent)
            patterns = pattern === nothing ? placed.placement.patterns :
                       (placed.placement.patterns...,
                (pattern = pattern, member = member, pose = pose))
            paths = group.path === nothing ? source.paths :
                    (source.paths...,
                (
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

function resolve(context::AbstractShape, group::Group)
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
