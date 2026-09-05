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
        boundary isa Sector && compact !== nothing && throw(ArgumentError(
            "sector bounded formations are intrinsically compacted"
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

function bounded_declarations(group::Group)
    group.item isa Stack || throw(ArgumentError(
        "a bounded formation must own an ordered Stack of strand declarations"
    ))
    parts = collect(group.item.items)
    for part in parts
        part isa Group || throw(ArgumentError(
            "a bounded formation Stack may contain only strand Groups"
        ))
        part.item isa Region || throw(ArgumentError(
            "each bounded strand declaration must own one Region"
        ))
        part.at == Pose2(0, 0, 0) || throw(ArgumentError(
            "bounded strand declarations must share the formation origin"
        ))
        part.boundary === nothing || throw(ArgumentError(
            "nested bounded formations are not supported"
        ))
        part.pattern === nothing || throw(ArgumentError(
            "a bounded strand inventory cannot declare a Ring or other placement"
        ))
        part.compact === nothing || throw(ArgumentError(
            "bounded compaction belongs to the containing formation"
        ))
    end
    return parts
end

function strand_path(paths, course::Int, maximum_course::Int)
    iszero(course) && return nothing
    paths === nothing && return nothing
    paths isa Tuple || return paths
    length(paths) == maximum_course || throw(DimensionMismatch(
        "lay schedule requires exactly $maximum_course inferred courses; got $(length(paths))"
    ))
    return paths[course]
end

function circular_members(boundary_shape::Disk, parts, compact::Bool)
    length(parts) == 2 || throw(ArgumentError(
        "a circular stranded core requires centre and outer strand declarations"
    ))
    centre_part, strand_part = parts
    centre = centre_part.item.primitive
    strand = strand_part.item.primitive
    centre isa Disk || throw(ArgumentError(
        "a circular stranded core requires a Disk centre wire"
    ))
    strand isa Disk || throw(ArgumentError(
        "circular strand packing requires Disk source wires"
    ))
    centre.r <= boundary_shape.r || throw(DomainError(
        centre.r, "the centre wire exceeds the circular core boundary"
    ))

    outer = circular_courses(boundary_shape, centre, strand, compact)
    sites = [(
        site = (boundary_shape.at.x, boundary_shape.at.y),
        course = 0,
        member = 1,
        angle = zero(boundary_shape.r)
    ); outer]
    maximum_course = maximum(getproperty.(sites, :course))
    members = NamedTuple[]
    for (member, candidate) in enumerate(sites)
        source = member == 1 ? centre_part.item : strand_part.item
        path = member == 1 ? nothing :
               strand_path(strand_part.path, candidate.course, maximum_course)
        angle = candidate.angle
        push!(members, (;
            source,
            candidate.course,
            member,
            path,
            angle,
            candidate.site
        ))
    end
    primitives = compact ? deform_disk_members(boundary_shape, members) :
                 [resolve(Pose2(member.site..., member.angle), member.source.primitive)
                  for member in members]
    return members, primitives
end

function sector_members(boundary_shape::SectorShape, parts)
    length(parts) == 1 || throw(ArgumentError(
        "a sector stranded core uses one declaration for its centre and course wires"
    ))
    part = only(parts)
    wire = part.item.primitive
    wire isa Disk || throw(ArgumentError(
        "a sector stranded core requires Disk source wires"
    ))
    courses, primitives, maximum_course = sector_courses(boundary_shape, wire)
    members = [(
                   source = part.item,
                   candidate.course,
                   member,
                   path = strand_path(part.path, candidate.course, maximum_course),
                   candidate.angle,
                   candidate.site
               ) for (member, candidate) in enumerate(courses)]
    return members, primitives
end

function rectangular_members(boundary_shape::Disk, parts, compact::Bool)
    compact || throw(ArgumentError(
        "rectangular strands require area-preserving bending"
    ))
    length(parts) == 2 || throw(ArgumentError(
        "a rectangular stranded core requires centre and strand declarations"
    ))
    centre_part, strand_part = parts
    centre = centre_part.item.primitive
    strand = strand_part.item.primitive
    centre isa Disk || throw(ArgumentError(
        "a rectangular stranded core requires a Disk centre wire"
    ))
    strand isa Rectangle || throw(ArgumentError(
        "rectangular strand packing requires Rectangle source strands"
    ))
    strips = rectangular_strands(boundary_shape, centre, strand)
    maximum_course = maximum(getproperty.(strips, :course))
    members = NamedTuple[(
        source = centre_part.item,
        course = 0,
        member = 1,
        path = nothing,
        angle = zero(boundary_shape.r),
        site = (boundary_shape.at.x, boundary_shape.at.y)
    )]
    primitives = AbstractShape[
        resolve(Pose2(boundary_shape.at.x, boundary_shape.at.y, boundary_shape.at.φ), centre)
    ]
    for strip in strips
        push!(members, (
            source = strand_part.item,
            strip.course,
            member = length(members) + 1,
            path = strand_path(strand_part.path, strip.course, maximum_course),
            angle = strip.primitive.at.φ,
            strip.site
        ))
        push!(primitives, strip.primitive)
    end
    return members, primitives
end

function bounded_members(boundary_shape::Disk, parts, compact::Bool)
    outer = last(parts).item.primitive
    return outer isa Disk ?
           circular_members(boundary_shape, parts, compact) :
           rectangular_members(boundary_shape, parts, compact)
end

bounded_members(boundary_shape::SectorShape, parts, ::Nothing) =
    sector_members(boundary_shape, parts)

function bounded_members(group::Group)
    boundary_shape = resolve(group.at, group.boundary)
    parts = bounded_declarations(group)
    members, primitives = bounded_members(boundary_shape, parts, group.compact)
    isempty(members) && throw(ArgumentError(
        "a bounded formation requires at least one strand"
    ))
    return boundary_shape, members, primitives
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
    outer, members, primitives = bounded_members(group)
    length(primitives) == length(members) || throw(DimensionMismatch(
        "bounded resolution must preserve the inferred strand count"
    ))
    formation_centre = centroid(outer)
    regions = PlacedRegion[]
    for (formation_member, (member, primitive)) in enumerate(zip(members, primitives))
        absolute_centre = centroid(primitive)
        pose = bounded_pose(member.site, group.at, member.angle)
        patterns = ((
            pattern = BoundedPlacement(outer, member.course),
            member = formation_member,
            pose = pose
        ),)
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
