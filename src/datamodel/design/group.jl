"""
$(TYPEDEF)

Represent one repeated-member coalescing scope.

Each conductive descendant resolves to terminal `name`. A group containing no
conductive descendant is a valid physical group and contributes no terminal.

$(TYPEDFIELDS)
"""
struct Group{
    A <: Pose2,
    E <: AbstractCablePart,
    P,
    H,
    C,
    B <: Union{Nothing, Disk, Sector}
} <: AbstractCablePart
    "Coalesced terminal name when the group contains conductive descendants."
    name::Symbol
    "Pose relative to the containing frame."
    at::A
    "Repeated physical member."
    member::E
    "Cross-sectional placement definition."
    pattern::P
    "Longitudinal path definition or `nothing`."
    path::H
    "Compaction definition or `nothing`."
    compact::C
    "Formation boundary, or rectangular-course packing limit. `nothing` for unbounded groups."
    boundary::B

    function Group(
            name::Symbol,
            at::A,
            member::E,
            pattern::P,
            path::H,
            compact::C,
            boundary::B
    ) where {A <: Pose2, E <: AbstractCablePart, P, H, C, B <: Union{Nothing, Disk}}
        return validate(new{A, E, P, H, C, B}(
            name, at, member, pattern, path, compact, boundary
        ))
    end

    # The sector outline compacts a sector formation, so it takes no compaction.
    function Group(
            name::Symbol,
            at::A,
            member::E,
            pattern::P,
            path::H,
            compact::Nothing,
            boundary::B
    ) where {A <: Pose2, E <: AbstractCablePart, P, H, B <: Sector}
        return validate(new{A, E, P, H, Nothing, B}(
            name, at, member, pattern, path, compact, boundary
        ))
    end
end

function validate(group::Group)
    group.name === Symbol("") && throw(ArgumentError("group name cannot be empty"))
    group.boundary === nothing || group.pattern === nothing || throw(ArgumentError(
        "a bounded formation owns its member placements and cannot carry a group pattern"
    ))
    group.boundary === nothing || group.path === nothing || throw(ArgumentError(
        "a bounded formation cannot carry a group-level longitudinal path"
    ))
    return group
end

Group(name::Symbol, at, member, pattern, path, compact) =
    Group(name, at, member, pattern, path, compact, nothing)

"""
Declare a repeated-member scope at the local origin, or a finite space of such scopes.
"""
function Group(
        name::Union{Symbol, AbstractGrid, Gridspace},
        member::Union{AbstractCablePart, AbstractGrid, Gridspace};
        at = Pose2(0, 0, 0),
        pattern = nothing,
        path = nothing,
        compact = nothing,
        boundary = nothing,
        combine::Symbol = :product
)
    values = (name, at, member, pattern, path, compact, boundary)
    return parameterize(Group, Group, values; combine)
end

function Base.:(==)(left::Group, right::Group)
    left.name == right.name && left.at == right.at && left.member == right.member &&
        left.pattern == right.pattern && left.path == right.path &&
        left.compact == right.compact && left.boundary == right.boundary
end

function _path_radius(pattern::Ring, pose::Pose2, primitive::Annulus)
    return iszero(pattern.r) ? (r_in(primitive) + r_ex(primitive)) / 2 : pattern.r
end
_path_radius(pattern::Ring, pose::Pose2, primitive::AbstractShape) = pattern.r
function _path_radius(::Nothing, pose::Pose2, primitive::Union{Annulus, BentStrip})
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
    center = centroid(primitive)
    return hypot(center...)
end
function _resolved_path_radius(::FillFactor, pattern, pose, primitive::Annulus)
    (r_in(primitive) + r_ex(primitive)) / 2
end

# A pattern places a Region by its primitive and any other part by its resolved geometry.
_member_definition(region::Region, child) = region.primitive
_member_definition(::AbstractCablePart, child) = child

# A cross-section that the course derives, such as a bent strip, replaces the primitive
# of a Region member. A declared primitive or a part's geometry keeps each resolved
# region primitive.
_member_primitive(shape::AbstractPrimitive, placed) = placed.primitive
_member_primitive(shape::CableGeometry, placed) = placed.primitive
_member_primitive(shape::AbstractShape, placed) = shape

function bounded_declarations(group::Group{<:Pose2, <:Stack})
    parts = collect(group.member.items)
    for part in parts
        validate(part, group)
    end
    return parts
end

function bounded_declarations(::Group)
    throw(ArgumentError(
        "a bounded formation must own an ordered Stack of strand declarations"
    ))
end

# A strand declaration of a bounded formation is a Group at the formation origin that owns
# one Region and leaves placement and compaction to the formation.
function validate(part::Group{<:Pose2, <:Region}, ::Group)
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
    return part
end

function validate(part::Group, ::Group)
    throw(ArgumentError(
        "each bounded strand declaration must own one Region; " *
        ":$(part.name) owns $(nameof(typeof(part.member)))"
    ))
end

function validate(part::AbstractCablePart, ::Group)
    throw(ArgumentError(
        "a bounded formation Stack may contain only strand Groups; " *
        "got $(nameof(typeof(part)))"
    ))
end

# A course follows no path at the center. Otherwise it follows the one path of every
# course, or its own entry of a lay schedule with one path per inferred course.
strand_path(path, course::Int, maximum_course::Int) = iszero(course) ? nothing : path
strand_path(::Nothing, course::Int, maximum_course::Int) = nothing
function strand_path(paths::Tuple, course::Int, maximum_course::Int)
    iszero(course) && return nothing
    length(paths) == maximum_course || throw(DimensionMismatch(
        "lay schedule requires exactly $maximum_course inferred courses; got $(length(paths))"
    ))
    return paths[course]
end

"""
$(TYPEDSIGNATURES)

Construct the members and placed primitives of a bounded wire formation from
its physical boundary, ordered part declarations, and compaction selection.
Coordinates and dimensions are in meters. Return `(boundary, members, primitives)` in
construction order. The center and strand declarations select the formation. Disk wires
pack in circular courses, and a Rectangle strand bends around a Disk center wire.
"""
function bounded_members(boundary_shape::Disk, parts, compact::Bool)
    length(parts) == 2 || throw(ArgumentError(
        "a disk-bounded stranded core requires center and strand declarations"
    ))
    center_part, strand_part = parts
    return bounded_members(boundary_shape, center_part, strand_part, compact)
end

function bounded_members(
        boundary_shape::Disk,
        center_part::Group{<:Pose2, <:Region{<:Disk}},
        strand_part::Group{<:Pose2, <:Region{<:Disk}},
        compact::Bool
)
    center = center_part.member.primitive
    strand = strand_part.member.primitive
    center.r <= boundary_shape.r || throw(DomainError(
        center.r, "the center wire exceeds the circular core boundary"
    ))

    outer = circular_courses(boundary_shape, center, strand, compact)
    sites = [(
        site = (boundary_shape.at.x, boundary_shape.at.y),
        course = 0,
        member = 1,
        angle = zero(boundary_shape.r)
    ); outer]
    maximum_course = maximum(getproperty.(sites, :course))
    members = NamedTuple[]
    for (member, placement) in enumerate(sites)
        source = member == 1 ? center_part.member : strand_part.member
        path = member == 1 ? nothing :
               strand_path(strand_part.path, placement.course, maximum_course)
        angle = placement.angle
        push!(members, (;
            source,
            placement.course,
            member,
            path,
            angle,
            placement.site
        ))
    end
    primitives = compact ? deform_disk_members(boundary_shape, members) :
                 [resolve(Pose2(member.site..., member.angle), member.source.primitive)
                  for member in members]
    return boundary_shape, members, primitives
end

function bounded_members(
        boundary_shape::Disk,
        center_part::Group{<:Pose2, <:Region{<:Disk}},
        strand_part::Group{<:Pose2, <:Region{<:Rectangle}},
        compact::Bool
)
    compact || throw(ArgumentError(
        "rectangular strands require area-preserving bending"
    ))
    center = center_part.member.primitive
    strand = strand_part.member.primitive
    strips = rectangular_strands(boundary_shape, center, strand)
    maximum_course = maximum(getproperty.(strips, :course))
    members = NamedTuple[(
        source = center_part.member,
        course = 0,
        member = 1,
        path = nothing,
        angle = zero(boundary_shape.r),
        site = (boundary_shape.at.x, boundary_shape.at.y)
    )]
    primitives = AbstractShape[
        resolve(Pose2(boundary_shape.at.x, boundary_shape.at.y, boundary_shape.at.φ), center)
    ]
    for strip in strips
        push!(members, (
            source = strand_part.member,
            strip.course,
            member = length(members) + 1,
            path = strand_path(strand_part.path, strip.course, maximum_course),
            angle = strip.primitive.at.φ,
            strip.site
        ))
        push!(primitives, strip.primitive)
    end
    # Complete bent courses occupy exactly this disk. The requested packing limit is not an
    # additional physical layer.
    return Disk(r_ex(last(primitives)), boundary_shape.at), members, primitives
end

function sector_members(boundary_shape::SectorShape, parts)
    length(parts) == 1 || throw(ArgumentError(
        "a sector stranded core uses one declaration for its center and course wires"
    ))
    part = only(parts)
    wire = part.member.primitive
    courses, primitives, maximum_course = sector_courses(boundary_shape, wire)
    members = [(
                   source = part.member,
                   placement.course,
                   member,
                   path = strand_path(part.path, placement.course, maximum_course),
                   placement.angle,
                   placement.site
               ) for (member, placement) in enumerate(courses)]
    return members, primitives
end

function bounded_members(boundary_shape::SectorShape, parts, ::Nothing)
    return boundary_shape, sector_members(boundary_shape, parts)...
end

"""
$(TYPEDSIGNATURES)

Resolve a bounded wire group into `(boundary_shape, members, primitives)`.
The returned geometric boundary follows the group's pose and the occupied rectangular
courses when present. Coordinates and dimensions are in meters. Member and
primitive order follows construction. An empty formation raises `ArgumentError`.
"""
function bounded_members(group::Group)
    boundary_shape = resolve(group.at, group.boundary)
    parts = bounded_declarations(group)
    outer, members, primitives = bounded_members(boundary_shape, parts, group.compact)
    isempty(members) && throw(ArgumentError(
        "a bounded formation requires at least one strand"
    ))
    return outer, members, primitives
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
    formation_center = centroid(outer)
    regions = PlacedRegion[]
    for (formation_member, (member, primitive)) in enumerate(zip(members, primitives))
        absolute_center = centroid(primitive)
        pose = bounded_pose(member.site, group.at, member.angle)
        patterns = ((
            owner = Group,
            pattern = BoundedPlacement(outer, member.course),
            member = formation_member,
            pose = pose
        ),)
        paths = member.path === nothing ? () :
                ((
                    path = member.path,
                    radius = _path_radius(
                        nothing,
                        Pose2(absolute_center[1] - formation_center[1],
                            absolute_center[2] - formation_center[2]),
                        primitive
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
_radial_half_extent(geometry::CableGeometry) = support(boundary(geometry))

# Members on a contextual radius start at the center of an empty context, or at the
# support radius of the enclosing geometry.
_context_radius(::EmptyBoundary, child) = zero(support(boundary(child)))
_context_radius(context::AbstractShape, child) = support(context)

_contextual_pattern(pattern, item, child, compact, context) = pattern

# A contextual radius lays the member centers one radial half extent outside the context.
function _contextual_pattern(
        pattern::Ring{<:Any, <:Real, Nothing},
        item,
        child,
        compact,
        context
)
    radial = _radial_half_extent(_member_definition(item, child))
    ring = Ring(
        pattern.n;
        r = _context_radius(context, child) + radial,
        φ0 = pattern.φ0,
        span = pattern.span,
        gap_frac = pattern.gap_frac
    )
    return _contextual_pattern(ring, item, child, compact, context)
end

# `capacity()` places as many members as fit on the resolved radius.
function _contextual_pattern(
        pattern::Ring{_DeferredCardinality, <:Real, <:Real},
        item,
        child,
        compact,
        context
)
    count = capacity(pattern, _member_definition(item, child), compact)
    count > 0 || throw(ArgumentError(
        "the group geometry cannot admit one member"
    ))
    return Ring(
        count;
        r = pattern.r,
        φ0 = pattern.φ0,
        span = pattern.span,
        gap_frac = pattern.gap_frac
    )
end

function _minimum_radius(primitive::Union{Annulus, BentStrip})
    iszero(primitive.at.x) && iszero(primitive.at.y) && return r_in(primitive)
    return invoke(_minimum_radius, Tuple{AbstractShape}, primitive)
end
function _minimum_radius(primitive::AbstractShape)
    center = centroid(primitive)
    center_radius = hypot(center...)
    iszero(center_radius) && return zero(eltype(primitive))
    φ = atan(center[2], center[1])
    return -support(primitive, φ + pi)
end

# A bounded formation is the innermost physical formation.
function _resolve_group(
        ::EmptyBoundary,
        group::Group{<:Pose2, <:Any, <:Any, <:Any, <:Any, <:Union{Disk, Sector}}
)
    return resolve_bounded(group)
end
function _resolve_group(
        ::AbstractShape,
        ::Group{<:Pose2, <:Any, <:Any, <:Any, <:Any, <:Union{Disk, Sector}}
)
    throw(ArgumentError(
        "a bounded stranded formation must be the innermost physical formation"
    ))
end

function _resolve_group(
        context::Union{EmptyBoundary, AbstractShape},
        group::Group{<:Pose2, <:Any, <:Any, <:Any, <:Any, Nothing}
)
    if group.pattern === nothing
        child = resolve(context, group.member)
        regions = PlacedRegion[]
        for placed in child.regions
            resolved = resolve(group.at, placed)
            terminal = placed.source.material.kind === :conductor ? group.name :
                       placed.terminal
            # This scope coalesces even terminal-preserving child assemblies.
            patterns = map(resolved.placement.patterns) do entry
                entry.owner === Assembly ? merge(entry, (owner = Group,)) : entry
            end
            center = centroid(placed.primitive)
            paths = group.path === nothing ? placed.paths :
                    (placed.paths...,
                (
                    path = group.path,
                    radius = _path_radius(
                        nothing,
                        Pose2(center[1], center[2], 0),
                        placed.primitive
                    )
                ))
            push!(regions, PlacedRegion(
                placed.source,
                resolved.primitive,
                terminal,
                (patterns = patterns,),
                paths
            ))
        end
        outer = group.path === nothing ? boundary(child) : Disk(support(boundary(child)))
        return CableGeometry(regions, resolve(group.at, outer))
    end

    child = resolve(EmptyBoundary(), group.member)
    pattern = _contextual_pattern(group.pattern, group.member, child, group.compact, context)
    course = placements(pattern, _member_definition(group.member, child), group.compact)
    # The local primitives of the members are resolved before the regions are built, so the
    # enclosing disk is known before the loop.
    local_primitives = [resolve(pose, _member_primitive(course.shape, placed))
                        for placed in child.regions, pose in course.poses]
    isempty(local_primitives) && throw(ArgumentError("group placement cannot be empty"))
    local_boundary = Disk(foldl(max, support(primitive) for primitive in local_primitives))

    regions = PlacedRegion[]
    for (member, pose) in enumerate(course.poses)
        placed_at = group.at * pose
        for (index, placed) in enumerate(child.regions)
            terminal = placed.source.material.kind === :conductor ? group.name :
                       placed.terminal
            local_primitive = local_primitives[index, member]
            primitive = resolve(group.at, local_primitive)
            resolved = resolve(placed_at, placed)
            patterns = map(resolved.placement.patterns) do entry
                entry.owner === Assembly ? merge(entry, (owner = Group,)) : entry
            end
            patterns = (patterns...,
                (owner = Group, pattern = pattern, member = member, pose = pose))
            paths = group.path === nothing ? placed.paths :
                    (placed.paths...,
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
                placed.source,
                primitive,
                terminal,
                (patterns = patterns,),
                paths
            ))
        end
    end
    return CableGeometry(regions, resolve(group.at, local_boundary))
end

resolve(context::EmptyBoundary, group::Group) = _resolve_group(context, group)

function resolve(context::AbstractShape, group::Group)
    result = _resolve_group(context, group)
    current_radius = support(context)
    tolerance = sqrt(eps(typeof(float(current_radius)))) *
                max(one(current_radius), current_radius)
    for placed in result.regions
        minimum_radius = _minimum_radius(placed.primitive)
        minimum_radius + tolerance >= current_radius || throw(DomainError(
            minimum_radius,
            "group geometry overlaps the current stack boundary at radius $current_radius"
        ))
    end
    return result
end
