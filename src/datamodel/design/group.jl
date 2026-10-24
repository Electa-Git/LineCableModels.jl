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
    B <: Union{Nothing, Disk, Sector, Rectangle, Ellipse, Polygon}
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
    "Formation boundary, the finished outline. `nothing` for unbounded groups."
    boundary::B

    function Group(
            name::Symbol,
            at::A,
            member::E,
            pattern::P,
            path::H,
            compact::C,
            boundary::B
    ) where {
            A <: Pose2, E <: AbstractCablePart,
            P <: PlacementPattern,
            H <: LongitudinalPath,
            C <: Compaction, B <: Union{Nothing, Disk}
    }
        return validate(new{A, E, P, H, C, B}(
            name, at, member, pattern, path, compact, boundary
        ))
    end

    # A non-circular outline deforms its formation, so it takes no compaction.
    function Group(
            name::Symbol,
            at::A,
            member::E,
            pattern::P,
            path::H,
            compact::Nothing,
            boundary::B
    ) where {
            A <: Pose2, E <: AbstractCablePart,
            P <: PlacementPattern,
            H <: LongitudinalPath,
            B <: Union{Sector, Rectangle, Ellipse, Polygon}
    }
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

function Group(name::Symbol, at::Pose2, member::AbstractCablePart,
        pattern::PlacementPattern,
        path::LongitudinalPath,
        compact::Compaction)
    return Group(name, at, member, pattern, path, compact, nothing)
end

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
# one Region. Its pattern declares its course counts, and the formation owns placement and
# compaction.
function validate(part::Group{<:Pose2, <:Region}, ::Group)
    part.at == Pose2(0, 0, 0) || throw(ArgumentError(
        "bounded strand declarations must share the formation origin"
    ))
    part.boundary === nothing || throw(ArgumentError(
        "nested bounded formations are not supported"
    ))
    validate(part.pattern, bounded_members)
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
# course, or its own entry of a lay schedule with one path per course.
strand_path(path, course::Int, maximum_course::Int) = iszero(course) ? nothing : path
strand_path(::Nothing, course::Int, maximum_course::Int) = nothing
function strand_path(paths::Tuple, course::Int, maximum_course::Int)
    iszero(course) && return nothing
    length(paths) == maximum_course || throw(DimensionMismatch(
        "lay schedule requires exactly $maximum_course courses; got $(length(paths))"
    ))
    return paths[course]
end

"""
$(TYPEDSIGNATURES)

Construct the members and placed primitives of a bounded wire formation from its
geometric context, physical boundary, compaction selection, and ordered part
declarations. Coordinates and dimensions are in meters. Return `(boundary, members,
primitives)` in construction order. A disk formation declares its strand, after an
optional center wire. At the cable center, the courses surround the declared center wire,
or by default a copy of the course wire. After a part, they surround its outline, and the
formation has no center. Disk wires pack in circular courses, and a Rectangle strand bends
around the inner disk.
"""
function bounded_members(
        ::EmptyBoundary, boundary_shape::Disk, compact::Bool, strand_part::Group)
    return bounded_members(boundary_shape, strand_part, strand_part, compact)
end

function bounded_members(::EmptyBoundary, boundary_shape::Disk, compact::Bool,
        center_part::Group, strand_part::Group)
    center_part.pattern === nothing || throw(ArgumentError(
        "a center wire declares no course counts"
    ))
    return bounded_members(boundary_shape, center_part, strand_part, compact)
end

function bounded_members(::Disk, ::Disk, ::Bool, ::Group, ::Group)
    throw(ArgumentError(
        "a stranded formation after a part has no center of its own; omit center"
    ))
end

# Around a part, Disk wires stay natural. Compaction in a ring waits for a real case.
function bounded_members(
        context::Disk,
        boundary_shape::Disk,
        compact::Bool,
        strand_part::Group{<:Pose2, <:Region{<:Disk}}
)
    compact && throw(ArgumentError(
        "a stranded formation around a part takes natural wires; omit compact"
    ))
    return bounded_members(boundary_shape, context, NamedTuple[], strand_part, compact)
end

function bounded_members(
        context::Disk,
        boundary_shape::Disk,
        compact::Bool,
        strand_part::Group{<:Pose2, <:Region{<:Rectangle}}
)
    return bounded_members(boundary_shape, context, NamedTuple[], strand_part, compact)
end

function bounded_members(
        boundary_shape::Disk,
        center_part::Group{<:Pose2, <:Region{<:Disk}},
        strand_part::Group{<:Pose2, <:Region{<:Union{Disk, Rectangle}}},
        compact::Bool
)
    center = center_part.member.primitive
    center.r <= boundary_shape.r || throw(DomainError(
        center.r, "the center wire exceeds the circular core boundary"
    ))
    inside = NamedTuple[(
        source = center_part.member,
        course = 0,
        member = 1,
        path = nothing,
        angle = zero(boundary_shape.r),
        site = (boundary_shape.at.x, boundary_shape.at.y)
    )]
    return bounded_members(boundary_shape, center, inside, strand_part, compact)
end

# At the cable center, rectangular strands without a declared Disk center have no wire to
# bend around.
function bounded_members(
        boundary_shape::Disk,
        center_part::Group{<:Pose2, <:Region{<:Rectangle}},
        strand_part::Group{<:Pose2, <:Region{<:Rectangle}},
        compact::Bool
)
    throw(ArgumentError("a rectangular stranded core requires center=Disk(...)"))
end

# Circular courses of Disk wires surround the inner disk and follow the members inside it.
function bounded_members(
        boundary_shape::Disk,
        inner::Disk,
        inside::AbstractVector,
        strand_part::Group{<:Pose2, <:Region{<:Disk}},
        compact::Bool
)
    sites = circular_courses(boundary_shape, inner, strand_part.member.primitive, compact,
        strand_part.pattern)
    maximum_course = maximum(site.course for site in sites)
    members = NamedTuple[inside...]
    for placement in sites
        push!(members, (;
            source = strand_part.member,
            placement.course,
            member = length(members) + 1,
            path = strand_path(strand_part.path, placement.course, maximum_course),
            placement.angle,
            placement.site
        ))
    end
    primitives = compact ? deform_disk_members(boundary_shape, members) :
                 [resolve(Pose2(member.site..., member.angle), member.source.primitive)
                  for member in members]
    return boundary_shape, members, primitives
end

# Bent courses of a Rectangle strand surround the inner disk, after the members inside it at
# the formation center.
function bounded_members(
        boundary_shape::Disk,
        inner::Disk,
        inside::AbstractVector,
        strand_part::Group{<:Pose2, <:Region{<:Rectangle}},
        compact::Bool
)
    compact || throw(ArgumentError(
        "rectangular strands require area-preserving bending"
    ))
    strips = rectangular_strands(boundary_shape, inner, strand_part.member.primitive,
        strand_part.pattern)
    maximum_course = maximum(strip.course for strip in strips)
    members = NamedTuple[inside...]
    primitives = AbstractShape[
        resolve(Pose2(boundary_shape.at.x, boundary_shape.at.y, boundary_shape.at.φ),
            member.source.primitive) for member in inside
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
    # The outermost bent course ends on the disk boundary, the finished outline.
    return boundary_shape, members, primitives
end

# A non-circular boundary maps one center wire, a copy of the course wire, and the inferred
# or declared courses into its outline by swept area.
function bounded_members(::EmptyBoundary,
        boundary_shape::Union{SectorShape, Rectangle, Ellipse, Polygon}, ::Nothing,
        part::Group)
    courses, primitives, maximum_course = mapped_courses(boundary_shape,
        part.member.primitive, part.pattern)
    members = [(
                   source = part.member,
                   placement.course,
                   member,
                   path = strand_path(part.path, placement.course, maximum_course),
                   placement.angle,
                   placement.site
               ) for (member, placement) in enumerate(courses)]
    return boundary_shape, members, primitives
end

function bounded_members(::Union{EmptyBoundary, Disk},
        ::Union{Disk, SectorShape, Rectangle, Ellipse, Polygon}, compact, parts::Group...)
    throw(ArgumentError(
        "a stranded core declares one strand, after at most one center wire in a disk"
    ))
end

"""
$(TYPEDSIGNATURES)

Resolve a bounded wire group in its geometric context into `(boundary_shape, members,
primitives)`. The context is `EmptyBoundary()` at the cable center, or the circular outline
of the preceding part. The returned geometric boundary is the declared formation
boundary, which follows the group's pose. Coordinates and dimensions are in meters.
Member and primitive order follows construction. An empty formation raises
`ArgumentError`.
"""
function bounded_members(context::Union{EmptyBoundary, Disk}, group::Group)
    boundary_shape = resolve(group.at, group.boundary)
    outer, members, primitives = bounded_members(context, boundary_shape, group.compact,
        bounded_declarations(group)...)
    isempty(members) && throw(ArgumentError(
        "a bounded formation requires at least one strand"
    ))
    return outer, members, primitives
end

# The courses of a bounded strand are inferred, or declared by a positive base count (course
# k holds k times it) or by a nonempty schedule of positive counts, from the inside out.
validate(counts::Nothing, ::typeof(bounded_members)) = counts
function validate(counts::Integer, ::typeof(bounded_members))
    counts > 0 || throw(ArgumentError(
        "a base course count must be positive; got $counts"
    ))
    return counts
end
function validate(counts::Bool, ::typeof(bounded_members))
    throw(ArgumentError("a course count must be an integer; got $counts"))
end
function validate(counts::AbstractVector{<:Integer}, ::typeof(bounded_members))
    !isempty(counts) && all(>(0), counts) || throw(ArgumentError(
        "a course count schedule must be nonempty and positive; got $counts"
    ))
    return counts
end
function validate(counts, ::typeof(bounded_members))
    throw(ArgumentError(
        "a bounded strand inventory declares course counts, not a " *
        "$(nameof(typeof(counts))) or other placement"
    ))
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

function resolve_bounded(context::Union{EmptyBoundary, Disk}, group::Group)
    outer, members, primitives = bounded_members(context, group)
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
                    ),
                    axis = formation_center
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
# Holes only remove area. The outer primitive of a fill bounds its minimum radius from below.
_minimum_radius(primitive::DifferenceShape) = _minimum_radius(primitive.outer)
function _minimum_radius(primitive::AbstractShape)
    center = centroid(primitive)
    center_radius = hypot(center...)
    iszero(center_radius) && return zero(eltype(primitive))
    φ = atan(center[2], center[1])
    return -support(primitive, φ + pi)
end

# At the cable center, a bounded formation fills its whole geometric boundary. After a
# part, a formation in a disk fills the ring between the circular outline of that part and
# the disk.
function _resolve_group(
        context::EmptyBoundary,
        group::Group{<:Pose2, <:Any, <:Any, <:Any, <:Any,
            <:Union{Disk, Sector, Rectangle, Ellipse, Polygon}}
)
    return resolve_bounded(context, group)
end
function _resolve_group(
        context::Disk,
        group::Group{<:Pose2, <:Any, <:Any, <:Any, <:Any, <:Disk}
)
    return resolve_bounded(validate(context, group), group)
end

# The preceding circular outline of a bounded formation shares its center and lies inside
# its disk.
function validate(
        context::Disk, group::Group{<:Pose2, <:Any, <:Any, <:Any, <:Any, <:Disk})
    boundary_shape = resolve(group.at, group.boundary)
    isapprox(context.at.x, boundary_shape.at.x) &&
        isapprox(context.at.y, boundary_shape.at.y) && context.r < boundary_shape.r ||
        throw(DomainError(
            context,
            "the formation boundary must be concentric with and outside the preceding part"
        ))
    return context
end
function _resolve_group(
        ::AbstractShape,
        ::Group{<:Pose2, <:Any, <:Any, <:Any, <:Any,
            <:Union{Disk, Sector, Rectangle, Ellipse, Polygon}}
)
    throw(ArgumentError(
        "after a part, a stranded formation takes a Disk boundary around a circular " *
        "outline; a Sector, Rectangle, Ellipse or Polygon boundary starts at the cable " *
        "center"
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
            paths = group.path === nothing ? resolved.paths :
                    (resolved.paths...,
                (
                    path = group.path,
                    radius = _path_radius(
                        nothing,
                        Pose2(center[1], center[2], 0),
                        placed.primitive
                    ),
                    axis = (group.at.x, group.at.y)
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
            paths = group.path === nothing ? resolved.paths :
                    (resolved.paths...,
                (
                    path = group.path,
                    radius = _resolved_path_radius(
                        group.compact,
                        pattern,
                        pose,
                        local_primitive
                    ),
                    axis = (group.at.x, group.at.y)
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
