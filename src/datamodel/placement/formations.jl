# One member's entry of a connection schedule. A per-member schedule gives the member's
# value. Any other value serves every member.
function _schedule_entry(schedule::Union{Tuple, AbstractVector}, member::Int, count::Int)
    length(schedule) == count || throw(DimensionMismatch(
        "connection schedules must contain $count values"
    ))
    return schedule[member]
end
_schedule_entry(value, ::Int, ::Int) = value

function _connection_for_member(connections::NamedTuple, member::Int, count::Int)
    resolved = map(value -> _schedule_entry(value, member, count), Base.values(connections))
    return NamedTuple{keys(connections)}(resolved)
end

function _connection_for_member(connections::AbstractDict, member::Int, count::Int)
    return Dict(key => _schedule_entry(value, member, count) for (key, value) in connections)
end

function _connection_for_member(connections::Union{Tuple, AbstractVector}, member, count)
    length(connections) == count || throw(DimensionMismatch(
        "formation connections must contain $count declarations"
    ))
    return connections[member]
end

# The placed `(design, pose, connections)` records of a formation, in member order.
function _formation(design, local_poses, center, connections)
    count = length(local_poses)
    return [(
                design = design,
                pose = center * pose,
                connections = _connection_for_member(connections, member, count)
            )
            for (member, pose) in enumerate(local_poses)]
end

# A spacing below the outer diameter plus 1 μm is raised to it, with the outer radius of the
# design that this realization receives. No formation overlaps.
function _trefoil(design, center, spacing, connections, φ0)
    spacing > zero(spacing) || throw(DomainError(
        spacing, "trefoil spacing must be positive"
    ))
    side = max(spacing, 2outer_radius(design) + _CABLE_CLEARANCE)
    # The centers lie on a circle through the triangle's vertices, at 90°, 210° and 330°.
    distance = side / 2 / cosd(30)
    poses = Pose2[Pose2(zero(distance), distance, 0),
        Pose2(distance * cosd(210), distance * sind(210), 0),
        Pose2(distance * cosd(330), distance * sind(330), 0)]
    origin = center * Pose2(0, 0, φ0)
    return _formation(design, poses, origin, connections)
end

"""
$(TYPEDSIGNATURES)

Place three copies of one completed cable in an equilateral trefoil formation.

# Arguments

- `design`: completed cable design reused by all three members.

# Keywords

- `center=Pose2(0, 0)`: formation-center pose \\[m, m, rad\\].
- `spacing`: cable center-to-center distance \\[m\\]. A spacing below the outer diameter
  plus 1 μm is raised to it for each realization of `design`. The cables then touch,
  1 μm apart, and never overlap.
- `connections`: scalar or three-member terminal connection schedules.
- `φ0=0`: formation rotation \\[rad\\].
- `combine=:product`: gridspace composition rule.

# Returns

- Three placed-design records, or a `Gridspace{Vector}` when a direct argument
  varies.
"""
function trefoil(
        design;
        center = Pose2(0, 0),
        spacing,
        connections,
        φ0 = 0,
        combine::Symbol = :product
)
    values = (design, center, spacing, connections, φ0)
    return parameterize(Vector, _trefoil, values; combine)
end

function _flat(design, center, spacing, connections, vertical::Bool)
    spacing > zero(spacing) || throw(DomainError(
        spacing, "flat-formation spacing must be positive"
    ))
    adjacent = max(spacing, 2outer_radius(design) + _CABLE_CLEARANCE)
    offsets = (-adjacent, zero(adjacent), adjacent)
    poses = Pose2[vertical ? Pose2(0, -offset, 0) : Pose2(offset, 0, 0)
                  for offset in offsets]
    return _formation(design, poses, center, connections)
end

"""
$(TYPEDSIGNATURES)

Place three copies of one completed cable in a horizontal flat formation.

`spacing` is the adjacent center-to-center distance \\[m\\]. A spacing below the outer
diameter plus 1 μm is raised to it for each realization of `design`. The cables then
touch, 1 μm apart, and never overlap. Scalar connection entries apply to every cable.
Three-element entries distribute by member.
"""
function hflat(
        design;
        center = Pose2(0, 0),
        spacing,
        connections,
        combine::Symbol = :product
)
    caller = (resolved...) -> _flat(resolved..., false)
    return parameterize(
        Vector, caller, (design, center, spacing, connections); combine
    )
end

"""
$(TYPEDSIGNATURES)

Place three copies of one completed cable in a vertical flat formation.

`spacing` is the adjacent center-to-center distance \\[m\\]. A spacing below the outer
diameter plus 1 μm is raised to it for each realization of `design`. The cables then
touch, 1 μm apart, and never overlap. Scalar connection entries apply to every cable.
Three-element entries distribute by member.
"""
function vflat(
        design;
        center = Pose2(0, 0),
        spacing,
        connections,
        combine::Symbol = :product
)
    caller = (resolved...) -> _flat(resolved..., true)
    return parameterize(
        Vector, caller, (design, center, spacing, connections); combine
    )
end
