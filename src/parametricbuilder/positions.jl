_pose(x, y, φ) = DataModel.Pose2(x, y, φ)

struct _ConnectionsAbsent end

const _CONNECTIONS_ABSENT = _ConnectionsAbsent()

_parameter_type(value) = typeof(value)
_parameter_type(::DeterministicGrid{Values}) where {Values} = eltype(Values)
_parameter_type(::Union{RelativeGrid, AbsoluteGrid}) = Real
_parameter_type(::Gridspace{Target}) where {Target} = Target

function _at_kind(left, right, connections)
    left_type = _parameter_type(left)
    right_type = _parameter_type(right)
    (left_type === Union{} || right_type === Union{}) && throw(ArgumentError(
        "at does not accept an empty finite source"
    ))
    return _at_kind(left_type, right_type, connections)
end

_at_kind(::Type{<:Real}, ::Type{<:Real}, ::_ConnectionsAbsent) = Val(:pose)
function _at_kind(
        ::Type{<:DataModel.AbstractCablePart},
        ::Type{<:DataModel.Pose2},
        ::_ConnectionsAbsent
)
    Val(:member)
end
function _at_kind(
        ::Type{<:DataModel.CableDesign},
        ::Type{<:DataModel.Pose2},
        connections
)
    Val(:design)
end
function _at_kind(
        ::Type{<:Union{Tuple, AbstractVector}},
        ::Type{<:DataModel.Pose2},
        ::_ConnectionsAbsent
)
    Val(:collection)
end

at(; x = 0, y = 0, φ = 0, combine::Symbol = :product) = at(x, y; φ, combine)

at(pose::DataModel.Pose2) = pose

_placed_member(part, pose) = DataModel.AssemblyMember(part, pose)

function _placed_design(design, pose, connections)
    (
        design = design,
        pose = pose,
        connections = connections
    )
end

_compose_pose(outer::DataModel.Pose2, inner::DataModel.Pose2) = outer * inner
function _compose_member(outer, member::DataModel.AssemblyMember)
    DataModel.AssemblyMember(member.item, _compose_pose(outer, member.at))
end
_compose_member(outer, pose::DataModel.Pose2) = _compose_pose(outer, pose)
function _compose_member(outer, placement::NamedTuple)
    merge(
        placement,
        (pose = _compose_pose(outer, placement.pose),)
    )
end

function _compose_collection(placements, pose)
    placements isa Union{Tuple, AbstractVector} || throw(ArgumentError(
        "outer placement requires a tuple or vector"
    ))
    return map(member -> _compose_member(pose, member), placements)
end

_collection_target(::Tuple) = Tuple
_collection_target(::AbstractVector) = Vector
_collection_target(::Gridspace{Target}) where {Target} = Target <: Tuple ? Tuple : Vector

function _at_pair(
        ::Val{:pose}, x, y, φ, ::_ConnectionsAbsent, combine::Symbol
)
    return parameterize(DataModel.Pose2, _pose, (x, y, φ); combine)
end

function _at_pair(
        ::Val{:member}, part, pose, φ, ::_ConnectionsAbsent, combine::Symbol
)
    iszero(φ) || throw(ArgumentError(
        "at(part, pose) does not accept a second rotation"
    ))
    return parameterize(
        DataModel.AssemblyMember, _placed_member, (part, pose); combine
    )
end

function _at_pair(
        ::Val{:design}, design, pose, φ, connections, combine::Symbol
)
    connections isa _ConnectionsAbsent && throw(ArgumentError(
        "at(design, pose) requires connections"
    ))
    iszero(φ) || throw(ArgumentError(
        "at(design, pose) does not accept a second rotation"
    ))
    return parameterize(
        NamedTuple, _placed_design, (design, pose, connections); combine
    )
end

function _at_pair(
        ::Val{:collection}, placements, pose, φ, ::_ConnectionsAbsent,
        combine::Symbol
)
    iszero(φ) || throw(ArgumentError(
        "at(placements, pose) does not accept a second rotation"
    ))
    return parameterize(
        _collection_target(placements),
        _compose_collection,
        (placements, pose);
        combine
    )
end

"""
$(TYPEDSIGNATURES)

Construct a pose, place a local cable part, place a completed cable with its
connections, or compose an outer transform onto a placement collection.

# Arguments

- `x`: Horizontal translation \\[m\\].
- `y`: Vertical translation \\[m\\].
- `subject`: Physical part, completed design, or placement collection.
- `pose`: Existing `Pose2` declaration.

# Keywords

- `φ=0`: Counter-clockwise rotation \\[rad\\].
- `connections`: Terminal connection declaration required for a completed
  cable design.
- `combine=:product`: Gridspace composition rule.

# Returns

- A `Pose2`, placed member, placed-design record, transformed collection, or
  the corresponding `Gridspace` when a direct argument varies.
"""
function at(
        left,
        right;
        φ = 0,
        connections = _CONNECTIONS_ABSENT,
        combine::Symbol = :product
)
    kind = _at_kind(left, right, connections)
    return _at_pair(kind, left, right, φ, connections, combine)
end

function at(
        subject,
        x,
        y;
        φ = 0,
        connections = _CONNECTIONS_ABSENT,
        combine::Symbol = :product
)
    pose = at(x, y; φ, combine)
    return at(subject, pose; connections, combine)
end

function _connection_for_member(connections::NamedTuple, member::Int, count::Int)
    resolved = map(Base.values(connections)) do value
        if value isa Union{Tuple, AbstractVector}
            length(value) == count || throw(DimensionMismatch(
                "connection schedules must contain $count values"
            ))
            return value[member]
        end
        return value
    end
    return NamedTuple{keys(connections)}(resolved)
end

function _connection_for_member(connections::AbstractDict, member::Int, count::Int)
    return Dict(key => begin
                    if value isa Union{Tuple, AbstractVector}
                        length(value) == count || throw(DimensionMismatch(
                            "connection schedules must contain $count values"
                        ))
                        value[member]
                    else
                        value
                    end
                end for (key, value) in connections)
end

function _connection_for_member(connections::Union{Tuple, AbstractVector}, member, count)
    length(connections) == count || throw(DimensionMismatch(
        "formation connections must contain $count declarations"
    ))
    return connections[member]
end

function _formation(design, local_poses, center, connections)
    count = length(local_poses)
    return [_placed_design(
                design,
                center * pose,
                _connection_for_member(connections, member, count)
            )
            for (member, pose) in enumerate(local_poses)]
end

function _trefoil(design, center, spacing, connections, φ0)
    spacing > zero(spacing) || throw(DomainError(
        spacing, "trefoil spacing must be positive"
    ))
    coordinates = DataModel.trefoil_formation(0, 0, spacing / 2)
    poses = DataModel.Pose2[DataModel.Pose2(coordinates[index], coordinates[index + 1], 0)
                            for index in 1:2:6]
    origin = center * DataModel.Pose2(0, 0, φ0)
    return _formation(design, poses, origin, connections)
end

"""
$(TYPEDSIGNATURES)

Place three copies of one completed cable in an equilateral trefoil formation.

# Arguments

- `design`: Completed cable design reused by all three members.

# Keywords

- `center=at(0, 0)`: Formation-centre pose \\[m, m, rad\\].
- `spacing`: Cable centre-to-centre distance \\[m\\].
- `connections`: Scalar or three-member terminal connection schedules.
- `φ0=0`: Formation rotation \\[rad\\].
- `combine=:product`: Gridspace composition rule.

# Returns

- Three placed-design records, or a `Gridspace{Vector}` when a direct argument
  varies.
"""
function trefoil(
        design;
        center = at(0, 0),
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
    offsets = (-spacing, zero(spacing), spacing)
    poses = DataModel.Pose2[vertical ? DataModel.Pose2(0, -offset, 0) :
                            DataModel.Pose2(offset, 0, 0)
                            for offset in offsets]
    return _formation(design, poses, center, connections)
end

"""
$(TYPEDSIGNATURES)

Place three copies of one completed cable in a horizontal flat formation.

`spacing` is the adjacent centre-to-centre distance \\[m\\]. Scalar connection
entries apply to every cable; three-element entries distribute by member.
"""
function hflat(
        design;
        center = at(0, 0),
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

`spacing` is the adjacent centre-to-centre distance \\[m\\]. Scalar connection
entries apply to every cable; three-element entries distribute by member.
"""
function vflat(
        design;
        center = at(0, 0),
        spacing,
        connections,
        combine::Symbol = :product
)
    caller = (resolved...) -> _flat(resolved..., true)
    return parameterize(
        Vector, caller, (design, center, spacing, connections); combine
    )
end
