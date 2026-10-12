struct _ConnectionsAbsent end

const _CONNECTIONS_ABSENT = _ConnectionsAbsent()

_parameter_type(value) = typeof(value)
_parameter_type(grid::DeterministicGrid) = eltype(grid)
_parameter_type(::Union{RelativeGrid, AbsoluteGrid}) = Real
_parameter_type(::Gridspace{Target}) where {Target} = Target

at(; x = 0, y = 0, φ = 0, combine::Symbol = :product) = at(x, y; φ, combine)

at(pose::DataModel.Pose2) = pose

function _compose_member(outer, member::DataModel.AssemblyMember)
    DataModel.AssemblyMember(member.item, outer * member.at)
end
_compose_member(outer, pose::DataModel.Pose2) = outer * pose
function _compose_member(outer, placement::NamedTuple)
    merge(
        placement,
        (pose = outer * placement.pose,)
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
        ::Type{<:Real}, ::Type{<:Real},
        x, y, φ, ::_ConnectionsAbsent, combine::Symbol
)
    return parameterize(DataModel.Pose2, DataModel.Pose2, (x, y, φ); combine)
end

function _at_pair(
        ::Type{<:DataModel.AbstractCablePart}, ::Type{<:DataModel.Pose2},
        part, pose, φ, ::_ConnectionsAbsent, combine::Symbol
)
    iszero(φ) || throw(ArgumentError(
        "at(part, pose) does not accept a second rotation"
    ))
    return parameterize(
        DataModel.AssemblyMember, DataModel.AssemblyMember, (part, pose); combine
    )
end

function _at_pair(
        ::Type{<:DataModel.CableDesign}, ::Type{<:DataModel.Pose2},
        design, pose, φ, connections, combine::Symbol
)
    connections isa _ConnectionsAbsent && throw(ArgumentError(
        "at(design, pose) requires connections"
    ))
    iszero(φ) || throw(ArgumentError(
        "at(design, pose) does not accept a second rotation"
    ))
    # The placed record, as the formations write it.
    caller = (design, pose, connections) -> (; design, pose, connections)
    return parameterize(NamedTuple, caller, (design, pose, connections); combine)
end

function _at_pair(
        ::Type{<:Union{Tuple, AbstractVector}}, ::Type{<:DataModel.Pose2},
        placements, pose, φ, ::_ConnectionsAbsent,
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

- `x`: horizontal translation \\[m\\].
- `y`: vertical translation \\[m\\].
- `subject`: physical part, completed design or placement collection.
- `pose`: existing `Pose2` declaration.

# Keywords

- `φ=0`: counter-clockwise rotation \\[rad\\].
- `connections`: terminal-to-active-phase declaration required for a completed
  cable design. Use one-based active phase IDs and `0` for a grounded or eliminated
  conductor.
- `combine=:product`: gridspace composition rule.

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
    left_type = _parameter_type(left)
    right_type = _parameter_type(right)
    (left_type === Union{} || right_type === Union{}) && throw(ArgumentError(
        "at does not accept an empty finite source"
    ))
    return _at_pair(left_type, right_type, left, right, φ, connections, combine)
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
