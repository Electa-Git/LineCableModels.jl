"""
    at(x, y; φ = 0, combine = :product)
    at(; x = 0, y = 0, φ = 0, combine = :product)
    at(pose)
    at(part, pose; combine = :product)
    at(part, x, y; φ = 0, combine = :product)
    at(design, pose; connections, combine = :product)
    at(design, x, y; φ = 0, connections, combine = :product)
    at(placements, pose; combine = :product)
    at(placements, x, y; φ = 0, combine = :product)

Construct a pose, place a local cable part, place a completed cable with its
connections, or compose an outer transform onto a placement collection. A placement at
`x` and `y` lifts over them, and over `φ`, as direct parameters of the one call, and it
builds its pose from the scalars of each point.

# Arguments

- `x`: horizontal translation \\[m\\].
- `y`: vertical translation \\[m\\].
- `part`: physical cable part, placed as a member of an assembly.
- `design`: completed cable design, placed in a line system.
- `placements`: tuple or vector of members, poses or placed-design records.
- `pose`: existing `Pose2` declaration.

# Keywords

- `φ=0`: counter-clockwise rotation \\[rad\\] of a pose built from `x` and `y`.
- `connections`: terminal-to-active-phase declaration of a completed cable design. Use
  one-based active phase IDs and `0` for a grounded or eliminated conductor.
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- A `Pose2`, placed member, placed-design record, transformed collection, or
  the corresponding `Gridspace` when an argument or a keyword is a `Grid` or `Gridspace`.
"""
function at end

# One value of `T`, or a finite space of them: the subjects and the poses that `at` places.
const Placed{T} = Union{T, Gridspace{<:T}, DeterministicGrid{<:Tuple{T, Vararg{T}}}}

function at(
        x::Union{Real, AbstractUncertainGrid, DeterministicGrid{<:Tuple{Vararg{Real}}},
            Gridspace{<:Real}},
        y::Union{Real, AbstractUncertainGrid, DeterministicGrid{<:Tuple{Vararg{Real}}},
            Gridspace{<:Real}};
        φ = 0,
        combine::Symbol = :product
)
    return DataModel.Pose2(x, y, φ; combine)
end

at(; x = 0, y = 0, φ = 0, combine::Symbol = :product) = at(x, y; φ, combine)

at(pose::Placed{DataModel.Pose2}) = pose

function at(part::Placed{DataModel.AbstractCablePart}, pose::Placed{DataModel.Pose2};
        combine::Symbol = :product)
    return parameterize(
        DataModel.AssemblyMember, DataModel.AssemblyMember, (part, pose); combine
    )
end
function at(part::Placed{DataModel.AbstractCablePart}, x, y; φ = 0,
        combine::Symbol = :product)
    caller = (part, x, y, φ) -> DataModel.AssemblyMember(part, DataModel.Pose2(x, y, φ))
    return parameterize(DataModel.AssemblyMember, caller, (part, x, y, φ); combine)
end

function at(design::Placed{DataModel.CableDesign}, pose::Placed{DataModel.Pose2};
        connections, combine::Symbol = :product)
    # The placed record, as the formations write it.
    caller = (design, pose, connections) -> (; design, pose, connections)
    return parameterize(NamedTuple, caller, (design, pose, connections); combine)
end
function at(design::Placed{DataModel.CableDesign}, x, y; φ = 0, connections,
        combine::Symbol = :product)
    caller = (design, x, y, φ, connections) ->
        (; design, pose = DataModel.Pose2(x, y, φ), connections)
    return parameterize(NamedTuple, caller, (design, x, y, φ, connections); combine)
end

function at(placements::Placed{Tuple}, pose::Placed{DataModel.Pose2};
        combine::Symbol = :product)
    return parameterize(Tuple, _compose_collection, (placements, pose); combine)
end
function at(placements::Placed{Tuple}, x, y; φ = 0, combine::Symbol = :product)
    caller = (placements, x, y, φ) ->
        _compose_collection(placements, DataModel.Pose2(x, y, φ))
    return parameterize(Tuple, caller, (placements, x, y, φ); combine)
end

function at(placements::Placed{AbstractVector}, pose::Placed{DataModel.Pose2};
        combine::Symbol = :product)
    return parameterize(Vector, _compose_collection, (placements, pose); combine)
end
function at(placements::Placed{AbstractVector}, x, y; φ = 0, combine::Symbol = :product)
    caller = (placements, x, y, φ) ->
        _compose_collection(placements, DataModel.Pose2(x, y, φ))
    return parameterize(Vector, caller, (placements, x, y, φ); combine)
end

_compose_collection(placements, pose) = map(member -> _compose_member(pose, member), placements)

function _compose_member(outer, member::DataModel.AssemblyMember)
    DataModel.AssemblyMember(member.part, outer * member.at)
end
_compose_member(outer, pose::DataModel.Pose2) = outer * pose
function _compose_member(outer, placement::NamedTuple)
    merge(
        placement,
        (pose = outer * placement.pose,)
    )
end
