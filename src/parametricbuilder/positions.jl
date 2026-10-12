"""
    at(x, y; φ = 0, combine = :product)
    at(; x = 0, y = 0, φ = 0, combine = :product)
    at(part, pose; combine = :product)
    at(design, pose; connections, combine = :product)
    at(placements, pose; combine = :product)
    at(subject, x, y; φ = 0, combine = :product, keywords...)

Construct a pose, place a local cable part, place a completed cable with its
connections, or compose an outer transform onto a placement collection. The last form
places `subject` at the pose `at(x, y; φ, combine)`, and passes the other keywords to the
placement.

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
- `combine=:product`: gridspace composition rule.

# Returns

- A `Pose2`, placed member, placed-design record, transformed collection, or
  the corresponding `Gridspace` when a direct argument is a `Grid` or `Gridspace`.
"""
function at end

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

at(pose::DataModel.Pose2) = pose

function at(
        part::Union{DataModel.AbstractCablePart, Gridspace{<:DataModel.AbstractCablePart},
            DeterministicGrid{<:Tuple{DataModel.AbstractCablePart,
                Vararg{DataModel.AbstractCablePart}}}},
        pose::Union{DataModel.Pose2, Gridspace{<:DataModel.Pose2},
            DeterministicGrid{<:Tuple{DataModel.Pose2, Vararg{DataModel.Pose2}}}};
        combine::Symbol = :product
)
    return parameterize(
        DataModel.AssemblyMember, DataModel.AssemblyMember, (part, pose); combine
    )
end

function at(
        design::Union{DataModel.CableDesign, Gridspace{<:DataModel.CableDesign},
            DeterministicGrid{<:Tuple{DataModel.CableDesign, Vararg{DataModel.CableDesign}}}},
        pose::Union{DataModel.Pose2, Gridspace{<:DataModel.Pose2},
            DeterministicGrid{<:Tuple{DataModel.Pose2, Vararg{DataModel.Pose2}}}};
        connections,
        combine::Symbol = :product
)
    # The placed record, as the formations write it.
    caller = (design, pose, connections) -> (; design, pose, connections)
    return parameterize(NamedTuple, caller, (design, pose, connections); combine)
end

function at(
        placements::Union{Tuple, Gridspace{<:Tuple},
            DeterministicGrid{<:Tuple{Tuple, Vararg{Tuple}}}},
        pose::Union{DataModel.Pose2, Gridspace{<:DataModel.Pose2},
            DeterministicGrid{<:Tuple{DataModel.Pose2, Vararg{DataModel.Pose2}}}};
        combine::Symbol = :product
)
    return parameterize(Tuple, _compose_collection, (placements, pose); combine)
end

function at(
        placements::Union{AbstractVector, Gridspace{<:AbstractVector},
            DeterministicGrid{<:Tuple{AbstractVector, Vararg{AbstractVector}}}},
        pose::Union{DataModel.Pose2, Gridspace{<:DataModel.Pose2},
            DeterministicGrid{<:Tuple{DataModel.Pose2, Vararg{DataModel.Pose2}}}};
        combine::Symbol = :product
)
    return parameterize(Vector, _compose_collection, (placements, pose); combine)
end

function at(subject, x, y; φ = 0, combine::Symbol = :product, keywords...)
    return at(subject, at(x, y; φ, combine); combine, keywords...)
end

_compose_collection(placements, pose) = map(member -> _compose_member(pose, member), placements)

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
