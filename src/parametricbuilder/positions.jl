struct AtBuilder <: Builder{DataModel.Pose2} end

"""
```julia
at(x, y; φ = 0, combine = :product)
at(; x = 0, y = 0, φ = 0, combine = :product)
at(pose)
at(part, pose; combine = :product)
at(part, x, y; φ = 0, combine = :product)
at(design, pose; connections, combine = :product)
at(design, x, y; φ = 0, connections, combine = :product)
at(placements, pose; combine = :product)
at(placements, x, y; φ = 0, combine = :product)
```

Construct a pose, place a local cable part, place a completed cable with its
connections, or compose an outer transform onto a placement collection. A placement at
`x` and `y` builds its pose from them and `φ`.

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

- A `Pose2`, a placed member, a placed-design record or a transformed collection. When an
  argument or a keyword varies, a `Gridspace` of the type that its subject places.
"""
const at = AtBuilder()

# What `at` places follows its subject. It places a pose by default. A part becomes a member,
# a design a placed record, and a collection a collection of the same type.
target(::AtBuilder, ::Type{<:DataModel.AbstractCablePart}) = DataModel.AssemblyMember
target(::AtBuilder, ::Type{<:DataModel.CableDesign}) = NamedTuple
target(::AtBuilder, ::Type{<:Tuple}) = Tuple
target(::AtBuilder, ::Type{<:AbstractVector}) = Vector

build(::AtBuilder, x::Real, y::Real; φ = 0) = DataModel.Pose2(x, y, φ)
build(::AtBuilder; x = 0, y = 0, φ = 0) = DataModel.Pose2(x, y, φ)
build(::AtBuilder, pose::DataModel.Pose2) = pose
build(::AtBuilder, part::DataModel.AbstractCablePart, pose::DataModel.Pose2) =
    DataModel.AssemblyMember(part, pose)
# The placed record, as the formations write it.
build(::AtBuilder, design::DataModel.CableDesign, pose::DataModel.Pose2; connections) =
    (; design, pose, connections)
build(::AtBuilder, placements::Union{Tuple, AbstractVector}, pose::DataModel.Pose2) =
    map(member -> _compose_member(pose, member), placements)
# A subject placed at `x` and `y` is placed at the pose they declare.
build(builder::AtBuilder, subject, x::Real, y::Real; φ = 0, keywords...) =
    build(builder, subject, DataModel.Pose2(x, y, φ); keywords...)

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
