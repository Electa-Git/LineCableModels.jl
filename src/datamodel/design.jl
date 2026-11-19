"""
$(TYPEDEF)

Describe one cable by its resolved parts and the regions that they place. The resolved
parts stay in `tree`, where a copy of a repeat and a member of an assembly lie in their own
frames. Its tables map the regions to terminals and to the copies of repeats. Each table
row holds resolved values in the floating type `T`, to which numbers from shapes and
materials promote. The resolved shapes of all regions are in `T`.

$(TYPEDFIELDS)
"""
struct CableDesign{
    T <: Real,
    O <: AbstractCablePart,
    V <: AbstractVector{<:CablePart},
    B <: AbstractShape
}
    "Name of the cable."
    id::String
    "Resolved parts: each part given to [`build`](@ref), with its geometry."
    tree::O
    "Placed regions, from the inside out. Together they cover the cross-section."
    regions::V
    "Outline of the cable."
    outer::B
    """
    Terminals, each as `(name, scope, regions)`, in the order of their first region. A row
    gives the declared name and the regions of the terminal. Its scope lists the sub-cables
    that contain the terminal, from the outermost, each by its number in the order of the
    parts. The cable is number 1, and the scope of a terminal of the cable itself is empty.
    """
    terminals::Vector{@NamedTuple{
        name::Symbol, scope::Vector{Int},
        regions::Vector{Int}}}
    "Repeats, each as `(placements, copies)`: the placement and the regions of each copy."
    repeats::Vector{@NamedTuple{
        placements::Vector{Placement{T}}, copies::Vector{UnitRange{Int}}}}
    function CableDesign{T}(id::String, tree::O, regions::V, outer::B,
            terminals::AbstractVector{<:NamedTuple},
            repeats::AbstractVector{<:NamedTuple}) where {
            T <: Real, O <: AbstractCablePart, V <: AbstractVector{<:CablePart},
            B <: AbstractShape}
        return validate(new{T, O, V, B}(id, tree, regions, outer, terminals, repeats))
    end
end

Base.eltype(::CableDesign{T}) where {T} = T

# The regions cover the outline of the cable, and each conductive region has a terminal.
function validate(design::CableDesign)
    covered = sum(placed -> area(placed.shape), design.regions)
    outline = area(design.outer)
    isapprox(covered, outline; rtol = sqrt(eps(float(outline)))) ||
        throw(DomainError(covered,
            "the regions of the cable $(design.id) must cover its outline of " *
            "$outline m²: a gap or an overlap lies between its parts"))
    assigned = falses(length(design.regions))
    for row in design.terminals
        assigned[row.regions] .= true
    end
    loose = unique([design.regions[index].source.tag
                    for index in findall(!, assigned)
                    if design.regions[index].source.material isa CONDUCTOR])
    isempty(loose) || throw(ArgumentError(
        "the conductive regions $loose of the cable $(design.id) must lie inside a " *
        "Terminal, which names their terminal"))
    return design
end

"""
$(TYPEDSIGNATURES)

Build the cable `id` from its parts. The parts resolve from the center of the cable outward
into placed regions. The design keeps the resolved parts, and its tables list the
conductive regions of each terminal and the copies of each repeat. The cable is the first
sub-cable of its design, at the origin. Each member of an assembly and each copy of a
repeated part other than a region is a sub-cable inside it.

# Arguments

- `id`: name of the cable.
- `parts`: declared parts, from the inside out.

# Returns

- A [`CableDesign`](@ref).

# Errors

- Throws `ArgumentError` when no part is given or when a conductive region lies outside
  every [`Terminal`](@ref).
- Throws `ArgumentError` when a [`Group`](@ref) holds terminals of different names.
- Throws `DomainError` when a part overlaps the part inside it or leaves a gap around it, or
  when a region leaves its group or its enclosure.
"""
function build(::Type{CableDesign}, id::AbstractString, parts::AbstractCablePart...)
    regions = CablePart[]
    resolved = resolve(EmptyBoundary(), Stack(parts...), regions)
    return CableDesign(String(id), resolved.first, resolved.second, regions)
end

"""
$(TYPEDSIGNATURES)

Describe the cable `id` by its resolved parts `tree`, its outline `outer` and its placed
regions `regions`, from the inside out. The rows of the tables come from one walk of `tree`.
The cable is the first sub-cable.
"""
function CableDesign(id::String, tree::Stack, outer::AbstractShape,
        regions::AbstractVector{<:CablePart})
    rows = resolve(
        (placed = regions, regions = sizehint!(CablePart[], length(regions)),
            terminals = NamedTuple[], scopes = [1], repeats = NamedTuple[]),
        [1], tree)
    T = mapreduce(
        placed -> promote_type(eltype(placed.shape), eltype(placed.source.material)),
        promote_type, regions)
    return CableDesign{T}(id, tree,
        [CablePart(T, placed) for placed in regions],
        map(x -> convert(T, x), outer),
        rows.terminals, rows.repeats)
end

function Base.map(f, design::CableDesign)
    return CableDesign(design.id, map(f, design.tree), map(f, design.outer),
        [map(f, placed) for placed in design.regions])
end
