"""
$(TYPEDEF)

Represent one region as a cable or a system places it: its declaration and its resolved
cross-section.

$(TYPEDFIELDS)
"""
struct PlacedRegion{R <: Region, S <: AbstractShape}
    "Declared region."
    source::R
    "Resolved cross-section in the frame of the cable or the system that places the region."
    shape::S
    function PlacedRegion(source::R, shape::S) where {R <: Region, S <: AbstractShape}
        return new{R, S}(source, shape)
    end
end

"""
$(TYPEDEF)

Describe one cable by its declared parts and the regions that they resolve to, with tables
that map the regions to terminals, sub-cables, repeats, lays, compactions and fills. Each
table row holds resolved values in the floating type `T`, to which numbers from shapes and
materials promote. A row exists only for a fact that holds. The declarations stay in `tree`.
The resolved shapes of all regions are in `T`.

$(TYPEDFIELDS)
"""
struct CableDesign{
    T <: Real,
    O <: AbstractCablePart,
    V <: AbstractVector{<:PlacedRegion},
    B <: AbstractShape
}
    "Name of the cable."
    id::String
    "Declared parts, as given to [`build`](@ref)."
    tree::O
    "Placed regions, from the inside out. Together they cover the cross-section."
    regions::V
    "Outline of the cable."
    outer::B
    """
    Terminals, each as `(name, scope, regions)`, in the order of their first region. A row
    gives the declared name and the regions of the terminal. Its scope lists the rows of the
    sub-cables that contain the terminal, from the outermost. It is empty for a terminal of
    the cable's own scope.
    """
    terminals::Vector{@NamedTuple{
        name::Symbol, scope::Vector{Int},
        regions::Vector{Int}}}
    """
    Sub-cables, each as `(at, regions, members)`. A row gives the placement of the sub-cable
    in the cable, its own regions and the rows of the sub-cables inside it. Row 1 is the
    cable, and each copy of a repeated part other than a region is a sub-cable. Each other
    row is a member of one earlier row, and each region belongs to one row.
    """
    assemblies::Vector{@NamedTuple{
        at::Placement{T}, regions::Vector{Int}, members::Vector{Int}}}
    "Repeats, each as `(placements, copies)`: the placement and the regions of each copy."
    repeats::Vector{@NamedTuple{
        placements::Vector{Placement{T}}, copies::Vector{UnitRange{Int}}}}
    """
    Lays, each as `(regions, radii, pitches, overlengths, axis)`. A row gives the regions of
    each laid copy and the point that the copies wind about \\[m\\]. For each copy, it gives
    the lay radius \\[m\\], the pitch signed by the handedness \\[m\\] and the overlength
    \\[dimensionless\\].
    """
    lays::Vector{@NamedTuple{
        regions::Vector{UnitRange{Int}}, radii::Vector{T}, pitches::Vector{T}, overlengths::Vector{T},
        axis::Tuple{T, T}}}
    """
    Compactions, each as `(repeat, factor)`: the row of a compacted repeat and its area
    factor \\[dimensionless\\].
    """
    compactions::Vector{@NamedTuple{repeat::Int, factor::T}}
    """
    Fills, each as `(region, contents)`: the fill region of an enclosure and the regions
    inside it. A fill whose conductive contents all belong to one terminal is interstitial
    to that conductor.
    """
    fills::Vector{@NamedTuple{region::Int, contents::Vector{Int}}}
    function CableDesign{T}(id::String, tree::O, regions::V, outer::B,
            terminals::AbstractVector{<:NamedTuple},
            assemblies::AbstractVector{<:NamedTuple},
            repeats::AbstractVector{<:NamedTuple}, lays::AbstractVector{<:NamedTuple},
            compactions::AbstractVector{<:NamedTuple},
            fills::AbstractVector{<:NamedTuple}) where {
            T <: Real, O <: AbstractCablePart, V <: AbstractVector{<:PlacedRegion},
            B <: AbstractShape}
        return validate(new{T, O, V, B}(id, tree, regions, outer, terminals, assemblies,
            repeats, lays, compactions, fills))
    end
end

Base.eltype(::CableDesign{T}) where {T} = T

# The regions cover the outline of the cable, each conductive region has a terminal, and the
# sub-cables nest with each region in one of them.
function validate(design::CableDesign)
    covered = sum(placed -> area(placed.shape), design.regions)
    outline = area(design.outer)
    isapprox(covered, outline; rtol = sqrt(eps(float(outline)))) ||
        throw(DomainError(covered,
            "the regions of the cable $(design.id) must cover its outline of " *
            "$outline m²: a gap or an overlap lies between its parts"))
    conductive = [index
                  for (index, placed) in pairs(design.regions)
                  if placed.source.material isa Material{<:Real, :conductor}]
    assigned = [index for row in design.terminals for index in row.regions]
    loose = unique([design.regions[index].source.tag
                    for index in setdiff(conductive, assigned)])
    isempty(loose) || throw(ArgumentError(
        "the conductive regions $loose of the cable $(design.id) must lie inside a " *
        "Terminal, which names their terminal"))
    rows = design.assemblies
    members = sort!([member for row in rows for member in row.members])
    !isempty(rows) && members == 2:length(rows) &&
    all(k -> all(>(k), rows[k].members), eachindex(rows)) ||
        throw(ArgumentError(
            "the sub-cables of the cable $(design.id) must nest: row 1 of its assemblies " *
            "is the cable, and each other row is a member of exactly one earlier row"))
    sort!([index for row in rows for index in row.regions]) == eachindex(design.regions) ||
        throw(ArgumentError(
            "each region of the cable $(design.id) must lie in exactly one row of its " *
            "assemblies"))
    return design
end

"""
$(TYPEDSIGNATURES)

Build the cable `id` from its parts. The parts resolve from the center of the cable outward
into placed regions. Each node of the declaration records its facts in the tables of the
design, such as the conductive regions of a terminal or the copies of a repeat. The cable is
the first sub-cable of its design, at the origin. Each copy of a repeated part other than a
region is a sub-cable inside it.

# Arguments

- `id`: name of the cable.
- `parts`: declared parts, from the inside out.

# Returns

- A [`CableDesign`](@ref).

# Errors

- Throws `ArgumentError` when no part is given or when a conductive region lies outside
  every [`Terminal`](@ref).
- Throws `DomainError` when a part overlaps the part inside it or leaves a gap around it, or
  when a region leaves its enclosure.
"""
function build(::Type{CableDesign}, id::AbstractString, parts::AbstractCablePart...)
    tree = Stack(parts...)
    described = resolve(EmptyBoundary(), tree)
    T = mapreduce(
        placed -> promote_type(eltype(placed.shape), eltype(placed.source.material)),
        promote_type, described.regions)
    regions = [PlacedRegion(placed.source, map(x -> convert(T, x), placed.shape))
               for placed in described.regions]
    cable = (
        at = Placement(zero(T), zero(T)),
        regions = Int[index
                      for index in eachindex(regions)
                      if all(row -> index ∉ row.regions, described.assemblies)],
        members = Int[index + 1
                      for index in eachindex(described.assemblies)
                      if all(row -> index ∉ row.members, described.assemblies)])
    assemblies = NamedTuple[cable;
                            NamedTuple[(at = row.at, regions = row.regions,
                                           members = row.members .+ 1)
                                       for row in described.assemblies]]
    terminals = NamedTuple[(name = row.name, scope = row.scope .+ 1, regions = row.regions)
                           for row in described.terminals]
    return CableDesign{T}(
        String(id), tree, regions, map(x -> convert(T, x), described.outer),
        terminals, assemblies, described.repeats, described.lays, NamedTuple[],
        described.fills)
end

function Base.map(f, design::CableDesign)
    return build(CableDesign, design.id, map.(f, design.tree.items)...)
end
