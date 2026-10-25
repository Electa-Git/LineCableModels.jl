"""
$(TYPEDEF)

Describe the cables of a line in one cross-section: the design of each cable, its placement
and the conductor index of each of its terminals. A conductor index `k ≥ 1` connects a
terminal to the conductor `k` of the circuit, terminals with one index form one conductor,
and the index zero eliminates a terminal. The system places the regions of each cable in its
own frame. Cables can touch, and the system keeps their placements as declared.

$(TYPEDFIELDS)
"""
struct CableSystem{
    T <: Real,
    D <: AbstractVector{<:CableDesign},
    V <: AbstractVector{<:PlacedRegion}
}
    "Name of the system."
    id::String
    "Length of the line \\[m\\]."
    line_length::T
    "Designs of the cables, one per cable."
    designs::D
    "Placement of each cable in the frame of the system."
    placements::Vector{Placement{T}}
    "Placed regions of the cables in the frame of the system, cable by cable."
    regions::V
    """
    Terminals, each as `(cable, terminal, connection, regions)`, cable by cable in the
    order of the terminals of each design. A row gives the cable, the row of the terminal in
    the terminals of its design, its conductor index and its regions in the system.
    """
    terminals::Vector{@NamedTuple{
        cable::Int, terminal::Int, connection::Int, regions::Vector{Int}}}
    function CableSystem{T}(id::String, line_length::T, designs::D,
            placements::Vector{Placement{T}}, regions::V,
            terminals::AbstractVector{<:NamedTuple}
    ) where {
            T <: Real,
            D <: AbstractVector{<:CableDesign},
            V <: AbstractVector{<:PlacedRegion}
    }
        return validate(new{T, D, V}(
            id, line_length, designs, placements, regions, terminals))
    end
end

Base.eltype(::CableSystem{T}) where {T} = T

# A connection specifies a conductor index, or zero, for each terminal name of its cable
# and for no other name.
function validate(
        connection::NamedTuple{<:Any, <:Tuple{Vararg{Integer}}}, design::CableDesign)
    declared = unique([row.name for row in design.terminals])
    issetequal(keys(connection), declared) || throw(ArgumentError(
        "the connection $connection of the cable $(design.id) must specify an index for " *
        "each of its terminals $declared, and for no other name"))
    all(>=(0), values(connection)) || throw(DomainError(connection,
        "the connection of the cable $(design.id) must give each terminal a conductor " *
        "index k ≥ 1, or 0 to eliminate it"))
    return connection
end

# The line has a length, and the cables do not overlap: the circles about their placements
# that contain them can touch.
function validate(system::CableSystem)
    isfinite(system.line_length) && system.line_length > zero(system.line_length) ||
        throw(DomainError(system.line_length,
            "CableSystem.line_length must be positive and finite"))
    for i in eachindex(system.designs), j in firstindex(system.designs):(i - 1)

        reach = support(system.designs[i].outer) + support(system.designs[j].outer)
        tolerance = sqrt(eps(float(reach))) * max(one(reach), reach)
        a, b = system.placements[i], system.placements[j]
        d = hypot(a.x - b.x, a.y - b.y)
        d + tolerance >= reach || throw(DomainError(d,
            "the cables $j and $i of $(system.id) overlap: their placements lie $d m " *
            "apart, closer than the $reach m that their outlines need"))
    end
    return system
end

"""
$(TYPEDSIGNATURES)

Build the system `id` of the cables of a line in one cross-section. Each cable is a design
at a placement in the frame of the system, which places the regions of the cable. The
connection of a cable maps each terminal name of its design to a conductor index. Every
terminal of that name takes the index, in each copy that declares the name, so copies that
share a name form one conductor. The system keeps the placements as declared, so cables can
touch.

# Arguments

- `designs`: designs of the cables.
- `placements`: placement of each cable in the frame of the system, its translation in
  \\[m\\] and its rotation in \\[rad\\].

# Keywords

- `connections`: one `NamedTuple` per cable that maps each terminal name of its design to a
  conductor index `k ≥ 1`, or to zero to eliminate it, such as
  `(core = 1, sheath = 2, armour = 0)`.
- `id`: name of the system.
- `line_length`: length of the line \\[m\\].

# Returns

- A [`CableSystem`](@ref).

# Errors

- Throws `DimensionMismatch` when the designs, placements and connections differ in number.
- Throws `ArgumentError` when a connection omits a terminal name of its cable, or
  specifies a name that the cable lacks.
- Throws `DomainError` when a conductor index is negative, when the line length is not
  positive and finite, or when two cables overlap.
"""
function build(::Type{CableSystem}, designs::AbstractVector{<:CableDesign},
        placements::AbstractVector{<:Placement};
        connections::AbstractVector{<:NamedTuple}, id::AbstractString, line_length::Real)
    cables = eachindex(designs, placements, connections)
    T = mapreduce(eltype, promote_type, designs;
        init = mapreduce(eltype, promote_type, placements; init = typeof(float(line_length))))
    frames = [convert(Placement{T}, placements[k]) for k in cables]
    links = [validate(connections[k], designs[k]) for k in cables]
    starts = cumsum([0; [length(designs[k].regions) for k in cables]])
    regions = [PlacedRegion(placed.source, resolve(frames[k], placed.shape))
               for k in cables for placed in designs[k].regions]
    terminals = [(cable = k, terminal = index, connection = links[k][row.name],
                     regions = row.regions .+ starts[k])
                 for k in cables for (index, row) in enumerate(designs[k].terminals)]
    return CableSystem{T}(
        String(id), convert(T, line_length), designs, frames, regions, terminals)
end
