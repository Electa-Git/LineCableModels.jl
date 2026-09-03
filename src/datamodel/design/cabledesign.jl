"""
$(TYPEDEF)

Store one completed cable design built from an authoritative physical root.

The root is the serialized declaration. Geometry and terminal indexing are
derived together by [`build`](@ref) and cannot be supplied independently.

$(TYPEDFIELDS)
"""
struct CableDesign{
    T <: Real,
    R <: AbstractCablePart,
    G <: CableGeometry,
    N <: NamedTuple
}
    "Stable cable identifier."
    cable_id::String
    "Authoritative physical declaration."
    root::R
    "Descriptive catalogue data supplied with the physical declaration."
    nominal_data::N
    "Resolved physical geometry."
    geometry::G
    "Retained terminals in physical order."
    terminal_order::Vector{Symbol}
    "Terminal index for every resolved region; zero denotes no terminal."
    terminal_map::Vector{Int}

    function CableDesign{T, R, G, N}(
            cable_id::String,
            root::R,
            nominal_data::N,
            geometry::G,
            terminal_order::Vector{Symbol},
            terminal_map::Vector{Int}
    ) where {
            T <: Real, R <: AbstractCablePart, G <: CableGeometry, N <: NamedTuple
    }
        return validate(new{T, R, G, N}(
            cable_id,
            root,
            nominal_data,
            geometry,
            terminal_order,
            terminal_map
        ))
    end
end

Base.eltype(::CableDesign{T}) where {T} = T
Base.eltype(::Type{<:CableDesign{T}}) where {T} = T

function validate(design::CableDesign)
    isempty(design.cable_id) && throw(ArgumentError(
        "CableDesign.cable_id cannot be empty"
    ))
    isempty(design.geometry.regions) && throw(ArgumentError(
        "CableDesign.geometry.regions must contain at least one resolved region"
    ))
    isempty(design.terminal_order) && throw(ArgumentError(
        "CableDesign.terminal_order must contain at least one retained terminal"
    ))
    all(name -> !isempty(String(name)), design.terminal_order) || throw(ArgumentError(
        "CableDesign.terminal_order cannot contain an empty terminal name"
    ))
    allunique(design.terminal_order) || throw(ArgumentError(
        "CableDesign.terminal_order must contain unique terminal names; " *
        "received $(repr(design.terminal_order))"
    ))
    length(design.terminal_map) == length(design.geometry.regions) ||
        throw(DimensionMismatch(
            "CableDesign.terminal_map must contain one entry per resolved region; " *
            "received $(length(design.terminal_map)) entries for " *
            "$(length(design.geometry.regions)) regions"
        ))
    reference = first(design.geometry.regions).source.material.T0
    for (index, placed) in pairs(design.geometry.regions)
        validate(placed.source.material)
        expected = if placed.terminal === nothing
            0
        else
            something(findfirst(==(placed.terminal), design.terminal_order), 0)
        end
        design.terminal_map[index] == expected || throw(DimensionMismatch(
            "CableDesign.terminal_map[$index] must be $expected for resolved " *
            "terminal $(repr(placed.terminal)); received $(design.terminal_map[index])"
        ))
        placed.source.material.kind === :conductor && placed.terminal === nothing &&
            throw(ArgumentError(
                "CableDesign.geometry.regions[$index] is conductive but has no terminal"
            ))
        placed.source.material.kind !== :conductor && placed.terminal !== nothing &&
            throw(ArgumentError(
                "CableDesign.geometry.regions[$index] is nonconductive but carries " *
                "terminal $(repr(placed.terminal))"
            ))
        isapprox(placed.source.material.T0, reference) || throw(ArgumentError(
            "CableDesign.geometry.regions[$index].source.material.T0 must match " *
            "the design reference temperature $reference °C; received " *
            "$(placed.source.material.T0) °C"
        ))
        bounded_index = findall(
            entry -> entry.pattern isa BoundedPlacement,
            placed.placement.patterns
        )
        length(bounded_index) <= 1 || throw(ArgumentError(
            "CableDesign.geometry.regions[$index] belongs to more than one " *
            "bounded formation"
        ))
        isempty(bounded_index) && continue
        placed.source.material.kind === :conductor || throw(ArgumentError(
            "CableDesign.geometry.regions[$index] is a bounded strand but its " *
            "material is not conductive"
        ))
        placed.source.primitive isa Union{Disk, Rectangle} || throw(ArgumentError(
            "CableDesign.geometry.regions[$index] has an unsupported bounded " *
            "strand declaration $(nameof(typeof(placed.source.primitive)))"
        ))
        placed.primitive isa Union{Disk, Rectangle, Polygon, BentStrip} ||
            throw(ArgumentError(
                "CableDesign.geometry.regions[$index] has unsupported resolved " *
                "bounded geometry $(nameof(typeof(placed.primitive)))"
            ))
        declared_area = area(placed.source.primitive)
        resolved_area = area(placed.primitive)
        isapprox(
            resolved_area,
            declared_area;
            rtol = 5.0e-6,
            atol = 0
        ) || throw(ArgumentError(
            "CableDesign.geometry.regions[$index] does not preserve its " *
            "declared strand area"
        ))
        entry_index = only(bounded_index)
        entry = placed.placement.patterns[entry_index]
        entry.member == 1 || continue
        enclosing_placements = placed.placement.patterns[(entry_index + 1):end]
        formation_indices = [peer_index
                             for (peer_index, peer) in pairs(design.geometry.regions)
                             if begin
                                 peer_entries = findall(
                                     candidate -> candidate.pattern isa BoundedPlacement,
                                     peer.placement.patterns
                                 )
                                 if length(peer_entries) != 1
                                     false
                                 else
                                     peer_entry_index = only(peer_entries)
                                     peer_entry = peer.placement.patterns[peer_entry_index]
                                     isequal(
                                         peer_entry.pattern.boundary,
                                         entry.pattern.boundary
                                     ) && peer.terminal === placed.terminal &&
                                         isequal(
                                         peer.placement.patterns[
                                             (peer_entry_index + 1):end
                                         ],
                                         enclosing_placements
                                     )
                                 end
                             end]
        members = sort([
            begin
                peer = design.geometry.regions[peer_index]
                peer_entry = findfirst(
                    candidate -> candidate.pattern isa BoundedPlacement,
                    peer.placement.patterns
                )
                peer.placement.patterns[peer_entry].member
            end for peer_index in formation_indices
        ])
        members == collect(1:length(members)) ||
            throw(ArgumentError(
                "bounded-formation member identities must be contiguous from one"
            ))
        source_area = sum(
            peer_index -> area(
                design.geometry.regions[peer_index].source.primitive
            ),
            formation_indices
        )
        formation_area = sum(
            peer_index -> area(design.geometry.regions[peer_index].primitive),
            formation_indices
        )
        isapprox(
            source_area,
            formation_area;
            rtol = 5.0e-6,
            atol = 0
        ) || throw(ArgumentError(
                "a bounded formation must preserve its total declared strand area"
            ))
        boundary = entry.pattern.boundary
        boundary isa Union{Disk, SectorShape} ||
            throw(ArgumentError(
                "a bounded formation has an unsupported authoritative boundary"
            ))
        boundary_area = area(boundary)
        formation_area <= boundary_area * (1 + 5.0e-6) ||
            throw(ArgumentError(
                "bounded strand area cannot exceed its authoritative boundary"
            ))
    end
    return design
end

function Base.:(==)(left::CableDesign, right::CableDesign)
    left.cable_id == right.cable_id && left.root == right.root &&
        left.nominal_data == right.nominal_data && left.geometry == right.geometry &&
        left.terminal_order == right.terminal_order &&
        left.terminal_map == right.terminal_map
end

"""
$(TYPEDSIGNATURES)

Build a completed cable design from one v1 physical root.

The method validates physical invariants, resolves contextual geometry,
assigns retained terminals, and freezes the resulting [`CableGeometry`](@ref).
It performs no formulation calculation.

# Arguments

- `CableDesign`: Completed target type.
- `cable_id`: Stable cable identifier.
- `root`: Physical `Region`, `Stack`, `Group`, `Assembly`, or `Enclosure` root.
# Keywords

- `combine`: Gridspace composition mode. It is validated here for a common
  scalar and parametric surface; scalar construction uses one value.

# Returns

- One completed `CableDesign`.
"""
function build(
        ::Type{CableDesign},
        cable_id::AbstractString,
        parts::Tuple{Vararg{AbstractCablePart}},
        nominal_data::Union{Nothing, NamedTuple};
        combine::Symbol = :product
)
    isempty(parts) && throw(ArgumentError("a cable design requires one physical part"))
    root = length(parts) == 1 ? only(parts) : Stack(parts...)
    # 1. Public conveniences have already lowered into the physical grammar.
    # A conductive Region does not acquire a terminal from its material class.
    normalized = root

    # 2. Validate formulation-independent declarations.
    combine in (:product, :zip) || throw(ArgumentError(
        "combine must be :product or :zip"
    ))
    identifier = String(cable_id)
    isempty(identifier) && throw(ArgumentError("cable_id cannot be empty"))
    catalogue = nominal_data === nothing ? (;) : nominal_data
    catalogue isa NamedTuple || throw(ArgumentError(
        "nominal_data must be a named tuple or nothing"
    ))

    # 3. Resolve intrinsic primitives against their contextual boundaries.
    # 4. Resolve pattern placements and compaction through `placements` dispatch.
    # 5. Resolve longitudinal path radii while traversing the same physical tree.
    # 6. Collect the resolved primitives as one CableGeometry in physical order.
    geometry = resolve(EmptyBoundary(), normalized)
    isempty(geometry.regions) && throw(ArgumentError(
        "a cable design requires one region"
    ))

    # 7. Assign conductive primitives to the terminal declared by their Group.
    # 8. Establish terminal order and the region-to-terminal map.
    terminal_order = Symbol[]
    for placed in geometry.regions
        if placed.source.material.kind === :conductor && placed.terminal === nothing
            throw(ArgumentError(
                "conductive region :$(placed.source.tag) is not owned by a Group terminal"
            ))
        end
        placed.terminal === nothing || placed.terminal in terminal_order ||
            push!(terminal_order, placed.terminal)
    end
    isempty(terminal_order) && throw(ArgumentError(
        "a cable design requires at least one retained terminal"
    ))
    terminal_map = Int[placed.terminal === nothing ? 0 :
                       something(findfirst(==(placed.terminal), terminal_order), 0)
                       for placed in geometry.regions]

    reference = first(geometry.regions).source.material.T0
    all(placed -> isapprox(placed.source.material.T0, reference), geometry.regions) ||
        throw(ArgumentError("all cable materials must share one reference temperature"))

    # 9. Freeze the authoritative root and its completed geometry together.
    T = promote_type(
        (eltype(placed.primitive) for placed in geometry.regions)...,
        (eltype(placed.source.material) for placed in geometry.regions)...
    )
    return CableDesign{T, typeof(normalized), typeof(geometry), typeof(catalogue)}(
        identifier,
        normalized,
        catalogue,
        geometry,
        terminal_order,
        terminal_map
    )
end

function build(
        ::Type{CableDesign},
        cable_id::AbstractString,
        parts::AbstractCablePart...;
        nominal_data = nothing,
        combine::Symbol = :product
)
    return build(CableDesign, cable_id, parts, nominal_data; combine)
end
