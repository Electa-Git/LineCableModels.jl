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
        N
}
    "Stable cable identifier."
    cable_id::String
    "Optional catalogue data."
    nominal_data::N
    "Authoritative physical declaration."
    root::R
    "Resolved physical geometry."
    geometry::G
    "Retained terminals in physical order."
    terminal_order::Vector{Symbol}
    "Terminal index for every resolved region; zero denotes no terminal."
    terminal_map::Vector{Int}

    function CableDesign{T, R, G, N}(
            cable_id::String,
            nominal_data::N,
            root::R,
            geometry::G,
            terminal_order::Vector{Symbol},
            terminal_map::Vector{Int}
    ) where {T <: Real, R <: AbstractCablePart, G <: CableGeometry, N}
        return new{T, R, G, N}(
            cable_id,
            nominal_data,
            root,
            geometry,
            terminal_order,
            terminal_map
        )
    end
end

Base.eltype(::CableDesign{T}) where {T} = T
Base.eltype(::Type{<:CableDesign{T}}) where {T} = T

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
- `nominal_data`: Optional catalogue data.

# Keywords

- `combine`: Gridspace composition mode. It is validated here for a common
  scalar and parametric surface; scalar construction uses one value.

# Returns

- One completed `CableDesign`.
"""
function build(
        ::Type{CableDesign},
        cable_id::AbstractString,
        root::AbstractCablePart,
        nominal_data::Union{Nothing, NominalData};
        combine::Symbol = :product
)
    # 1. Normalize the public single-conductor shorthand into the v1 grammar.
    normalized = if root isa Region && root.material.kind === :conductor
        Group(root.tag, Pose2(0, 0, 0), root, nothing, nothing, nothing)
    else
        root
    end

    # 2. Validate formulation-independent declarations.
    combine in (:product, :zip) || throw(ArgumentError(
        "combine must be :product or :zip"
    ))
    identifier = String(cable_id)
    isempty(identifier) && throw(ArgumentError("cable_id cannot be empty"))

    # 3. Resolve intrinsic primitives against their contextual boundaries.
    # 4. Resolve pattern placements and compaction through `placements` dispatch.
    # 5. Resolve longitudinal path radii while traversing the same physical tree.
    # 6. Collect the resolved shapes as one CableGeometry in physical order.
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
    terminal_map = Int[
        placed.terminal === nothing ? 0 :
        something(findfirst(==(placed.terminal), terminal_order), 0)
        for placed in geometry.regions
    ]

    reference = first(geometry.regions).source.material.T0
    all(placed -> isapprox(placed.source.material.T0, reference), geometry.regions) ||
        throw(ArgumentError("all cable materials must share one reference temperature"))

    # 9. Freeze the authoritative root and its completed geometry together.
    T = promote_type(
        (eltype(placed.shape) for placed in geometry.regions)...,
        (eltype(placed.source.material) for placed in geometry.regions)...
    )
    nominal = nominal_data === nothing ? nothing : convert(NominalData{T}, nominal_data)
    return CableDesign{T, typeof(normalized), typeof(geometry), typeof(nominal)}(
        identifier,
        nominal,
        normalized,
        geometry,
        terminal_order,
        terminal_map
    )
end

function Base.show(io::IO, ::MIME"text/plain", design::CableDesign)
    print(
        io,
        "CableDesign \"$(design.cable_id)\": [regions=$(length(design.geometry.regions)), " *
        "terminals=($(join(design.terminal_order, ", "))), " *
        "outer_radius=$(round(outer_radius(design), sigdigits=5))]"
    )
end
