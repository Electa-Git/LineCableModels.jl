"""
$(TYPEDEF)

Store one frequency-independent equivalent coaxial conductor row.

$(TYPEDFIELDS)
"""
struct BlueprintConductor{T <: Real}
    "Retained terminal name."
    terminal::Symbol
    "Concentric assembly containing the terminal."
    assembly::Int
    "Equivalent inner radius [m]."
    r_in::T
    "Equivalent outer radius [m]."
    r_ex::T
    "Artificial homogeneous conductor material."
    material::Material{T, :conductor}
end

"""
$(TYPEDEF)

Store one physical dielectric layer owned by a coaxial conductor interval.

$(TYPEDFIELDS)
"""
struct BlueprintDielectric{T <: Real}
    "Index of the conductor immediately inside this radial interval."
    conductor::Int
    "Layer inner radius [m]."
    r_in::T
    "Layer outer radius [m]."
    r_ex::T
    "Unmodified physical material."
    material::Material{T}
end

"""
$(TYPEDEF)

Store the frequency-independent, unreduced numerical description of one cable
design.

The blueprint is the result of [`flatten`](@ref). It retains equivalent conductor
annuli and every physical dielectric layer, and it depends on the design and the
scalar type only. Frequency-dependent constitutive evaluation, conductor
temperature correction, earth return and matrix reductions remain calculation work.

`blueprint[i]` returns the row of conductor `i`, `(conductor, layers, mu_r_ins)`, with
`layers` a view of the dielectric layers of its interval. Iteration gives the rows in
terminal order. `eltype` stays the scalar type `T`.

$(TYPEDFIELDS)
"""
struct CableBlueprint{T <: Real}
    "Source cable identifier."
    cable_id::String
    "Equivalent conductors in terminal order."
    conductors::Vector{BlueprintConductor{T}}
    "Physical dielectric layers in radial order."
    dielectrics::Vector{BlueprintDielectric{T}}
    "Dielectric-layer range owned by every conductor row."
    dielectric_ranges::Vector{UnitRange{Int}}
    "Equivalent relative permeability of the dielectric interval of every conductor row."
    mu_r_ins::Vector{T}
    "Contiguous conductor ranges for independent concentric assemblies."
    assembly_ranges::Vector{UnitRange{Int}}
    "Center of each concentric assembly [m], in the frame of the design."
    positions::Vector{Tuple{T, T}}

    function CableBlueprint{T}(
            cable_id::String,
            conductors::Vector{BlueprintConductor{T}},
            dielectrics::Vector{BlueprintDielectric{T}},
            dielectric_ranges::Vector{UnitRange{Int}},
            mu_r_ins::Vector{T},
            assembly_ranges::Vector{UnitRange{Int}},
            positions::Vector{Tuple{T, T}}
    ) where {T <: Real}
        return new{T}(
            cable_id,
            conductors,
            dielectrics,
            dielectric_ranges,
            mu_r_ins,
            assembly_ranges,
            positions
        )
    end
end

Base.eltype(::CableBlueprint{T}) where {T} = T
Base.eltype(::Type{<:CableBlueprint{T}}) where {T} = T
Base.length(blueprint::CableBlueprint) = length(blueprint.conductors)
Base.firstindex(blueprint::CableBlueprint) = firstindex(blueprint.conductors)
Base.lastindex(blueprint::CableBlueprint) = lastindex(blueprint.conductors)
Base.eachindex(blueprint::CableBlueprint) = eachindex(blueprint.conductors)
function Base.getindex(blueprint::CableBlueprint, index::Integer)
    return (conductor = blueprint.conductors[index],
        layers = view(blueprint.dielectrics, blueprint.dielectric_ranges[index]),
        mu_r_ins = blueprint.mu_r_ins[index])
end
function Base.iterate(blueprint::CableBlueprint, index::Int = firstindex(blueprint))
    index > lastindex(blueprint) && return nothing
    return blueprint[index], index + 1
end
# `eltype` is the scalar type, as for a design, and iteration does not declare an element
# type. `IteratorEltype(Any)` gives the unknown trait through a public name.
Base.IteratorEltype(::Type{<:CableBlueprint}) = Base.IteratorEltype(Any)
