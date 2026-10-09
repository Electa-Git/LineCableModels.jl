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
    material::Material{T}
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
        return validate(new{T}(
            cable_id,
            conductors,
            dielectrics,
            dielectric_ranges,
            mu_r_ins,
            assembly_ranges,
            positions
        ))
    end
end

Base.eltype(::CableBlueprint{T}) where {T} = T
Base.eltype(::Type{<:CableBlueprint{T}}) where {T} = T
Base.length(blueprint::CableBlueprint) = length(blueprint.conductors)

function validate(blueprint::CableBlueprint)
    isempty(blueprint.cable_id) && throw(ArgumentError(
        "CableBlueprint.cable_id cannot be empty"
    ))
    count = length(blueprint.conductors)
    count > 0 || throw(ArgumentError(
        "CableBlueprint.conductors must contain at least one conductor"
    ))
    length(blueprint.dielectric_ranges) == count || throw(DimensionMismatch(
        "CableBlueprint.dielectric_ranges must contain one range per conductor; " *
        "received $(length(blueprint.dielectric_ranges)) ranges for $count conductors",
    ))
    length(blueprint.mu_r_ins) == count || throw(DimensionMismatch(
        "CableBlueprint.mu_r_ins must contain one value per conductor; " *
        "received $(length(blueprint.mu_r_ins)) values for $count conductors",
    ))
    isempty(blueprint.assembly_ranges) && throw(ArgumentError(
        "CableBlueprint.assembly_ranges must contain at least one assembly",
    ))
    collect(Iterators.flatten(blueprint.assembly_ranges)) == collect(1:count) ||
        throw(DimensionMismatch(
            "CableBlueprint.assembly_ranges must partition conductor indices 1:$count in order",
        ))
    length(blueprint.positions) == length(blueprint.assembly_ranges) ||
        throw(DimensionMismatch(
            "CableBlueprint.positions must contain one center per assembly; received " *
            "$(length(blueprint.positions)) for $(length(blueprint.assembly_ranges)) assemblies",
        ))
    @inbounds for (assembly, position) in pairs(blueprint.positions)
        all(isfinite, position) || throw(DomainError(
            position,
            "CableBlueprint.positions[$assembly] must be finite"
        ))
    end
    layer_count = length(blueprint.dielectrics)
    collected_layers = collect(Iterators.flatten(blueprint.dielectric_ranges))
    collected_layers == collect(1:layer_count) || throw(DimensionMismatch(
        "CableBlueprint.dielectric_ranges must partition dielectric indices " *
        "1:$layer_count in order",
    ))
    @inbounds for (index, conductor) in pairs(blueprint.conductors)
        isempty(String(conductor.terminal)) && throw(ArgumentError(
            "CableBlueprint.conductors[$index].terminal cannot be empty"
        ))
        conductor.assembly in eachindex(blueprint.assembly_ranges) || throw(DomainError(
            conductor.assembly,
            "CableBlueprint.conductors[$index].assembly must index assembly_ranges"
        ))
        index in blueprint.assembly_ranges[conductor.assembly] || throw(DimensionMismatch(
            "CableBlueprint.conductors[$index].assembly does not own conductor $index",
        ))
        isfinite(conductor.r_in) && conductor.r_in >= zero(conductor.r_in) ||
            throw(DomainError(
                conductor.r_in,
                "CableBlueprint.conductors[$index].r_in must be nonnegative and finite"
            ))
        isfinite(conductor.r_ex) && conductor.r_ex > conductor.r_in ||
            throw(DomainError(
                conductor.r_ex,
                "CableBlueprint.conductors[$index].r_ex must be finite and greater than r_in"
            ))
        validate(conductor.material)
        isfinite(conductor.material.rho) || throw(DomainError(
            conductor.material.rho,
            "CableBlueprint.conductors[$index].material.rho must be finite"
        ))
        mu_r = blueprint.mu_r_ins[index]
        isfinite(mu_r) && mu_r > zero(mu_r) || throw(DomainError(
            mu_r,
            "CableBlueprint.mu_r_ins[$index] must be positive and finite"
        ))
        for layer_index in blueprint.dielectric_ranges[index]
            layer = blueprint.dielectrics[layer_index]
            layer.conductor == index || throw(DimensionMismatch(
                "CableBlueprint.dielectrics[$layer_index].conductor must be $index; " *
                "received $(layer.conductor)",
            ))
            isfinite(layer.r_in) && layer.r_in >= zero(layer.r_in) ||
                throw(DomainError(
                    layer.r_in,
                    "CableBlueprint.dielectrics[$layer_index].r_in must be nonnegative and finite"
                ))
            isfinite(layer.r_ex) && layer.r_ex > layer.r_in ||
                throw(DomainError(
                    layer.r_ex,
                    "CableBlueprint.dielectrics[$layer_index].r_ex must be finite and greater than r_in"
                ))
            validate(layer.material)
        end
    end
    return blueprint
end
