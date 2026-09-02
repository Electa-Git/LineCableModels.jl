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
    "Physical conductor cross-section [m²]."
    cross_section::T
    "Number of explicitly represented wires."
    num_wires::Int
    "Equivalent helical turns per unit length [1/m]."
    num_turns::T
    "Equivalent resistance at the material reference temperature [Ω/m]."
    resistance::T
    "Equivalent temperature coefficient [1/°C]."
    alpha::T
    "Equivalent geometric-mean radius [m]."
    gmr::T
    "Assembly-local centre in the design frame [m]."
    position::Tuple{T, T}
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
design consumed by the coaxial backend.

The blueprint is the result of computational flattening. It retains equivalent
conductor annuli and every physical dielectric layer, but contains no evaluated
constitutive relation, frequency, temperature correction, earth property, or
matrix result.

$(TYPEDFIELDS)
"""
struct CableBlueprint{T <: Real}
    "Source cable identifier."
    cable_id::String
    "Equivalent conductors in DataModel terminal order."
    conductors::Vector{BlueprintConductor{T}}
    "Physical dielectric layers in radial order."
    dielectrics::Vector{BlueprintDielectric{T}}
    "Dielectric-layer range owned by every conductor row."
    dielectric_ranges::Vector{UnitRange{Int}}
    "Contiguous conductor ranges for independent concentric assemblies."
    assembly_ranges::Vector{UnitRange{Int}}

    function CableBlueprint{T}(
            cable_id::String,
            conductors::Vector{BlueprintConductor{T}},
            dielectrics::Vector{BlueprintDielectric{T}},
            dielectric_ranges::Vector{UnitRange{Int}},
            assembly_ranges::Vector{UnitRange{Int}}
    ) where {T <: Real}
        return validate(new{T}(
            cable_id,
            conductors,
            dielectrics,
            dielectric_ranges,
            assembly_ranges
        ))
    end
end

Base.eltype(::CableBlueprint{T}) where {T} = T
Base.eltype(::Type{<:CableBlueprint{T}}) where {T} = T
Base.length(blueprint::CableBlueprint) = length(blueprint.conductors)

function _assembly_ranges(components)
    isempty(components) && throw(ArgumentError(
        "a coaxial blueprint requires at least one retained terminal",
    ))
    starts = Int[1]
    positions = [first(components).conductor.position]
    @inbounds for index in 2:length(components)
        position = components[index].conductor.position
        if !DataModel.same_radial_position(position, last(positions))
            any(reference -> DataModel.same_radial_position(position, reference), positions) &&
                throw(ArgumentError(
                    "a concentric assembly cannot reappear after another assembly",
                ))
            push!(starts, index)
            push!(positions, position)
        end
    end
    return UnitRange{Int}[start:(index == length(starts) ? length(components) :
                                 starts[index + 1] - 1)
                          for (index, start) in pairs(starts)]
end

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
    isempty(blueprint.assembly_ranges) && throw(ArgumentError(
        "CableBlueprint.assembly_ranges must contain at least one assembly",
    ))
    collect(Iterators.flatten(blueprint.assembly_ranges)) == collect(1:count) ||
        throw(DimensionMismatch(
            "CableBlueprint.assembly_ranges must partition conductor indices 1:$count in order",
        ))
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
        isfinite(conductor.cross_section) &&
        conductor.cross_section > zero(conductor.cross_section) || throw(DomainError(
            conductor.cross_section,
            "CableBlueprint.conductors[$index].cross_section must be positive and finite"
        ))
        conductor.num_wires >= 0 || throw(DomainError(
            conductor.num_wires,
            "CableBlueprint.conductors[$index].num_wires must be nonnegative"
        ))
        isfinite(conductor.num_turns) || throw(DomainError(
            conductor.num_turns,
            "CableBlueprint.conductors[$index].num_turns must be finite"
        ))
        isfinite(conductor.resistance) &&
        conductor.resistance > zero(conductor.resistance) || throw(DomainError(
            conductor.resistance,
            "CableBlueprint.conductors[$index].resistance must be positive and finite"
        ))
        isfinite(conductor.alpha) || throw(DomainError(
            conductor.alpha,
            "CableBlueprint.conductors[$index].alpha must be finite"
        ))
        isfinite(conductor.gmr) && conductor.gmr > zero(conductor.gmr) ||
            throw(DomainError(
                conductor.gmr,
                "CableBlueprint.conductors[$index].gmr must be positive and finite"
            ))
        all(isfinite, conductor.position) || throw(DomainError(
            conductor.position,
            "CableBlueprint.conductors[$index].position must be finite"
        ))
        validate(conductor.material)
        conductor.material.kind === :conductor || throw(ArgumentError(
            "CableBlueprint.conductors[$index].material.kind must be :conductor; " *
            "received $(repr(conductor.material.kind))",
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
            layer.material.kind in (:insulator, :semicon) || throw(ArgumentError(
                "CableBlueprint.dielectrics[$layer_index].material.kind must be " *
                ":insulator or :semicon; received $(repr(layer.material.kind))",
            ))
        end
    end
    return blueprint
end

"""
$(TYPEDSIGNATURES)

Flatten one completed cable design into the frequency-independent numerical
description consumed by the coaxial backend.

# Arguments

- `engine`: Coaxial backend identity.
- `design`: Completed physical cable design.
- `T`: Scalar type used by the numerical payload.

# Returns

- A validated [`CableBlueprint`](@ref) with conductor annuli, physical
  dielectric layers, and assembly partitions.
"""
function flatten(
        ::LineCableModelsCoaxial,
        design::CableDesign,
        ::Type{T}
) where {T <: Real}
    components = DataModel.radial_components(design, T)
    ranges = _assembly_ranges(components)
    assembly_by_conductor = Vector{Int}(undef, length(components))
    @inbounds for (assembly, indices) in pairs(ranges), index in indices

        assembly_by_conductor[index] = assembly
    end

    conductors = Vector{BlueprintConductor{T}}(undef, length(components))
    layer_count = sum(component -> length(component.dielectric.layers), components)
    dielectrics = Vector{BlueprintDielectric{T}}(undef, layer_count)
    dielectric_ranges = Vector{UnitRange{Int}}(undef, length(components))
    layer_index = 0
    @inbounds for (index, component) in pairs(components)
        conductor = component.conductor
        conductors[index] = BlueprintConductor{T}(
            component.name,
            assembly_by_conductor[index],
            conductor.r_in,
            conductor.r_ex,
            conductor.cross_section,
            conductor.num_wires,
            conductor.num_turns,
            conductor.resistance,
            conductor.alpha,
            conductor.gmr,
            conductor.position,
            conductor.material
        )
        first_layer = layer_index + 1
        for layer in component.dielectric.layers
            layer_index += 1
            dielectrics[layer_index] = BlueprintDielectric{T}(
                index,
                layer.r_in,
                layer.r_ex,
                layer.material
            )
        end
        dielectric_ranges[index] = first_layer:layer_index
    end
    return CableBlueprint{T}(
        design.cable_id,
        conductors,
        dielectrics,
        dielectric_ranges,
        ranges
    )
end

function flatten(
        engine::LineCableModelsCoaxial,
        design::CableDesign
)
    T = eltype(design)
    return flatten(engine, design, T)
end

"""
$(TYPEDEF)

Store the concrete array representation shared by the coaxial local primitive
impedance and potential-coefficient assemblers.

$(TYPEDFIELDS)
"""
struct LocalCableData{T <: Real}
    "Terminal names in DataModel order."
    terminals::Vector{Symbol}
    "Contiguous conductor-index ranges for concentric assemblies."
    assemblies::Vector{UnitRange{Int}}
    "Source-design index for each concentric assembly."
    assembly_designs::Vector{Int}
    "Assembly-local conductor centres [m]."
    positions::Vector{Tuple{T, T}}
    "Equivalent conductor inner radii [m]."
    r_in::Vector{T}
    "Equivalent conductor outer radii [m]."
    r_ext::Vector{T}
    "Inner radii of conductor-owned dielectric intervals [m]."
    r_ins_in::Vector{T}
    "Outer radii of conductor-owned dielectric intervals [m]."
    r_ins_ext::Vector{T}
    "Conductor resistivities at their material reference temperatures [Ω·m]."
    rho0_cond::Vector{T}
    "Conductor material reference temperatures [°C]."
    T0_cond::Vector{T}
    "Conductor temperature coefficients [1/°C]."
    alpha_cond::Vector{T}
    "Conductor relative permeabilities."
    mu_cond::Vector{T}
    "Equivalent relative permeabilities of conductor-owned dielectric intervals."
    mu_ins::Vector{T}
    "Physical dielectric-layer range owned by each conductor."
    dielectric_ranges::Vector{UnitRange{Int}}
    "Physical dielectric-layer inner radii [m]."
    r_layer_in::Vector{T}
    "Physical dielectric-layer outer radii [m]."
    r_layer_ext::Vector{T}
    "Physical dielectric materials in radial order."
    dielectric_materials::Vector{Material{T}}
    "Indices of dielectric layers classified as insulation."
    insulation_indices::Vector{Int}
    "Indices of dielectric layers classified as semiconducting material."
    semicon_indices::Vector{Int}
end

function LocalCableData(blueprints::AbstractVector{<:CableBlueprint{T}}) where {T <: Real}
    isempty(blueprints) && throw(ArgumentError(
        "local cable data require at least one blueprint",
    ))
    conductor_count = sum(length, blueprints)
    layer_count = sum(blueprint -> length(blueprint.dielectrics), blueprints)
    assembly_count = sum(blueprint -> length(blueprint.assembly_ranges), blueprints)

    terminals = Vector{Symbol}(undef, conductor_count)
    assemblies = Vector{UnitRange{Int}}(undef, assembly_count)
    assembly_designs = Vector{Int}(undef, assembly_count)
    positions = Vector{Tuple{T, T}}(undef, conductor_count)
    r_in_values = Vector{T}(undef, conductor_count)
    r_ext_values = Vector{T}(undef, conductor_count)
    r_ins_in = Vector{T}(undef, conductor_count)
    r_ins_ext = Vector{T}(undef, conductor_count)
    rho0_cond = Vector{T}(undef, conductor_count)
    T0_cond = Vector{T}(undef, conductor_count)
    alpha_cond = Vector{T}(undef, conductor_count)
    mu_cond = Vector{T}(undef, conductor_count)
    mu_ins = Vector{T}(undef, conductor_count)
    dielectric_ranges = Vector{UnitRange{Int}}(undef, conductor_count)
    r_layer_in = Vector{T}(undef, layer_count)
    r_layer_ext = Vector{T}(undef, layer_count)
    dielectric_materials = Vector{Material{T}}(undef, layer_count)
    insulation_indices = Int[]
    semicon_indices = Int[]
    sizehint!(insulation_indices, layer_count)
    sizehint!(semicon_indices, layer_count)

    conductor_offset = 0
    layer_offset = 0
    assembly_offset = 0
    @inbounds for (design_index, blueprint) in pairs(blueprints)
        for local_range in blueprint.assembly_ranges
            assembly_offset += 1
            assemblies[assembly_offset] = (
                first(local_range) + conductor_offset
            ):(last(local_range) + conductor_offset)
            assembly_designs[assembly_offset] = design_index
        end
        for local_index in eachindex(blueprint.conductors)
            index = conductor_offset + local_index
            conductor = blueprint.conductors[local_index]
            terminals[index] = conductor.terminal
            positions[index] = conductor.position
            r_in_values[index] = conductor.r_in
            r_ext_values[index] = conductor.r_ex
            rho0_cond[index] = conductor.material.rho
            T0_cond[index] = conductor.material.T0
            alpha_cond[index] = conductor.material.alpha
            mu_cond[index] = conductor.material.mu_r

            local_layers = blueprint.dielectric_ranges[local_index]
            first_layer = layer_offset + 1
            for local_layer in local_layers
                layer_offset += 1
                layer = blueprint.dielectrics[local_layer]
                r_layer_in[layer_offset] = layer.r_in
                r_layer_ext[layer_offset] = layer.r_ex
                dielectric_materials[layer_offset] = layer.material
                if layer.material.kind === :insulator
                    push!(insulation_indices, layer_offset)
                elseif layer.material.kind === :semicon
                    push!(semicon_indices, layer_offset)
                else
                    throw(ArgumentError(
                        "unsupported coaxial dielectric kind :$(layer.material.kind)",
                    ))
                end
            end
            dielectric_ranges[index] = first_layer:layer_offset
            if isempty(local_layers)
                r_ins_in[index] = conductor.r_ex
                r_ins_ext[index] = conductor.r_ex
                mu_ins[index] = one(T)
            else
                r_ins_in[index] = blueprint.dielectrics[first(local_layers)].r_in
                r_ins_ext[index] = blueprint.dielectrics[last(local_layers)].r_ex
                layers = @view blueprint.dielectrics[local_layers]
                mu_ins[index] = DataModel.equivalent_dielectric_permeability(
                    layers,
                    conductor.num_turns,
                    conductor.r_ex,
                    r_ins_ext[index]
                )
            end
        end
        conductor_offset += length(blueprint.conductors)
    end

    return LocalCableData{T}(
        terminals,
        assemblies,
        assembly_designs,
        positions,
        r_in_values,
        r_ext_values,
        r_ins_in,
        r_ins_ext,
        rho0_cond,
        T0_cond,
        alpha_cond,
        mu_cond,
        mu_ins,
        dielectric_ranges,
        r_layer_in,
        r_layer_ext,
        dielectric_materials,
        insulation_indices,
        semicon_indices
    )
end

function LocalCableData(blueprint::CableBlueprint{T}) where {T <: Real}
    LocalCableData(CableBlueprint{T}[blueprint])
end
