"""
$(TYPEDEF)

Represent a closed metallic reference at the last real domain terminal.

$(TYPEDFIELDS)
"""
struct ShuntShield{T <: Real}
    "Inner radius of the shield \\[m\\]."
    radius::T
end

"""
$(TYPEDEF)

Represent the external equipotential at the outer radius of the jacket.
This reference does not own a terminal.

$(TYPEDFIELDS)
"""
struct ShuntJacket{T <: Real}
    "Outer radius of the jacket \\[m\\]."
    radius::T
end

Base.convert(::Type{Union{ShuntShield{T}, ShuntJacket{T}}}, reference::ShuntShield) where {T} =
    ShuntShield{T}(reference.radius)
Base.convert(::Type{Union{ShuntShield{T}, ShuntJacket{T}}}, reference::ShuntJacket) where {T} =
    ShuntJacket{T}(reference.radius)
nominal(reference::ShuntShield) = ShuntShield(nominal(reference.radius))
nominal(reference::ShuntJacket) = ShuntJacket(nominal(reference.radius))

"""
$(TYPEDEF)

Store the internal capacitance of one dielectric group in its owner's conductor indices:
the lossless terminal capacitance `C` \\[F/m\\] relative to the group's reference and its
potential coefficients `P` \\[m/F\\]. The shunt stage remaps the cable-local indices to the
system.

$(TYPEDFIELDS)
"""
struct InternalCapacitance{T <: Real}
    "Conductor range of the containing concentric assembly."
    assembly::UnitRange{Int}
    "Real conductors of the group, with the shield as the last terminal when present."
    terminals::UnitRange{Int}
    "Closed metallic shield or external jacket equipotential."
    reference::Union{ShuntShield{T}, ShuntJacket{T}}
    "Terminal capacitance \\[F/m\\]."
    C::Matrix{T}
    "Terminal charge-potential coefficients \\[m/F\\]."
    P::Matrix{T}
end

function Base.convert(::Type{InternalCapacitance{T}}, block::InternalCapacitance) where {T <: Real}
    InternalCapacitance{T}(block.assembly, block.terminals,
        convert(Union{ShuntShield{T}, ShuntJacket{T}}, block.reference), block.C, block.P)
end
Base.convert(::Type{InternalCapacitance{T}}, block::InternalCapacitance{T}) where {T <: Real} = block

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
    "Assembly-local conductor centers [m]."
    positions::Vector{Tuple{T, T}}
    "Equivalent conductor inner radii [m]."
    r_in::Vector{T}
    "Equivalent conductor outer radii [m]."
    r_ext::Vector{T}
    "Inner radii of conductor-owned dielectric intervals [m]."
    r_ins_in::Vector{T}
    "Outer radii of conductor-owned dielectric intervals [m]."
    r_ins_ext::Vector{T}
    "Reference conductor materials, including their resistivity calibration."
    conductor_materials::Vector{Material{T, :conductor}}
    "Conductor relative permeabilities."
    mu_r_cond::Vector{T}
    "Equivalent relative permeabilities of conductor-owned dielectric intervals."
    mu_r_ins::Vector{T}
    "Physical dielectric-layer range owned by each conductor."
    dielectric_ranges::Vector{UnitRange{Int}}
    "Physical dielectric-layer inner radii [m]."
    r_layer_in::Vector{T}
    "Physical dielectric-layer outer radii [m]."
    r_layer_ext::Vector{T}
    "Physical dielectric materials in radial order. `flatten` admits insulators and semicons
    only."
    dielectric_materials::Vector{Union{Material{T, :insulator}, Material{T, :semicon}}}
end

# Widen the representation without recalculating the blueprints or detaching
# correlated material values. Native field conversion retains unchanged indices.
function Base.convert(::Type{LocalCableData{T}}, cable::LocalCableData) where {T <: Real}
    LocalCableData{T}(map(name -> getfield(cable, name), fieldnames(typeof(cable)))...)
end
Base.convert(::Type{LocalCableData{T}}, cable::LocalCableData{T}) where {T <: Real} = cable

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
    conductor_materials = Vector{Material{T, :conductor}}(undef, conductor_count)
    mu_r_cond = Vector{T}(undef, conductor_count)
    mu_r_ins = Vector{T}(undef, conductor_count)
    dielectric_ranges = Vector{UnitRange{Int}}(undef, conductor_count)
    r_layer_in = Vector{T}(undef, layer_count)
    r_layer_ext = Vector{T}(undef, layer_count)
    dielectric_materials =
        Vector{Union{Material{T, :insulator}, Material{T, :semicon}}}(undef, layer_count)

    conductor_offset = 0
    layer_offset = 0
    assembly_offset = 0
    @inbounds for (design_index, blueprint) in pairs(blueprints)
        for (local_assembly, local_range) in pairs(blueprint.assembly_ranges)
            assembly_offset += 1
            assemblies[assembly_offset] = (
                first(local_range) + conductor_offset
            ):(last(local_range) + conductor_offset)
            assembly_designs[assembly_offset] = design_index
            for local_index in local_range
                positions[conductor_offset + local_index] =
                    blueprint.positions[local_assembly]
            end
        end
        for (local_index, row) in enumerate(blueprint)
            index = conductor_offset + local_index
            conductor = row.conductor
            terminals[index] = conductor.terminal
            r_in_values[index] = conductor.r_in
            r_ext_values[index] = conductor.r_ex
            conductor_materials[index] = conductor.material
            mu_r_cond[index] = conductor.material.mu_r

            first_layer = layer_offset + 1
            for layer in row.layers
                layer_offset += 1
                r_layer_in[layer_offset] = layer.r_in
                r_layer_ext[layer_offset] = layer.r_ex
                dielectric_materials[layer_offset] = layer.material
            end
            dielectric_ranges[index] = first_layer:layer_offset
            mu_r_ins[index] = row.mu_r_ins
            if isempty(row.layers)
                r_ins_in[index] = conductor.r_ex
                r_ins_ext[index] = conductor.r_ex
            else
                r_ins_in[index] = first(row.layers).r_in
                r_ins_ext[index] = last(row.layers).r_ex
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
        conductor_materials,
        mu_r_cond,
        mu_r_ins,
        dielectric_ranges,
        r_layer_in,
        r_layer_ext,
        dielectric_materials
    )
end

function LocalCableData(blueprint::CableBlueprint{T}) where {T <: Real}
    LocalCableData(CableBlueprint{T}[blueprint])
end
