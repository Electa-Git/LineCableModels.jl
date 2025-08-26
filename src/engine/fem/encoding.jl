"""
Functions for physical group tag encoding and decoding in the FEMTools.jl module.
Implements the unified SCCCOOGMMM scheme for all entity types.
"""

"""
    encode_physical_group_tag(surface_type, entity_num, component_num, material_group, material_id)

Encode entity information into a single integer ID using the unified SCCCOOGMMM scheme.

# Arguments

- `surface_type`: Surface type (1=cable, 2=physical space, 3=infinite shell) \\[dimensionless\\].
- `entity_num`: Cable number or layer number (1-999) \\[dimensionless\\].
- `component_num`: Component number (0-99, 0 for spatial regions) \\[dimensionless\\].
- `material_group`: Material group (1=conductor, 2=insulator) \\[dimensionless\\].
- `material_id`: Material identifier (0-999) \\[dimensionless\\].

# Returns

- Encoded tag as an integer \\[dimensionless\\].

# Examples

```julia
# Cable core conductor with material ID 5
tag = encode_physical_group_tag(1, 1, 1, 1, 5)  

# Air region with material ID 1
air_tag = encode_physical_group_tag(2, 1, 0, 2, 1)

# Earth layer with material ID 3
earth_tag = encode_physical_group_tag(2, 2, 0, 1, 3)
```
"""
function encode_physical_group_tag(
    surface_type::Int,
    entity_num::Int,
    component_num::Int,
    material_group::Int,
    material_id::Int
)
    # Input validation with detailed error messages
    if !(1 <= surface_type <= 9)
        error("Invalid surface type: $surface_type. Must be between 1 and 9")
    end

    if !(0 <= entity_num <= 999)
        error("Invalid entity number: $entity_num. Must be between 0 and 999")
    end

    if !(0 <= component_num <= 99)
        error("Invalid component number: $component_num. Must be between 0 and 99")
    end

    if !(1 <= material_group <= 2)
        error("""
        Invalid material group: $material_group
        Material group must be either:
        - 1: Conductor (accounts for eddy currents)
        - 2: Insulator (no eddy currents)
        """)
    end

    if !(0 <= material_id <= 99)
        error("Invalid material ID: $material_id. Must be between 0 and 99")
    end

    # SCCCOOGMM encoding
    tag = (
        surface_type * 100_000_000 +
        entity_num * 100_000 +
        component_num * 1_000 +
        material_group * 100 +
        material_id
    )

    # Validate the generated tag
    tag_str = string(tag)
    expected_length = 9  # SCCCOOGMM = 9 digits

    if length(tag_str) != expected_length
        error("Generated tag $tag has invalid length ($(length(tag_str))) for inputs: surface_type=$surface_type, entity_num=$entity_num, component_num=$component_num, material_group=$material_group, material_id=$material_id")
    end

    return tag
end

"""
    decode_physical_group_tag(tag)

Decode a physical group tag into its component parts.

# Arguments

- `tag`: Encoded tag as an integer \\[dimensionless\\].

# Returns

- Tuple of (surface_type, entity_num, component_num, material_group, material_id) \\[dimensionless\\].

# Examples

```julia
surface_type, entity_num, component_num, material_group, material_id = decode_physical_group_tag(1001010005)
println((surface_type, entity_num, component_num, material_group, material_id))  
# Output: (1, 1, 1, 1, 5)
```
"""
function decode_physical_group_tag(tag::Int)
    tag_str = string(tag)

    # Validate format
    expected_length = 9  # SCCCOOGMM = 9 digits
    if length(tag_str) != expected_length
        error("Invalid tag format: $tag. Expected a $expected_length-digit number")
    end

    # Extract parts
    surface_type = parse(Int, tag_str[1:1])
    entity_num = parse(Int, tag_str[2:4])
    component_num = parse(Int, tag_str[5:6])
    material_group = parse(Int, tag_str[7:7])
    material_id = parse(Int, tag_str[8:9])

    return (surface_type, entity_num, component_num, material_group, material_id)
end

function encode_boundary_tag(
    curve_type::Int,
    layer_idx::Int,
    sequence_num::Int=1
)
    # Input validation
    if !(1 <= curve_type <= 3)
        error("Invalid curve type: $curve_type. Must be between 1 and 3:
              1 = domain boundary
              2 = domain -> infinity
              3 = layer interface")
    end

    if !(1 <= layer_idx <= 999)
        error("Invalid layer index: $layer_idx. Must be between 1 and 999")
    end

    if !(1 <= sequence_num <= 99)
        error("Invalid sequence number: $sequence_num. Must be between 1 and 99")
    end

    # Use entity_num format consistent with the cable parts encoding
    # Format: 1CCCLSS
    # 1: Fixed prefix for boundaries/interfaces
    # CCC: Layer index (1-999)
    # L: Curve type (1-3)
    # SS: Sequence number (1-99)
    tag = 1_000_000 +
          layer_idx * 1_000 +
          curve_type * 100 +
          sequence_num

    return tag
end

function decode_boundary_tag(tag::Int)
    tag_str = string(tag)

    # Validate format (should start with 1)
    if length(tag_str) != 7 || tag_str[1] != '1'
        error("Invalid boundary tag format: $tag. Expected a 7-digit number starting with 1")
    end

    # Extract parts
    layer_idx = parse(Int, tag_str[2:4])
    curve_type = parse(Int, tag_str[5])
    sequence_num = parse(Int, tag_str[6:7])

    return (curve_type, layer_idx, sequence_num)
end

"""
    get_material_group(part)

Get the material group (conductor or insulator) for a cable part based on its type.

# Arguments

- `part`: An AbstractCablePart instance.

# Returns

- Material group (1=conductor, 2=insulator) \\[dimensionless\\].

# Examples

```julia
group = get_material_group(wire_array)  # Returns 1 (conductor)
```
"""
function get_material_group(part::AbstractCablePart)
    if part isa AbstractConductorPart
        return 1  # Conductor
    elseif part isa AbstractInsulatorPart
        return 2  # Insulator
    else
        error("Unknown part type: $(typeof(part))")
    end
end

"""
    get_material_group(earth_model, layer_idx)

Get the material group (conductor or insulator) for an earth layer.

# Arguments

- `earth_model`: The EarthModel containing layer information.
- `layer_idx`: The layer index to check.

# Returns

- Material group (1=conductor, 2=insulator) \\[dimensionless\\].

# Examples

```julia
group = get_material_group(earth_model, 1)  # Layer 1 is air -> Returns 2 (insulator)
group = get_material_group(earth_model, 2)  # Layer 2 is earth -> Returns 1 (conductor)
```
"""
function get_material_group(earth_model::EarthModel, layer_idx::Int)
    # Layer 1 is always air (insulator)
    if layer_idx == 1
        return 2  # Insulator
    else
        # All other layers are earth (conductor)
        return 1  # Conductor
    end
end

"""
    get_or_register_material_id(workspace, material)

Find or create a unique ID for a material within the current workspace.

# Arguments

- `workspace`: The FEMWorkspace containing the material registry.
- `material`: The Material object to register.

# Returns

- A unique material ID (1-99) \\[dimensionless\\].

# Examples

```julia
material_id = get_or_register_material_id(workspace, copper_material)
```
"""
function get_or_register_material_id(workspace::FEMWorkspace, material::Material)
    # Create material_registry if it doesn't exist
    if !isdefined(workspace, :material_registry)
        workspace.material_registry = Dict{String,Int}()
    end

    # Get material name using existing function that checks library first
    material_name = get_material_name(material, workspace.formulation.materials)

    # Find or create the ID
    if !haskey(workspace.material_registry, material_name)
        # New material - assign next available ID
        material_id = length(workspace.material_registry) + 1
        if material_id > 99
            error("Material registry full: Maximum of 99 unique materials supported")
        end
        workspace.material_registry[material_name] = material_id
    else
        material_id = workspace.material_registry[material_name]
    end

    return material_id
end

function register_physical_group!(workspace::FEMWorkspace, physical_group_tag::Int, material::Material)

    # Create physical_groups if it doesn't exist
    if !isdefined(workspace, :physical_groups)
        workspace.physical_groups = Dict{Int,Material}()
    end

    # Find or create the ID
    if !haskey(workspace.physical_groups, physical_group_tag)
        # New material - assign next available ID
        workspace.physical_groups[physical_group_tag] = material
    end

end
"""
$(TYPEDSIGNATURES)

Generate a readable elementary name for a cable component.
Format: cable_X_<component>_<group>_layer_<Y>_<part>[_wire_N][_phase_M]

# Arguments

- `cable_idx`: Cable index \\[dimensionless\\].
- `component_id`: Component ID (e.g., "core", "sheath") \\[dimensionless\\].
- `group_type`: Group type (1=conductor, 2=insulator, 3=empty) \\[dimensionless\\].
- `part_type`: Part type (e.g., "wire", "strip", "tubular") \\[dimensionless\\].
- `layer_idx`: Layer index \\[dimensionless\\].
- `wire_idx`: Optional wire index \\[dimensionless\\].
- `phase`: Optional phase index \\[dimensionless\\].

# Returns

- Human-readable physical name as a string.

# Examples

```julia
name = $(FUNCTIONNAME)(
    cable_idx=1,
    component_id="core",
    group_type=1,
    part_type="wire",
    layer_idx=2,
    wire_idx=3,
    phase=1
)
println(name)  # Output: "cable_1_core_con_layer_2_wire_wire_3_phase_1"
```
"""
function create_cable_elementary_name(;
    cable_idx::Int,
    component_id::String,
    group_type::Int,  # 1=conductor, 2=insulator, 3=air gap
    part_type::String,
    layer_idx::Union{Int,Nothing}=nothing,
    wire_idx::Union{Int,Nothing}=nothing,
    phase::Union{Int,Nothing}=nothing
)
    # Convert group_type to string
    group_str = if group_type == 1
        "con"
    elseif group_type == 2
        "ins"
    else
        error("Invalid group_type: $group_type")
    end

    # Base name without optional parts
    name = "cable_$(cable_idx)_$(component_id)_$(group_str)"

    # Add layer index if provided
    if !isnothing(layer_idx)
        name *= "_layer_$(layer_idx)"
    end

    name *= "_$(part_type)"

    # Add wire index if provided
    if !isnothing(wire_idx)
        name *= "_wire_$(wire_idx)"
    end

    # Add phase if provided
    if !isnothing(phase) && phase > 0
        name *= "_phase_$(phase)"
    elseif !isnothing(phase) && phase == 0
        name *= "_ground"
    end

    return name
end

function create_physical_group_name(workspace::FEMWorkspace, tag::Int)
    # Determine tag type by length
    tag_str = string(tag)

    if length(tag_str) == 9
        # This is a physical group tag (SCCCOOGMM format)
        return _create_surface_physical_name(workspace, tag)
    elseif length(tag_str) == 7 && tag_str[1] == '1'
        # This is a boundary tag (1CCCLSS format)
        return _create_boundary_physical_name(workspace, tag)
    else
        # Unknown format - return generic name
        return "group_$(tag)"
    end
end

function _create_surface_physical_name(workspace::FEMWorkspace, tag::Int)
    # Decode the tag
    surface_type, entity_num, component_num, material_group, material_id = decode_physical_group_tag(tag)

    # Get material name if available
    material_name = "unknown"
    for (name, id) in workspace.material_registry
        if id == material_id
            material_name = name
            break
        end
    end

    # Create base string based on surface type
    base_str = if surface_type == 1
        # Cable component
        # Try to get component name
        component_name = "unknown"
        if 1 <= entity_num <= length(workspace.problem_def.system.cables)
            cable = workspace.problem_def.system.cables[entity_num]

            # Validate component_num is within range
            if 1 <= component_num <= length(cable.design_data.components)
                component = cable.design_data.components[component_num]
                component_name = component.id
            end
        end

        group_str = material_group == 1 ? "con" : "ins"
        "cable_$(entity_num)_$(component_name)_$(group_str)"
    elseif surface_type == 2
        # Physical domain
        layer_str = entity_num == 1 ? "air" : "earth"
        group_str = material_group == 1 ? "con" : "ins"
        "layer_$(entity_num)_$(layer_str)_$(group_str)"
    elseif surface_type == 3
        # Infinite shell
        layer_str = entity_num == 1 ? "air" : "earth"
        group_str = material_group == 1 ? "con" : "ins"
        "infshell_$(entity_num)_$(layer_str)_$(group_str)"
    else
        "surf_$(surface_type)"
    end

    # Add material information
    return "$(base_str)_$(material_name)"
end

function _create_boundary_physical_name(workspace::FEMWorkspace, tag::Int)
    # Decode the boundary tag
    curve_type, layer_idx, sequence_num = decode_boundary_tag(tag)

    # Create boundary name based on curve type
    base_str = if curve_type == 1
        # Domain boundary
        layer_str = layer_idx == 1 ? "air" : "earth"
        "boundary_domain_$(layer_str)"
    elseif curve_type == 2
        # Domain to infinity
        layer_str = layer_idx == 1 ? "air" : "earth"
        "boundary_infinity_$(layer_str)"
    elseif curve_type == 3
        # Layer interface
        if layer_idx == 1
            "interface_air_earth"
        else
            "interface_earth_layers_$(layer_idx)_$(layer_idx+1)"
        end
    else
        "boundary_unknown"
    end

    # Add sequence number if more than one of the same type
    if sequence_num > 1
        base_str *= "_$(sequence_num)"
    end

    return base_str
end
