"""
Cable geometry creation functions for the FEMTools.jl module.
These functions handle the creation of cable components.
"""

"""
$(TYPEDSIGNATURES)

Create the cable geometry for all cables in the system.

# Arguments

- `workspace`: The [`FEMWorkspace`](@ref) containing the model parameters.

# Returns

- Nothing. Updates the conductors and insulators vectors in the workspace.

# Examples

```julia
$(FUNCTIONNAME)(workspace)
```
"""
function make_cable_geometry(workspace::FEMWorkspace)

	# Get the cable system
	cable_system = workspace.problem_def.system

	# Process each cable in the system
	for (cable_idx, cable_position) in enumerate(cable_system.cables)
		@info "Processing cable $(cable_idx) at position ($(cable_position.horz), $(cable_position.vert))"

		# Get the cable design
		cable_design = cable_position.design_data

		# Get the phase assignments
		phase_assignments = cable_position.conn

		# Process each component in the cable
		for (comp_idx, component) in enumerate(cable_design.components)
			# Get the component ID
			comp_id = component.id

			# Get the phase assignment for this component
			phase = comp_idx <= length(phase_assignments) ? phase_assignments[comp_idx] : 0

			@debug "Processing component $(comp_id) (phase $(phase))"

			# Process conductor group
			if !isnothing(component.conductor_group)
				@debug "Processing conductor group for component $(comp_id)"

				# Process each layer in the conductor group
				for (layer_idx, layer) in enumerate(component.conductor_group.layers)
					@debug "Processing conductor layer $(layer_idx)"

					# Create the cable part
					_make_cablepart!(
						workspace,
						layer,
						cable_idx,
						comp_idx,
						comp_id,
						phase,
						layer_idx,
					)
				end
			end

			# Process insulator group
			if !isnothing(component.insulator_group)
				@debug "Processing insulator group for component $(comp_id)"

				# Process each layer in the insulator group
				for (layer_idx, layer) in enumerate(component.insulator_group.layers)
					@debug "Processing insulator layer $(layer_idx)"

					# Create the cable part
					_make_cablepart!(
						workspace,
						layer,
						cable_idx,
						comp_idx,
						comp_id,
						phase,
						layer_idx,
					)
				end
			end
		end
	end

	@info "Cable geometry created"
end

"""
$(TYPEDSIGNATURES)

Create a cable part entity for all tubular shapes.

# Arguments

- `workspace`: The [`FEMWorkspace`](@ref) containing the model parameters.
- `part`: The [`AbstractCablePart`](@ref) to create.
- `cable_idx`: The index of the cable.
- `comp_idx`: The index of the component.
- `comp_id`: The ID of the component.
- `phase`: The phase assignment.
- `layer_idx`: The index of the layer.

# Returns

- Nothing. Updates the conductors or insulators vector in the workspace.

# Examples

```julia
$(FUNCTIONNAME)(workspace, part, 1, 1, "core", 1, 1)
```
"""
function _make_cablepart!(workspace::FEMWorkspace, part::AbstractCablePart,
	cable_idx::Int, comp_idx::Int, comp_id::String,
	phase::Int, layer_idx::Int)

	# Get the cable definition
	cable_position = workspace.problem_def.system.cables[cable_idx]

	# Get the center coordinates
	x_center = to_nominal(cable_position.horz)
	y_center = to_nominal(cable_position.vert)

	# Determine material group directly from part type
	material_group = get_material_group(part)

	# Get or register material ID
	material_id = get_or_register_material_id(workspace, part.material_props)

	# Create physical tag with new encoding scheme
	physical_group_tag = encode_physical_group_tag(
		1,              # Surface type 1 = cable component
		cable_idx,      # Cable number
		comp_idx,       # Component number
		material_group, # Material group from part type
		material_id,     # Material ID from registry
	)

	# Create physical name
	part_type = lowercase(string(nameof(typeof(part))))
	elementary_name = create_cable_elementary_name(
		cable_idx = cable_idx,
		component_id = comp_id,
		group_type = material_group,
		part_type = part_type,
		layer_idx = layer_idx,
		phase = phase,
	)

	# Extract parameters
	radius_in = to_nominal(part.radius_in)
	radius_ext = to_nominal(part.radius_ext)

	# Calculate mesh size for this part
	if part isa AbstractConductorPart
		num_elements = workspace.formulation.elements_per_length_conductor
	elseif part isa Insulator
		num_elements = workspace.formulation.elements_per_length_insulator
	elseif part isa Semicon
		num_elements = workspace.formulation.elements_per_length_semicon
	end

	mesh_size_current =
		_calc_mesh_size(radius_in, radius_ext, part.material_props, num_elements, workspace)

	# Calculate mesh size for the next part
	num_layers =
		length(cable_position.design_data.components[comp_idx].conductor_group.layers)
	next_part =
		layer_idx < num_layers ?
		cable_position.design_data.components[comp_idx].conductor_group.layers[layer_idx+1] :
		nothing

	if !isnothing(next_part)
		next_radius_in = to_nominal(next_part.radius_in)
		next_radius_ext = to_nominal(next_part.radius_ext)
		mesh_size_next = _calc_mesh_size(
			next_radius_in,
			next_radius_ext,
			next_part.material_props,
			num_elements,
			workspace,
		)
		if next_part isa Insulator
			mesh_size = min(mesh_size_current, mesh_size_next)
		else
			mesh_size = max(mesh_size_current, mesh_size_next)
		end
	else
		mesh_size = mesh_size_current
	end

	num_points_circumference = workspace.formulation.points_per_circumference

	# Create annular shape and assign marker
	if radius_in ≈ 0
		# Solid disk
		_, _, marker, _ =
			draw_disk(x_center, y_center, radius_ext, mesh_size, num_points_circumference)
	else
		# Annular shape
		_, _, marker, _ = draw_annular(
			x_center,
			y_center,
			radius_in,
			radius_ext,
			mesh_size,
			num_points_circumference,
		)
	end

	# Create entity data
	core_data = CoreEntityData(physical_group_tag, elementary_name, mesh_size)
	entity_data = CablePartEntity(core_data, part)

	# Add to workspace in the unassigned container for subsequent processing
	workspace.unassigned_entities[marker] = entity_data

	# Add physical groups to the workspace
	register_physical_group!(workspace, physical_group_tag, part.material_props)



end

"""
$(TYPEDSIGNATURES)

Specialized method to create individual wire entities for `WireArray` parts.

# Arguments

- `workspace`: The [`FEMWorkspace`](@ref) containing the model parameters.
- `part`: The [`WireArray`](@ref)  to create.
- `cable_idx`: The index of the cable.
- `comp_idx`: The index of the component.
- `comp_id`: The ID of the component.
- `phase`: The phase assignment.
- `layer_idx`: The index of the layer.

# Returns

- Nothing. Updates the conductors vector in the workspace.

# Examples

```julia
$(FUNCTIONNAME)(workspace, part, 1, 1, "core", 1, 1)
```
"""
function _make_cablepart!(workspace::FEMWorkspace, part::WireArray,
	cable_idx::Int, comp_idx::Int, comp_id::String,
	phase::Int, layer_idx::Int)

	# Get the cable definition
	cable_position = workspace.problem_def.system.cables[cable_idx]

	# Get the center coordinates
	x_center = to_nominal(cable_position.horz)
	y_center = to_nominal(cable_position.vert)

	# Determine material group directly from part type
	material_group = get_material_group(part)

	# Get or register material ID
	material_id = get_or_register_material_id(workspace, part.material_props)

	# Create physical tag with new encoding scheme
	physical_group_tag = encode_physical_group_tag(
		1,              # Surface type 1 = cable component
		cable_idx,      # Cable number
		comp_idx,       # Component number
		material_group, # Material group from part type
		material_id,     # Material ID from registry
	)

	# -------- First handle the wires

	# Create physical name
	part_type = lowercase(string(nameof(typeof(part))))

	# Extract parameters
	radius_in = to_nominal(part.radius_in)
	radius_ext = to_nominal(part.radius_ext)

	radius_wire = to_nominal(part.radius_wire)
	num_wires = part.num_wires


	# Calculate mesh size for this part
	num_elements = workspace.formulation.elements_per_length_conductor
	mesh_size_current =
		_calc_mesh_size(radius_in, radius_ext, part.material_props, num_elements, workspace)

	# Calculate mesh size for the next part
	num_layers =
		length(cable_position.design_data.components[comp_idx].conductor_group.layers)
	next_part =
		layer_idx < num_layers ?
		cable_position.design_data.components[comp_idx].conductor_group.layers[layer_idx+1] :
		nothing

	if !isnothing(next_part)
		next_radius_in = to_nominal(next_part.radius_in)
		next_radius_ext = to_nominal(next_part.radius_ext)
		mesh_size_next = _calc_mesh_size(
			next_radius_in,
			next_radius_ext,
			next_part.material_props,
			num_elements,
			workspace,
		)
		mesh_size = max(mesh_size_current, mesh_size_next)
	else
		mesh_size = mesh_size_current
	end

	# A single wire without air gaps
	is_single_wire =
		(num_wires == 1) && (isnothing(next_part) || !(next_part isa WireArray))



	num_points_circumference = workspace.formulation.points_per_circumference

	# Calculate wire positions
	function _calc_wirearray_coords(
		num_wires::Number,
		# radius_wire::Number,
		radius_in::Number,
		radius_ext::Number;
		C = (0.0, 0.0),
	)
		wire_coords = []  # Global coordinates of all wires

		lay_radius = num_wires == 1 ? 0 : (radius_in + radius_ext) / 2

		# Calculate the angle between each wire
		angle_step = 2 * π / num_wires
		for i in 0:(num_wires-1)
			angle = i * angle_step
			x = C[1] + lay_radius * cos(angle)
			y = C[2] + lay_radius * sin(angle)
			push!(wire_coords, (x, y))  # Add wire center
		end
		return wire_coords
	end

	wire_positions =
		_calc_wirearray_coords(num_wires, radius_in, radius_ext, C = (x_center, y_center))

	# Create wires
	TOL = is_single_wire ? 0 : 5e-6 # Shrink the radius to avoid overlapping boundaries, this must be greater than Gmsh geometry tolerance
	for (wire_idx, (wx, wy)) in enumerate(wire_positions)

		_, _, marker, _ =
			draw_disk(wx, wy, radius_wire - TOL, mesh_size, num_points_circumference)

		# Create wire name
		elementary_name = create_cable_elementary_name(
			cable_idx = cable_idx,
			component_id = comp_id,
			group_type = material_group,
			part_type = part_type,
			layer_idx = layer_idx,
			phase = phase,
			wire_idx = wire_idx,
		)

		# Create entity data
		core_data = CoreEntityData(physical_group_tag, elementary_name, mesh_size)
		entity_data = CablePartEntity(core_data, part)

		# Add to workspace
		workspace.unassigned_entities[marker] = entity_data
	end
	# Add physical groups to the workspace
	register_physical_group!(workspace, physical_group_tag, part.material_props)

	# Handle WireArray outermost boundary
	mesh_size = (radius_ext - radius_in)
	if !(next_part isa WireArray) && !isnothing(next_part)
		# step_angle = 2 * pi / num_wires
		add_mesh_points(
			radius_in = radius_ext,
			radius_ext = radius_ext,
			theta_0 = 0,
			theta_1 = 2 * pi,
			mesh_size = mesh_size,
			num_points_ang = num_points_circumference,
			num_points_rad = 0,
			C = (x_center, y_center),
			theta_offset = 0, #step_angle / 2
		)
	end

	# Create air gaps for:
	# - Multiple wires (always)
	# - Single wire IF next part is a WireArray
	# Skip ONLY for single wire when next part is not a WireArray
	if !is_single_wire
		# Air gaps will be determined from the boolean fragmentation operation and do not need to be drawn. Only the markers are needed.
		markers_air_gap = get_air_gap_markers(num_wires, radius_wire, radius_in)

		# Adjust air gap markers to cable center
		for marker in markers_air_gap
			marker[1] += x_center
			marker[2] += y_center
		end

		# Determine material group - air gaps map to insulators
		material_group = 2

		# Get air material
		air_material = get_air_material(workspace)

		# Get or register material ID
		material_id = get_or_register_material_id(workspace, air_material)

		# Create physical tag with new encoding scheme
		physical_group_tag_air_gap = encode_physical_group_tag(
			1,              # Surface type 1 = cable component
			cable_idx,      # Cable number
			comp_idx,       # Component number
			material_group, # Material group from part type
			material_id,     # Material ID from registry
		)

		for marker in markers_air_gap
			# elementary names are not assigned to the air gaps because they are not drawn and appear as a result of the boolean operation
			core_data = CoreEntityData(physical_group_tag_air_gap, "", mesh_size)
			entity_data = SurfaceEntity(core_data, air_material)

			# Add to unassigned entities with type information
			workspace.unassigned_entities[marker] = entity_data
		end

		# Add physical groups to the workspace
		register_physical_group!(workspace, physical_group_tag_air_gap, air_material)
	end
end

"""
$(TYPEDSIGNATURES)

Create a cable part entity for all tubular shapes. This function is specialized for `Tubular` parts.

# Arguments

- `workspace`: The [`FEMWorkspace`](@ref) containing the model parameters.
- `part`: The [`Tubular`](@ref) to create.
- `cable_idx`: The index of the cable.
- `comp_idx`: The index of the component.
- `comp_id`: The ID of the component.
- `phase`: The phase assignment.
- `layer_idx`: The index of the layer.

# Returns

- Nothing. Updates the conductors or insulators vector in the workspace.

# Examples

```julia
$(FUNCTIONNAME)(workspace, part, 1, 1, "core", 1, 1)
```
"""
function _make_cablepart!(workspace::FEMWorkspace, part::Tubular,
	cable_idx::Int, comp_idx::Int, comp_id::String,
	phase::Int, layer_idx::Int)

	# Get the cable definition
	cable_position = workspace.problem_def.system.cables[cable_idx]

	# Get the center coordinates
	x_center = to_nominal(cable_position.horz)
	y_center = to_nominal(cable_position.vert)

	# Determine material group directly from part type
	material_group = get_material_group(part)

	# Get or register material ID
	material_id = get_or_register_material_id(workspace, part.material_props)

	# Create physical tag with new encoding scheme
	physical_group_tag = encode_physical_group_tag(
		1,              # Surface type 1 = cable component
		cable_idx,      # Cable number
		comp_idx,       # Component number
		material_group, # Material group from part type
		material_id,     # Material ID from registry
	)

	# Create physical name
	part_type = lowercase(string(nameof(typeof(part))))
	elementary_name = create_cable_elementary_name(
		cable_idx = cable_idx,
		component_id = comp_id,
		group_type = material_group,
		part_type = part_type,
		layer_idx = layer_idx,
		phase = phase,
	)

	# Extract parameters
	radius_in = to_nominal(part.radius_in)
	radius_ext = to_nominal(part.radius_ext)

	# Calculate mesh size for this part
	if part isa AbstractConductorPart
		num_elements = workspace.formulation.elements_per_length_conductor
	elseif part isa Insulator
		num_elements = workspace.formulation.elements_per_length_insulator
	elseif part isa Semicon
		num_elements = workspace.formulation.elements_per_length_semicon
	end

	mesh_size_current =
		_calc_mesh_size(radius_in, radius_ext, part.material_props, num_elements, workspace)

	# Calculate mesh size for the next part
	num_layers =
		length(cable_position.design_data.components[comp_idx].conductor_group.layers)
	next_part =
		layer_idx < num_layers ?
		cable_position.design_data.components[comp_idx].conductor_group.layers[layer_idx+1] :
		nothing

	if !isnothing(next_part)
		next_radius_in = to_nominal(next_part.radius_in)
		next_radius_ext = to_nominal(next_part.radius_ext)
		mesh_size_next = _calc_mesh_size(
			next_radius_in,
			next_radius_ext,
			next_part.material_props,
			num_elements,
			workspace,
		)
		if next_part isa Insulator
			mesh_size = min(mesh_size_current, mesh_size_next)
		else
			mesh_size = max(mesh_size_current, mesh_size_next)
		end
	else
		mesh_size = mesh_size_current
	end

	num_points_circumference = workspace.formulation.points_per_circumference

	# Create annular shape and assign marker
	if radius_in ≈ 0
		# Solid disk
		_, _, marker, _ =
			draw_disk(x_center, y_center, radius_ext, mesh_size, num_points_circumference)
	else
		# Annular shape
		_, _, marker, _ = draw_annular(
			x_center,
			y_center,
			radius_in,
			radius_ext,
			mesh_size,
			num_points_circumference,
		)
		
		# Define the inner region as an insulator (air)
        air_material = get_air_material(workspace)
        air_material_id = get_or_register_material_id(workspace, air_material)
        air_physical_group_tag = encode_physical_group_tag(1, cable_idx, comp_idx, 2, air_material_id)
        air_elementary_name = create_cable_elementary_name(
            cable_idx=cable_idx, component_id=comp_id, group_type=2,
            part_type="tubular_inner_air", layer_idx=layer_idx, phase=phase
        )
        
        # Place a marker in the inner region
        inner_marker_x = x_center
        inner_marker_y = y_center
        inner_marker = [inner_marker_x, inner_marker_y, 0.0]
        
        core_data_air = CoreEntityData(air_physical_group_tag, air_elementary_name, mesh_size)
        entity_data_air = SurfaceEntity(core_data_air, air_material)
        workspace.unassigned_entities[inner_marker] = entity_data_air
        register_physical_group!(workspace, air_physical_group_tag, air_material)

	end

	# Create entity data
	core_data = CoreEntityData(physical_group_tag, elementary_name, mesh_size)
	entity_data = CablePartEntity(core_data, part)

	# Add to workspace in the unassigned container for subsequent processing
	workspace.unassigned_entities[marker] = entity_data

	# Add physical groups to the workspace
	register_physical_group!(workspace, physical_group_tag, part.material_props)



end
