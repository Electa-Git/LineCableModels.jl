"""
Domain creation functions for the FEMTools.jl module.
These functions handle the creation of domain boundaries and earth interfaces.
"""

"""
$(TYPEDSIGNATURES)

Create the domain boundaries (inner solid disk and outer annular region) for the simulation.

# Arguments

- `workspace`: The [`FEMWorkspace`](@ref) containing the model parameters.

# Returns

- Nothing. Updates the boundaries vector in the workspace.

# Examples

```julia
$(FUNCTIONNAME)(workspace)
```
"""
function make_space_geometry(workspace::FEMWorkspace)
    @info "Creating domain boundaries..."

    # Extract parameters
    formulation = workspace.core.formulation
    domain_radius = formulation.domain_radius
    domain_radius_inf = formulation.domain_radius_inf  # External radius for boundary transform
    mesh_size_default = formulation.mesh_size_default
    mesh_size_domain = formulation.mesh_size_max
    mesh_size_inf = 1.25 * formulation.mesh_size_max

    # Center coordinates
    x_center = 0.0
    y_center = 0.0
    eps = 1e-6

    # Create inner domain disk
    num_points_circumference = formulation.points_per_circumference
    @debug "Creating inner domain disk with radius $(domain_radius) m"
    _, _, air_region_marker, domain_boundary_markers = draw_disk(
        x_center,
        y_center,
        domain_radius,
        mesh_size_domain,
        num_points_circumference,
    )

    # Create outer domain annular region
    @debug "Creating outer domain annular region with radius $(domain_radius_inf) m"
    _, _, air_infshell_marker, domain_infty_markers = draw_annular(
        x_center,
        y_center,
        domain_radius,
        domain_radius_inf,
        mesh_size_inf,
        num_points_circumference,
    )

    # Get earth model from workspace
    earth_props = workspace.core.earth_props
    air_layer_idx = 1 # air layer is 1 by default
    num_earth_layers = length(earth_props.layers) # Number of earth layers

    # Air layer (Layer 1)
    air_material = get_earth_model_material(workspace, air_layer_idx)
    air_material_id = get_or_register_material_id(workspace, air_material)
    air_material_group = get_material_group(earth_props, air_layer_idx) # Will return 2 (insulator)

    # Physical domain air tag 
    air_region_tag = encode_physical_group_tag(
        2,                  # Surface type 2 = physical domain
        air_layer_idx,      # Layer 1 = air
        0,                  # Component 0 (not a cable component)
        air_material_group, # Material group 2 (insulator)
        air_material_id,    # Material ID
    )
    air_region_name = create_physical_group_name(workspace, air_region_tag)

    # Infinite shell air tag
    air_infshell_tag = encode_physical_group_tag(
        3,                  # Surface type 3 = infinite shell
        air_layer_idx,      # Layer 1 = air
        0,                  # Component 0 (not a cable component)
        air_material_group, # Material group 2 (insulator)
        air_material_id,    # Material ID
    )
    air_infshell_name = create_physical_group_name(workspace, air_infshell_tag)

    # Create group tags for boundary curves - above ground (air) - inner domain
    air_boundary_tag = encode_boundary_tag(1, air_layer_idx, 1)
    air_boundary_name = create_physical_group_name(workspace, air_boundary_tag)
    air_boundary_marker = [0.0, domain_radius, 0.0]

    # Above ground (air) - domain -> infinity
    air_infty_tag = encode_boundary_tag(2, air_layer_idx, 1)
    air_infty_name = create_physical_group_name(workspace, air_infty_tag)
    air_infty_marker = [0.0, domain_radius_inf, 0.0]

    # Create boundary curves
    air_boundary_entity = CurveEntity(
        CoreEntityData(air_boundary_tag, air_boundary_name, mesh_size_domain),
        air_material,
    )

    air_infty_entity = CurveEntity(
        CoreEntityData(air_infty_tag, air_infty_name, mesh_size_inf),
        air_material,
    )
    # Add curves to the workspace
    workspace.core.unassigned_entities[air_boundary_marker] = air_boundary_entity
    workspace.core.unassigned_entities[air_infty_marker] = air_infty_entity


    # Create domain surfaces
    air_region_entity = SurfaceEntity(
        CoreEntityData(air_region_tag, air_region_name, mesh_size_default),
        air_material,
    )

    air_infshell_entity = SurfaceEntity(
        CoreEntityData(air_infshell_tag, air_infshell_name, mesh_size_default),
        air_material,
    )
    # Add surfaces to the workspace
    workspace.core.unassigned_entities[air_region_marker] = air_region_entity
    workspace.core.unassigned_entities[air_infshell_marker] = air_infshell_entity

    # Add physical groups to the workspace
	register_physical_group!(workspace, air_region_tag, air_material)
	register_physical_group!(workspace, air_infshell_tag, air_material)
   
    # Below ground (earth) - inner domain
	earth_material = get_earth_model_material(workspace, num_earth_layers)
    earth_boundary_tag = encode_boundary_tag(1, num_earth_layers, 1)
    earth_boundary_name = create_physical_group_name(workspace, earth_boundary_tag)
    earth_boundary_entity = CurveEntity(
        CoreEntityData(earth_boundary_tag, earth_boundary_name, mesh_size_domain),
        earth_material,
    )
    
    # Below ground (earth) - domain -> infinity
    earth_infty_tag = encode_boundary_tag(2, num_earth_layers, 1)
    earth_infty_name = create_physical_group_name(workspace, earth_infty_tag)
    earth_infty_entity = CurveEntity(
        CoreEntityData(earth_infty_tag, earth_infty_name, mesh_size_inf),
        earth_material,
    )
    
    @debug "Domain boundary markers:"
    for point_marker in domain_boundary_markers
        target_entity = point_marker[2] > 0 ? air_boundary_entity : earth_boundary_entity
        workspace.core.unassigned_entities[point_marker] = target_entity
        @debug "  Point $point_marker: ($(point_marker[1]), $(point_marker[2]), $(point_marker[3]))"
    end
    
    @debug "Domain -> infinity markers:"
    for point_marker in domain_infty_markers
        target_entity = point_marker[2] > 0 ? air_infty_entity : earth_infty_entity
        workspace.core.unassigned_entities[point_marker] = target_entity
        @debug "  Point $point_marker: ($(point_marker[1]), $(point_marker[2]), $(point_marker[3]))"
    end
    
    # Physical groups for Dirichlet boundary
    register_physical_group!(workspace, earth_infty_tag, earth_material)
    register_physical_group!(workspace, air_infty_tag, air_material)
    
    # Interface Earth Air
    num_elements = formulation.elements_per_length_interfaces
	
    # Earth layer
    current_x_start = -domain_radius
    current_y_start = 0.0

    for layer_idx in 2:num_earth_layers
        # Material properties
        earth_material = get_earth_model_material(workspace, layer_idx)
        earth_material_id = get_or_register_material_id(workspace, earth_material)
        earth_material_group = get_material_group(earth_props, layer_idx) # Will return 1 (conductor)

        ## Physical domain earth tag
        earth_region_tag = encode_physical_group_tag(
            2,                    # Surface type 2 = physical domain
            layer_idx,      # Layer 2 = first earth layer
            0,                    # Component 0 (not a cable component)
            earth_material_group, # Material group 1 (conductor)
            earth_material_id,    # Material ID
        )
        earth_region_name = create_physical_group_name(workspace, earth_region_tag)

        ## Infinite shell earth tag
        earth_infshell_tag = encode_physical_group_tag(
            3,                    # Surface type 3 = infinite shell
            layer_idx,      # Layer 2 = first earth layer
            0,                    # Component 0 (not a cable component)
            earth_material_group, # Material group 1 (conductor)
            earth_material_id,    # Material ID
        )
        earth_infshell_name = create_physical_group_name(workspace, earth_infshell_tag)

        # Add physical groups to the workspace
        register_physical_group!(workspace, earth_region_tag, earth_material)
        register_physical_group!(workspace, earth_infshell_tag, earth_material)

        # Determine layer thickness
        layer_thickness = earth_props.layers[layer_idx].t == Inf ? domain_radius : earth_props.layers[layer_idx].t

        if earth_props.vertical_layers

            next_x_start = clamp(current_x_start + layer_thickness, 0, domain_radius)
            @assert next_x_start <= domain_radius "Layer $layer_idx extends beyond domain radius $next_x_start"

            earth_boundary_marker = (layer_idx == 2) ?
                [next_x_start-eps, -sqrt(domain_radius^2 - (next_x_start-eps)^2), 0.0] :
                [next_x_start+eps, -sqrt(domain_radius^2 - (next_x_start-eps)^2), 0.0]

            earth_infty_marker = (layer_idx == 2) ?
                [next_x_start-eps, -sqrt(domain_radius_inf^2 - (next_x_start-eps)^2), 0.0] :
                [next_x_start+eps, -sqrt(domain_radius_inf^2 - (next_x_start-eps)^2), 0.0]

            earth_region_marker = [next_x_start-eps, -eps, 0.0]

            workspace.core.unassigned_entities[earth_boundary_marker] = earth_boundary_entity
            workspace.core.unassigned_entities[earth_infty_marker] = earth_infty_entity
            earth_infshell_marker = [[0.99*(next_x_start-eps), -0.99*sqrt(domain_radius_inf^2 - (next_x_start-eps)^2), 0.0]]

            interface_idx = layer_idx
            earth_interface_tag = encode_boundary_tag(3, interface_idx, 1)
            earth_interface_name = create_physical_group_name(workspace, earth_interface_tag)
            if layer_idx < num_earth_layers
                earth_interface_mesh_size = _calc_mesh_size(0, domain_radius, earth_material, num_elements, workspace)
                y_pos = sqrt(domain_radius^2 - (next_x_start)^2)
                y_pos_inf = sqrt(domain_radius_inf^2 - (next_x_start)^2)

                # Interface in earth layer 
                _, _, earth_interface_markers = draw_line(
                    next_x_start,
                    0.0,
                    next_x_start,
                    -y_pos,
                    mesh_size_domain,
                    round(Int, y_pos),
                )
                # Interface in infinite shell
                _, _, earth_inter_markers_2 = draw_line(
                    next_x_start,
                    -y_pos,
                    next_x_start,
                    -y_pos_inf,
                    mesh_size_domain,
                    round(Int, y_pos_inf-y_pos),
                )

                append!(earth_interface_markers, earth_inter_markers_2)
                earth_interface_entity = CurveEntity(
                    CoreEntityData(
                        earth_interface_tag,
                        earth_interface_name,
                        mesh_size_domain,
                    ),
                    get_earth_model_material(workspace, layer_idx),  # Earth material
                )

                _ = gmsh.model.occ.add_point(next_x_start, 0.0, 0.0, earth_interface_mesh_size)
                _ = gmsh.model.occ.add_point(next_x_start, -y_pos, 0.0, earth_interface_mesh_size)
                _ = gmsh.model.occ.add_point(next_x_start, -y_pos_inf, 0.0, earth_interface_mesh_size)

                @debug "Domain interface vertical layers markers:"
                for point_marker in earth_interface_markers
                    workspace.core.unassigned_entities[point_marker] = earth_interface_entity
                    @debug "  Point $point_marker: ($(point_marker[1]), $(point_marker[2]), $(point_marker[3]))"
                end
                current_x_start = next_x_start
            end

        else # Horizontal and Uniform layers

            next_y_start = clamp(current_y_start - layer_thickness, -domain_radius, 0.0)
            @assert next_y_start >= -domain_radius "Layer $layer_idx extends beyond domain radius"

            x_pos_inf_current = sqrt(domain_radius_inf^2 - (current_y_start - eps)^2)
            
            earth_boundary_marker = [[sqrt(domain_radius^2 - (next_y_start/2)^2), next_y_start/2, 0.0]]
            push!(earth_boundary_marker, [-sqrt(domain_radius^2 - (next_y_start/2)^2), (next_y_start/2), 0.0]) 

            earth_infty_marker = [[sqrt(domain_radius_inf^2 - (next_y_start/2)^2), next_y_start/2, 0.0]]
            push!(earth_infty_marker, [-sqrt(domain_radius_inf^2 - (next_y_start/2)^2), (next_y_start/2), 0.0])


            for marker in earth_boundary_marker
                workspace.core.unassigned_entities[marker] = earth_boundary_entity
            end
            for marker in earth_infty_marker
                workspace.core.unassigned_entities[marker] = earth_infty_entity
            end

            earth_region_marker = [0.0, current_y_start - eps, 0.0]
            
            # For the outer domain (the infinite shell)
            earth_infshell_marker = [[-x_pos_inf_current+eps, current_y_start - eps, 0.0]]
            push!(earth_infshell_marker, [x_pos_inf_current-eps, current_y_start - eps, 0.0])
            
            if layer_idx < num_earth_layers
                x_pos = sqrt(domain_radius^2 - next_y_start^2)
                x_pos_inf = sqrt(domain_radius_inf^2 - next_y_start^2)

                interface_idx = layer_idx # The interface belongs to the layer below it
                earth_interface_tag = encode_boundary_tag(3, interface_idx, 1)
                earth_interface_name = create_physical_group_name(workspace, earth_interface_tag)
                earth_interface_mesh_size = _calc_mesh_size(0, domain_radius, earth_material, num_elements, workspace)

                # Interface between -domain_radius to domain_radius
                _, _, earth_interface_markers = draw_line(
                    -x_pos,
                    next_y_start,
                    x_pos,
                    next_y_start,
                    earth_interface_mesh_size,
                    round(Int, 2 * 2*x_pos / 10),
                )
                # Interface between -domain_radius_inf to -domain_radius
                _, _, earth_inter_markers_2 = draw_line(
                    - x_pos_inf,
                    next_y_start,
                    - x_pos,
                    next_y_start,
                    earth_interface_mesh_size,
                    round(Int, 2 * (x_pos_inf - x_pos) / 10),
                )
                # Interface between +domain_radius_inf to +domain_radius
                _, _, earth_inter_markers_3 = draw_line(
                    + x_pos_inf,
                    next_y_start,
                    + x_pos,
                    next_y_start,
                    earth_interface_mesh_size,
                    round(Int, 2 * (x_pos_inf - x_pos) / 10),
                )
                append!(earth_interface_markers, earth_inter_markers_2)
                append!(earth_interface_markers, earth_inter_markers_3)
                marker_tag = gmsh.model.occ.add_point(-x_pos_inf, next_y_start, 0.0, earth_interface_mesh_size)
                marker_tag = gmsh.model.occ.add_point(-x_pos, next_y_start, 0.0, earth_interface_mesh_size)
                marker_tag = gmsh.model.occ.add_point( x_pos, next_y_start, 0.0, earth_interface_mesh_size)
                marker_tag = gmsh.model.occ.add_point( x_pos_inf, next_y_start, 0.0, earth_interface_mesh_size)

                # Creates the interface curve entity
                interface_entity = CurveEntity(
                    CoreEntityData(
                        earth_interface_tag,
                        earth_interface_name,
                        earth_interface_mesh_size,
                    ),
                    get_earth_model_material(workspace, layer_idx), # Material of the next layer
                )

                # Associates the line points with the interface entity
                @debug "Domain interface horizontal layer markers at y = $(next_y_start):"
                for point_marker in earth_interface_markers
                    workspace.core.unassigned_entities[point_marker] = interface_entity
                    @debug "  Point $point_marker: ($(point_marker[1]), $(point_marker[2]), $(point_marker[3]))"
                end
            end

            current_y_start = next_y_start
        end
        marker_tag = gmsh.model.occ.add_point(
            earth_region_marker[1],
            earth_region_marker[2],
            earth_region_marker[3],
            mesh_size_domain,
        )
        gmsh.model.set_entity_name(
            0,
            marker_tag,
            "marker_$(round(mesh_size_domain, sigdigits=6))",
        )
        # Create markers for all the infinite shell
        for point in earth_infshell_marker
            marker_tag = gmsh.model.occ.add_point(
                point[1],
                point[2],
                point[3],
                mesh_size_inf,
            )

            gmsh.model.set_entity_name(0, marker_tag, "marker_$(round(mesh_size_inf, sigdigits=6))")
        end
        earth_region_entity = SurfaceEntity(
            CoreEntityData(earth_region_tag, earth_region_name, mesh_size_default),
            earth_material,
        )

        earth_infshell_entity = SurfaceEntity(
            CoreEntityData(earth_infshell_tag, earth_infshell_name, mesh_size_default),
            earth_material,
        )

        workspace.core.unassigned_entities[earth_region_marker] = earth_region_entity
        for marker in earth_infshell_marker
            workspace.core.unassigned_entities[marker] = earth_infshell_entity
        end


    end
    earth_material = get_earth_model_material(workspace, num_earth_layers)
    earth_interface_mesh_size = _calc_mesh_size(0, domain_radius, earth_material, num_elements, workspace)
	_, _, earth_interface_markers = draw_line(
		-domain_radius_inf,
		0.0,
		domain_radius_inf,
		0.0,
		earth_interface_mesh_size,
		round(Int, domain_radius),
	)

	interface_idx = 1   # Earth interface index
	earth_interface_tag = encode_boundary_tag(3, interface_idx, 1)
	earth_interface_name = create_physical_group_name(workspace, earth_interface_tag)

	# Create domain entity
	earth_interface_entity = CurveEntity(
		CoreEntityData(
			earth_interface_tag,
			earth_interface_name,
			earth_interface_mesh_size,
		),
		get_earth_model_material(workspace, num_earth_layers),  # Earth material
	)


	# Create mesh transitions if specified
	if !isempty(workspace.core.formulation.mesh_transitions)
		@info "Creating $(length(workspace.core.formulation.mesh_transitions)) mesh transition regions"

		for (idx, transition) in enumerate(workspace.core.formulation.mesh_transitions)
			cx, cy = transition.center
			transition_radii =
				collect(LinRange(transition.r_min, transition.r_max, transition.n_regions))
			# Use provided layer or auto-detect
			layer_idx = if !isnothing(transition.earth_layer)
				transition.earth_layer
			else
				# Fallback auto-detection (should rarely happen due to constructor)
                get_cable_layer(cx, cy, transition_radii[end], earth_props, domain_radius)
			end

			# Validate layer index exists in earth model
			if layer_idx > num_earth_layers
				Base.error(
					"Earth layer $layer_idx does not exist in earth model (max: $(num_earth_layers))",
				)
			end

			# Get material for this earth layer
			transition_material = get_earth_model_material(workspace, layer_idx)
			material_id = get_or_register_material_id(workspace, transition_material)
			material_group = get_material_group(earth_props, layer_idx)

			# Create physical tag for this transition
			transition_tag = encode_physical_group_tag(
				2,                # Surface type 2 = physical domain
				layer_idx,        # Earth layer index
				0,                # Component 0 (not a cable component)
				material_group,   # Material group (1=conductor for earth, 2=insulator for air)
				material_id,       # Material ID
			)

			layer_name = layer_idx == 1 ? "air" : "earth_$(layer_idx-1)"
			transition_name = "mesh_transition_$(idx)_$(layer_name)"

			# Calculate radii and mesh sizes
			mesh_size_min = transition.mesh_factor_min * earth_interface_mesh_size
			mesh_size_max = transition.mesh_factor_max * earth_interface_mesh_size

			transition_mesh =
				collect(LinRange(mesh_size_min, mesh_size_max, transition.n_regions))
			@debug "Transition $(idx): radii=$(transition_radii), mesh sizes=$(transition_mesh)"

			# Draw the transition regions
			_, _, transition_markers = draw_transition_region(
				cx, cy,
				transition_radii,
				transition_mesh,
				num_points_circumference,
			)

			# Register each transition region
			for k in 1:transition.n_regions
				transition_region = SurfaceEntity(
					CoreEntityData(
						transition_tag,
						"$(transition_name)_region_$(k)",
						transition_mesh[k],
					),
					transition_material,
				)
				workspace.core.unassigned_entities[transition_markers[k]] = transition_region

				@debug "Created transition region $k at ($(cx), $(cy)) with radius $(transition_radii[k]) m in layer $layer_idx"
			end

			# Register physical group
			register_physical_group!(workspace, transition_tag, transition_material)
		end

		@info "Mesh transition regions created"
	else
		@debug "No mesh transitions specified"
	end

    # Add interface to the workspace
    @debug "Domain -> infinity markers:"
	for point_marker in earth_interface_markers
		workspace.core.unassigned_entities[point_marker] = earth_interface_entity
		@debug "  Point $point_marker: ($(point_marker[1]), $(point_marker[2]), $(point_marker[3]))"
	end

    @info "Earth interfaces created"
end

"""
    get_cable_layer(x::Number, y::Number, earth_model::EarthModel) -> Int

Determines the index of the soil layer from the `EarthModel` where the cable 
centroid located at `(x, y)` is situated.

The function assumes Layer 1 is always air for any `y >= 0`. It then checks if the
model uses vertical or horizontal layering to determine the correct layer index
for `y < 0`.

# Arguments
- `x::Number`: The horizontal coordinate of the cable's centroid.
- `y::Number`: The vertical coordinate of the cable's centroid (depth, negative).
- `outermost_radius::Float64`: The outermost radius of the cable.
- `earth_model::EarthModel`: The earth model structure, which must contain a 
  `layers` vector and a `vertical_layers` boolean.
- `domain_radius::Float64`: The radius of the simulation domain.

# Returns
- `Int`: The index of the layer (1 for air, 2 for the first earth layer, etc.).
"""
function get_cable_layer(x::Number, y::Number, outermost_radius::Float64, earth_model::EarthModel, domain_radius::Float64)
    
    number_of_points = round(Int, 2 * pi * outermost_radius / 0.01)
    layer_idx_in_circle = Int64[]
    for θ in LinRange(0, 2*pi, number_of_points)
        x_p = x + outermost_radius * cos(θ)
        y_p = y + outermost_radius * sin(θ)
    
        if (x_p^2 + y_p^2) > domain_radius^2
            # This point is irrelevant to the simulation. Skip to the next one.
            Base.error("Mesh Transition outside domain radius!") 
        end
        # Air layer (above ground)
        if y_p >= 0.0
            push!(layer_idx_in_circle, 1)
            continue
        end

        num_layers = length(earth_model.layers)

        # Homogeneous earth (only 1 earth layer)
        if num_layers == 2
            push!(layer_idx_in_circle, 2)
            continue
        end

        # Multi-layer earth
        if earth_model.vertical_layers
            # Vertical layers: extend horizontally
            current_x = 0.0  # Starting boundary for layer 2
            layer_found = false 
            for layer_idx in 2:num_layers
                layer_thickness = earth_model.layers[layer_idx].t

                # Determine next boundary
                next_x = isinf(layer_thickness) ? 0.0 : current_x + layer_thickness

                if layer_idx == 2
                    if x_p >= -domain_radius && x_p <= next_x
                        push!(layer_idx_in_circle, layer_idx)
                        layer_found = true
                        break
                    end
                else
                    if x_p > current_x && x_p <= next_x
                        push!(layer_idx_in_circle, layer_idx)
                        layer_found = true
                        break
                    end
                end

                # Update boundary for next iteration
                current_x = next_x
            end

            if !layer_found
                # Fallback: return last layer
                push!(layer_idx_in_circle, num_layers)
            end

        else
            # Horizontal layers: extend vertically downward
            current_y = 0.0  # Starting depth
            layer_found = false 
            for layer_idx in 2:num_layers
                layer_thickness = earth_model.layers[layer_idx].t
                
                # Determine next depth boundary
                next_y = isinf(layer_thickness) ? -domain_radius : current_y - layer_thickness

                if y_p <= current_y && y_p >= next_y
                    push!(layer_idx_in_circle, layer_idx)
                    layer_found = true
                    break
                end
                
                # Update boundary for next iteration
                current_y = next_y
            end
            if !layer_found
                # Fallback: return last layer
                push!(layer_idx_in_circle, num_layers)
            end

        end
    end
    if !isempty(layer_idx_in_circle) && all(==(layer_idx_in_circle[1]), layer_idx_in_circle)
        # The condition passed, do nothing or proceed.
        @debug "All mesh elements are self-contained in a single layer."
        return layer_idx_in_circle[1]
    else
        # The condition failed, throw an error.
        Base.error("Not all mesh elements are self-contained in a single layer!")
    end
end