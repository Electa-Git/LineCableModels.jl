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
    formulation = workspace.formulation
    domain_radius = formulation.domain_radius
    domain_radius_inf = formulation.domain_radius_inf  # External radius for boundary transform
    mesh_size_default = formulation.mesh_size_default
    mesh_size_domain = formulation.mesh_size_max
    mesh_size_inf = 1.25 * formulation.mesh_size_max

    # Center coordinates
    x_center = 0.0
    y_center = 0.0

    # Create inner domain disk
    num_points_circumference = formulation.points_per_circumference
    @debug "Creating inner domain disk with radius $(domain_radius) m"
    _, _, air_region_marker, domain_boundary_markers = draw_disk(x_center, y_center, domain_radius, mesh_size_domain, num_points_circumference)

    # Create outer domain annular region
    @debug "Creating outer domain annular region with radius $(domain_radius_inf) m"
    _, _, air_infshell_marker, domain_infty_markers = draw_annular(x_center, y_center, domain_radius, domain_radius_inf, mesh_size_inf, num_points_circumference)

    # Get earth model from workspace
    earth_props = workspace.problem_def.earth_props
    air_layer_idx = 1 # air layer is 1 by default
    num_earth_layers = length(earth_props.layers) # Number of earth layers
    earth_layer_idx = num_earth_layers

    # Air layer (Layer 1)
    air_material = get_earth_model_material(workspace, air_layer_idx)
    air_material_id = get_or_register_material_id(workspace, air_material)
    air_material_group = get_material_group(earth_props, air_layer_idx) # Will return 2 (insulator)

    # Physical domain air tag 
    air_region_tag = encode_physical_group_tag(
        2,                # Surface type 2 = physical domain
        air_layer_idx,    # Layer 1 = air
        0,                # Component 0 (not a cable component)
        air_material_group, # Material group 2 (insulator)
        air_material_id   # Material ID
    )
    air_region_name = create_physical_group_name(workspace, air_region_tag)

    # Infinite shell air tag
    air_infshell_tag = encode_physical_group_tag(
        3,                # Surface type 3 = infinite shell
        air_layer_idx,    # Layer 1 = air
        0,                # Component 0 (not a cable component)
        air_material_group, # Material group 2 (insulator)
        air_material_id   # Material ID
    )
    air_infshell_name = create_physical_group_name(workspace, air_infshell_tag)


    # Earth layer (Layer 2+)
    earth_material = get_earth_model_material(workspace, earth_layer_idx)
    earth_material_id = get_or_register_material_id(workspace, earth_material)
    earth_material_group = get_material_group(earth_props, earth_layer_idx) # Will return 1 (conductor)

    # Physical domain earth tag
    earth_region_tag = encode_physical_group_tag(
        2,                  # Surface type 2 = physical domain
        earth_layer_idx,    # Layer 2 = first earth layer
        0,                  # Component 0 (not a cable component)
        earth_material_group, # Material group 1 (conductor)
        earth_material_id   # Material ID
    )
    earth_region_name = create_physical_group_name(workspace, earth_region_tag)

    # Infinite shell earth tag
    earth_infshell_tag = encode_physical_group_tag(
        3,                  # Surface type 3 = infinite shell
        earth_layer_idx,    # Layer 2 = first earth layer
        0,                  # Component 0 (not a cable component)
        earth_material_group, # Material group 1 (conductor)
        earth_material_id   # Material ID
    )
    earth_infshell_name = create_physical_group_name(workspace, earth_infshell_tag)


    # Create group tags for boundary curves - above ground (air) - inner domain
    air_boundary_tag = encode_boundary_tag(1, air_layer_idx, 1)
    air_boundary_name = create_physical_group_name(workspace, air_boundary_tag)
    air_boundary_marker = [0.0, domain_radius, 0.0]

    # Below ground (earth) - inner domain
    earth_boundary_tag = encode_boundary_tag(1, earth_layer_idx, 1)
    earth_boundary_name = create_physical_group_name(workspace, earth_boundary_tag)
    earth_boundary_marker = [0.0, -domain_radius, 0.0]

    # Above ground (air) - domain -> infinity
    air_infty_tag = encode_boundary_tag(2, air_layer_idx, 1)
    air_infty_name = create_physical_group_name(workspace, air_infty_tag)
    air_infty_marker = [0.0, domain_radius_inf, 0.0]

    # Below ground (earth) - domain -> infinity
    earth_infty_tag = encode_boundary_tag(2, earth_layer_idx, 1)
    earth_infty_name = create_physical_group_name(workspace, earth_infty_tag)
    earth_infty_marker = [0.0, -domain_radius_inf, 0.0]

    # Create markers for the domain surfaces
    earth_region_marker = [0.0, -domain_radius * 0.99, 0.0]
    marker_tag = gmsh.model.occ.add_point(earth_region_marker[1], earth_region_marker[2], earth_region_marker[3], mesh_size_domain)
    gmsh.model.set_entity_name(0, marker_tag, "marker_$(round(mesh_size_domain, sigdigits=6))")

    earth_infshell_marker = [0.0, -(domain_radius + 0.99 * (domain_radius_inf - domain_radius)), 0.0]
    marker_tag = gmsh.model.occ.add_point(earth_infshell_marker[1], earth_infshell_marker[2], earth_infshell_marker[3], mesh_size_inf)
    gmsh.model.set_entity_name(0, marker_tag, "marker_$(round(mesh_size_inf, sigdigits=6))")

    # Create boundary curves
    air_boundary_entity = CurveEntity(
        CoreEntityData(air_boundary_tag, air_boundary_name, mesh_size_domain),
        air_material
    )

    earth_boundary_entity = CurveEntity(
        CoreEntityData(earth_boundary_tag, earth_boundary_name, mesh_size_domain),
        earth_material
    )

    air_infty_entity = CurveEntity(
        CoreEntityData(air_infty_tag, air_infty_name, mesh_size_inf),
        air_material
    )

    earth_infty_entity = CurveEntity(
        CoreEntityData(earth_infty_tag, earth_infty_name, mesh_size_inf),
        earth_material
    )

    # Add curves to the workspace
    workspace.unassigned_entities[air_boundary_marker] = air_boundary_entity
    workspace.unassigned_entities[air_infty_marker] = air_infty_entity
    workspace.unassigned_entities[earth_boundary_marker] = earth_boundary_entity
    workspace.unassigned_entities[earth_infty_marker] = earth_infty_entity

    @debug "Domain boundary markers:"
    for point_marker in domain_boundary_markers
        target_entity = point_marker[2] > 0 ? air_boundary_entity : earth_boundary_entity
        workspace.unassigned_entities[point_marker] = target_entity
        @debug "  Point $point_marker: ($(point_marker[1]), $(point_marker[2]), $(point_marker[3]))"
    end

    @debug "Domain -> infinity markers:"
    for point_marker in domain_infty_markers
        target_entity = point_marker[2] > 0 ? air_infty_entity : earth_infty_entity
        workspace.unassigned_entities[point_marker] = target_entity
        @debug "  Point $point_marker: ($(point_marker[1]), $(point_marker[2]), $(point_marker[3]))"
    end

    # Add physical groups to the workspace
    register_physical_group!(workspace, air_region_tag, air_material)
    register_physical_group!(workspace, earth_region_tag, earth_material)
    register_physical_group!(workspace, air_infshell_tag, air_material)
    register_physical_group!(workspace, earth_infshell_tag, earth_material)

    # Physical groups for Dirichlet boundary
    register_physical_group!(workspace, air_infty_tag, air_material)
    register_physical_group!(workspace, earth_infty_tag, earth_material)

    # Create domain surfaces
    air_region_entity = SurfaceEntity(
        CoreEntityData(air_region_tag, air_region_name, mesh_size_default),
        air_material
    )

    air_infshell_entity = SurfaceEntity(
        CoreEntityData(air_infshell_tag, air_infshell_name, mesh_size_default),
        air_material
    )

    # Earth regions will be created after boolean fragmentation
    earth_region_entity = SurfaceEntity(
        CoreEntityData(earth_region_tag, earth_region_name, mesh_size_default),
        earth_material
    )

    earth_infshell_entity = SurfaceEntity(
        CoreEntityData(earth_infshell_tag, earth_infshell_name, mesh_size_default),
        earth_material
    )

    # Add surfaces to the workspace
    workspace.unassigned_entities[air_region_marker] = air_region_entity
    workspace.unassigned_entities[air_infshell_marker] = air_infshell_entity
    workspace.unassigned_entities[earth_region_marker] = earth_region_entity
    workspace.unassigned_entities[earth_infshell_marker] = earth_infshell_entity

    @info "Domain boundaries created"

    # Create earth interface line (y=0)
    @debug "Creating earth interface line at y=0"

    # Create line from -domain_radius to +domain_radius at y=0
    num_elements = formulation.elements_per_length_interfaces
    earth_interface_mesh_size = _calc_mesh_size(0, domain_radius, earth_material, num_elements, workspace)

    _, _, earth_interface_markers = draw_line(-domain_radius_inf, 0.0, domain_radius_inf, 0.0, earth_interface_mesh_size, round(Int, domain_radius))

    # Create physical tag for the earth interface
    interface_idx = 1  # Earth interface index
    earth_interface_tag = encode_boundary_tag(3, interface_idx, 1)
    earth_interface_name = create_physical_group_name(workspace, earth_interface_tag)

    # Create domain entity
    earth_interface_entity = CurveEntity(
        CoreEntityData(earth_interface_tag, earth_interface_name, earth_interface_mesh_size),
        get_earth_model_material(workspace, earth_layer_idx)  # Earth material
    )

    # Create mesh transitions if specified
    if !isempty(workspace.formulation.mesh_transitions)
        @info "Creating $(length(workspace.formulation.mesh_transitions)) mesh transition regions"

        for (idx, transition) in enumerate(workspace.formulation.mesh_transitions)
            cx, cy = transition.center

            # Use provided layer or auto-detect
            layer_idx = if !isnothing(transition.earth_layer)
                transition.earth_layer
            else
                # Fallback auto-detection (should rarely happen due to constructor)
                cy >= 0 ? 1 : 2
            end

            # Validate layer index exists in earth model
            if layer_idx > num_earth_layers
                error("Earth layer $layer_idx does not exist in earth model (max: $(num_earth_layers))")
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
                material_id       # Material ID
            )

            layer_name = layer_idx == 1 ? "air" : "earth_$(layer_idx-1)"
            transition_name = "mesh_transition_$(idx)_$(layer_name)"

            # Calculate radii and mesh sizes
            mesh_size_min = transition.mesh_factor_min * earth_interface_mesh_size
            mesh_size_max = transition.mesh_factor_max * earth_interface_mesh_size

            transition_radii = collect(LinRange(transition.r_min, transition.r_max, transition.n_regions))
            transition_mesh = collect(LinRange(mesh_size_min, mesh_size_max, transition.n_regions))
            @debug "Transition $(idx): radii=$(transition_radii), mesh sizes=$(transition_mesh)"

            # Draw the transition regions
            _, _, transition_markers = draw_transition_region(
                cx, cy,
                transition_radii,
                transition_mesh,
                num_points_circumference
            )

            # Register each transition region
            for k in 1:transition.n_regions
                transition_region = SurfaceEntity(
                    CoreEntityData(transition_tag, "$(transition_name)_region_$(k)", transition_mesh[k]),
                    transition_material
                )
                workspace.unassigned_entities[transition_markers[k]] = transition_region

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
        workspace.unassigned_entities[point_marker] = earth_interface_entity
        @debug "  Point $point_marker: ($(point_marker[1]), $(point_marker[2]), $(point_marker[3]))"
    end

    @info "Earth interfaces created"

end