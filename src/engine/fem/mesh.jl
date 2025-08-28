"""
Mesh generation functions for the FEMTools.jl module.
These functions handle the configuration and generation of the mesh.
"""

"""
$(TYPEDSIGNATURES)

Calculate the skin depth for a conductive material.

# Arguments

- `rho`: Electrical resistivity \\[Ω·m\\].
- `mu_r`: Relative permeability \\[dimensionless\\].
- `freq`: Frequency \\[Hz\\].

# Returns

- Skin depth \\[m\\].

# Examples

```julia
depth = $(FUNCTIONNAME)(1.7241e-8, 1.0, 50.0)
```

# Notes

```math
\\delta = \\sqrt{\\frac{\\rho}{\\pi \\cdot f \\cdot \\mu_0 \\cdot \\mu_r}}
```

where \\(\\mu_0 = 4\\pi \\times 10^{-7}\\) H/m is the vacuum permeability.
"""
function calc_skin_depth(rho::Number, mu_r::Number, freq::Number)
    # Convert to nominal values in case of Measurement types
    rho = to_nominal(rho)
    mu_r = to_nominal(mu_r)

    # Constants
    mu_0 = 4e-7 * π  # Vacuum permeability

    # Calculate skin depth
    # δ = sqrt(ρ / (π * f * μ_0 * μ_r))
    return sqrt(rho / (π * freq * mu_0 * mu_r))
end

function _calc_mesh_size(part::AbstractCablePart, workspace::FEMWorkspace)

    # Extract geometric properties
    radius_in = to_nominal(part.radius_in)
    radius_ext = to_nominal(part.radius_ext)
    thickness = radius_ext - radius_in

    # Extract formulation parameters
    formulation = workspace.formulation

    # Calculate mesh size based on part type and properties
    scale_length = thickness
    if part isa WireArray
        # For wire arrays, consider the wire radius
        scale_length = to_nominal(part.radius_wire) * 2
        num_elements = formulation.elements_per_length_conductor
    elseif part isa AbstractConductorPart
        num_elements = formulation.elements_per_length_conductor
    elseif part isa Insulator
        num_elements = formulation.elements_per_length_insulator
    elseif part isa Semicon
        num_elements = formulation.elements_per_length_semicon
    end

    # Apply bounds from configuration
    mesh_size = scale_length / num_elements
    mesh_size = max(mesh_size, formulation.mesh_size_min)
    mesh_size = min(mesh_size, formulation.mesh_size_max)

    return mesh_size
end

function _calc_mesh_size(radius_in::Number, radius_ext::Number, material::Material, num_elements::Int, workspace::FEMWorkspace)
    # Extract geometric properties
    thickness = radius_ext - radius_in

    # Extract problem_def parameters
    formulation = workspace.formulation
    mesh_size = thickness / num_elements

    # Apply bounds from configuration
    mesh_size = max(mesh_size, formulation.mesh_size_min)
    mesh_size = min(mesh_size, formulation.mesh_size_max)

    return mesh_size
end

"""
$(TYPEDSIGNATURES)

Configure mesh sizes for all entities in the model.

# Arguments

- `workspace`: The [`FEMWorkspace`](@ref) containing the entities.

# Returns

- Nothing. Updates the mesh size map in the workspace.

# Examples

```julia
$(FUNCTIONNAME)(workspace)
```
"""
function config_mesh_options(workspace::FEMWorkspace)


    gmsh.option.set_number("General.InitialModule", 2)

    # Set mesh algorithm
    gmsh.option.set_number("Mesh.Algorithm", workspace.formulation.mesh_algorithm)
    gmsh.option.set_number("Mesh.AlgorithmSwitchOnFailure", 1)
    # Set mesh optimization parameters
    gmsh.option.set_number("Mesh.Optimize", 0)
    gmsh.option.set_number("Mesh.OptimizeNetgen", 0)

    # Set mesh globals
    gmsh.option.set_number("Mesh.SaveAll", 1)  # Mesh all regions
    gmsh.option.set_number("Mesh.MaxRetries", workspace.formulation.mesh_max_retries)
    gmsh.option.set_number("Mesh.MeshSizeMin", workspace.formulation.mesh_size_min)
    gmsh.option.set_number("Mesh.MeshSizeMax", workspace.formulation.mesh_size_max)
    gmsh.option.set_number("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.set_number("Mesh.MeshSizeFromParametricPoints", 0)

    gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 1)
    gmsh.option.set_number("Mesh.MeshSizeFromCurvature", workspace.formulation.points_per_circumference)


    @debug "Mesh algorithm: $(workspace.formulation.mesh_algorithm)"
    @debug "Mesh size range: [$(workspace.formulation.mesh_size_min), $(workspace.formulation.mesh_size_max)]"
end

"""
$(TYPEDSIGNATURES)

Generate the mesh.

# Arguments

- `workspace`: The [`FEMWorkspace`](@ref) containing the model.

# Returns

- Nothing. Generates the mesh in the Gmsh model.

# Examples

```julia
$(FUNCTIONNAME)(workspace)
```
"""
function generate_mesh(workspace::FEMWorkspace)
    # Generate 2D mesh
    gmsh.model.mesh.generate(2)

    # Get mesh statistics
    nodes = gmsh.model.mesh.get_nodes()
    elements = gmsh.model.mesh.get_elements()

    num_nodes = length(nodes[1])
    num_elements = sum(length.(elements[2]))

    @info "Mesh generation completed"
    @info "Created mesh with $(num_nodes) nodes and $(num_elements) elements"
end

"""
$(TYPEDSIGNATURES)

Initialize a Gmsh model with appropriate settings.

# Arguments

- `case_id`: Identifier for the model.
- `problem_def`: The [`FEMFormulation`](@ref) containing mesh parameters.
- `solver`: The [`FEMSolver`](@ref) containing visualization parameters.

# Returns

- Nothing. Initializes the Gmsh model.

# Examples

```julia
$(FUNCTIONNAME)("test_case", problem_def, solver)
```
"""
function initialize_gmsh(workspace::FEMWorkspace)
    # Create a new model
    system_id = workspace.problem_def.system.system_id
    gmsh.model.add(system_id)

    # Module launched on startup (0: automatic, 1: geometry, 2: mesh, 3: solver, 4: post-processing)
    gmsh.option.set_number("General.InitialModule", 0)
    gmsh.option.set_string("General.DefaultFileName", system_id * ".geo")

    # Define verbosity level
    gmsh_verbosity = map_verbosity_to_gmsh(workspace.opts.verbosity)
    gmsh.option.set_number("General.Verbosity", gmsh_verbosity)

    # Set OCC model healing options
    gmsh.option.set_number("Geometry.AutoCoherence", 1)
    gmsh.option.set_number("Geometry.OCCFixDegenerated", 1)
    gmsh.option.set_number("Geometry.OCCFixSmallEdges", 1)
    gmsh.option.set_number("Geometry.OCCFixSmallFaces", 1)
    gmsh.option.set_number("Geometry.OCCSewFaces", 1)
    gmsh.option.set_number("Geometry.OCCMakeSolids", 1)

    # Log settings based on verbosity
    @info "Initialized Gmsh model: $system_id"

end

function _do_make_mesh!(workspace::FEMWorkspace)

    # Initialize Gmsh model and set parameters
    initialize_gmsh(workspace)

    # Create geometry
    @info "Creating domain boundaries..."
    make_space_geometry(workspace)

    @info "Creating cable geometry..."
    make_cable_geometry(workspace)

    # Synchronize the model
    gmsh.model.occ.synchronize()

    # Boolean operations
    @info "Performing boolean operations..."
    process_fragments(workspace)

    # Entity identification and entity assignment
    @info "Identifying entities after fragmentation..."
    identify_by_marker(workspace)

    # Physical group assignment
    @info "Assigning physical groups..."
    assign_physical_groups(workspace)

    # Mesh sizing
    @info "Setting up mesh sizing..."
    config_mesh_options(workspace)

    # Mesh generation
    @info "Generating mesh..."
    generate_mesh(workspace)

    # Save mesh
    @info "Saving mesh to file: $(display_path(workspace.paths[:mesh_file]))"
    gmsh.write(workspace.paths[:mesh_file])

    # Save geometry
    @info "Saving geometry to file: $(display_path(workspace.paths[:geo_file]))"
    gmsh.write(workspace.paths[:geo_file])
end

function mesh_exists(workspace::FEMWorkspace)
    mesh_file = workspace.paths[:mesh_file]

    # Force remesh overrides everything
    if workspace.opts.force_remesh
        @debug "Force remesh requested"
        return false
    end

    # If workspace is empty (no entities), force remesh regardless of file existence
    if isempty(workspace.conductors) && isempty(workspace.insulators) && isempty(workspace.space_regions) && isempty(workspace.boundaries) && isempty(workspace.physical_groups) && isempty(workspace.material_registry)
        @warn "Empty workspace detected - forcing remesh"
        return false
    end

    # Check if mesh file exists
    if !isfile(mesh_file)
        @debug "No existing mesh file found"
        return false
    end

    # Mesh exists - can reuse
    @debug "Existing mesh found and will be reused"
    return true
end

function make_mesh!(workspace::FEMWorkspace)
    # If mesh exists and we are not forcing a remesh, do nothing and continue.
    if mesh_exists(workspace)
        @info "Using existing mesh"
        return false # Signal to continue to solver
    end

    # --- Mesh generation is required from this point on ---
    @info "Building mesh for system: $(workspace.problem_def.system.system_id)"

    try
        # Ensure Gmsh is initialized
        if gmsh.is_initialized() == 0
            gmsh.initialize()
        end

        # Perform the actual meshing
        _do_make_mesh!(workspace)
        @info "Mesh generation completed"

        # Handle mesh-only mode: preview the mesh and stop.
        # The Gmsh session is still active here.
        if workspace.opts.mesh_only
            @info "Mesh-only mode: Opening preview. Close the preview window to continue."
            preview_mesh(workspace)
            @info "Preview closed. Halting computation as per mesh_only=true."
            return true # Signal to stop computation
        end

    catch e
        @error "An error occurred during mesh generation or preview" exception = e
        rethrow(e)
    finally
        # CRITICAL: Finalize Gmsh only after all operations, including the
        # potential preview, are complete. This ensures the session is
        # always closed cleanly.
        if gmsh.is_initialized() == 1
            try
                gmsh.finalize()
            catch fin_err
                @warn "Gmsh finalization error" exception = fin_err
            end
        end
    end

    # If we are not in mesh_only mode, signal to continue to the solver.
    return false
end