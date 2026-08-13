"""
Visualization functions for the FEMTools.jl module.
These functions handle the visualization of the mesh and results.
"""

"""
$(TYPEDSIGNATURES)

Preview the mesh in the Gmsh GUI.

# Arguments

- `workspace`: The [`FEMWorkspace`](@ref) containing the model.

# Returns

- Nothing. Launches the Gmsh GUI.

# Examples

```julia
$(FUNCTIONNAME)(workspace)
```
"""
function preview_mesh(workspace::FEMWorkspace)

    if gmsh.is_initialized() == 0
        gmsh.initialize()
        @debug "Initialized Gmsh for mesh preview"
    else
        @debug "Gmsh already initialized"
    end

    try
        # Set visualization options
        gmsh.option.set_number("Geometry.SurfaceLabels", 0)  # Show surface labels
        gmsh.option.set_number("Geometry.PointNumbers", 0)
        gmsh.option.set_number("Geometry.CurveNumbers", 0)
        gmsh.option.set_number("Geometry.SurfaceNumbers", 0)
        gmsh.option.set_number("Geometry.NumSubEdges", 160)
        gmsh.option.set_number("Geometry.Points", 1)
        gmsh.option.set_number("Geometry.Curves", 1)
        gmsh.option.set_number("Geometry.Surfaces", 0)
        gmsh.option.set_number("Mesh.ColorCarousel", 2)  # Colors by physical group
        gmsh.option.set_number("Mesh.LineWidth", 1)
        gmsh.option.set_number("Mesh.SurfaceFaces", 1)

        # Initialize FLTK GUI
        gmsh.fltk.initialize()

        @info "Launching Gmsh GUI for mesh preview"
        @info "Close the Gmsh window to continue..."

        # Define event check function
        function check_for_event()
            action = gmsh.onelab.get_string("ONELAB/Action")
            if length(action) > 0 && action[1] == "check"
                gmsh.onelab.set_string("ONELAB/Action", [""])
                @debug "UI interaction detected"
                gmsh.graphics.draw()
            end
            return true
        end

        # Wait for user to close the window
        while gmsh.fltk.is_available() == 1 && check_for_event()
            gmsh.fltk.wait()
        end

        @info "Mesh preview closed"

    catch e
        @warn "Error during mesh preview: $e"
    end
end

"""
$(TYPEDSIGNATURES)

Preview a single electromagnetic field result file in Gmsh GUI.

# Arguments
- `workspace`: The [`FEMWorkspace`](@ref) containing the model.
- `pos_file`: Path to the .pos file to visualize.

# Examples
```julia
$(FUNCTIONNAME)(workspace, "/path/to/result.pos")
```
"""
function preview_results(workspace::FEMWorkspace, pos_file::String)
    # Validate inputs
    if !isfile(pos_file)
        @error "Result file not found: $pos_file"
        return
    end

    if !endswith(pos_file, ".pos")
        @error "File must be a .pos file: $pos_file"
        return
    end

    # Initialize Gmsh
    gmsh.initialize()

    try
        # Add single model
        gmsh.model.add("field_view")

        # Merge mesh file
        mesh_file = workspace.paths[:mesh_file]
        if isfile(mesh_file)
            gmsh.merge(abspath(mesh_file))
        else
            @error "Mesh file not found: $mesh_file"
            return
        end

        # Merge the single result file
        @info "Loading field data: $(display_path(pos_file))"
        gmsh.merge(abspath(pos_file))

        # Set mesh color to light gray
        gmsh.option.set_color("Mesh.Color.Lines", 240, 240, 240)
        gmsh.option.set_number("Mesh.ColorCarousel", 0)
        gmsh.option.set_number("Mesh.LineWidth", 1)
        gmsh.option.set_number("Mesh.SurfaceFaces", 0)
        gmsh.option.set_number("Mesh.Lines", 1)
        gmsh.option.set_number("Geometry.Points", 0)
        gmsh.option.set_number("General.InitialModule", 4)

        # Get view tags and configure
        view_tags = gmsh.view.getTags()

        if isempty(view_tags)
            @warn "No field views found in file"
            return
        end

        # Configure field visualization
        for view_tag in view_tags
            gmsh.view.option.set_number(view_tag, "IntervalsType", 2)
            gmsh.view.option.set_number(view_tag, "RangeType", 3)
            gmsh.view.option.set_number(view_tag, "ShowTime", 0)
        end

        @info "Launching Gmsh GUI with $(length(view_tags)) field view(s)"
        @info "Close the Gmsh window to continue..."

        # Launch GUI
        gmsh.fltk.run()

        @info "Field visualization closed"

    catch e
        @error "Error during field visualization" exception = e
    finally
        gmsh.finalize()
    end
end