"""
Utility functions for the FEMTools.jl module.
These functions provide various utilities for file management, logging, etc.
"""


"""
$(TYPEDSIGNATURES)

Set up directory structure and file paths for a FEM simulation.

# Arguments

- `solver`: The [`FEMSolver`](@ref) containing the base path.
- `cable_system`: The [`LineCableSystem`](@ref) containing the case ID.

# Returns

- A dictionary of paths for the simulation.

# Examples

```julia
paths = $(FUNCTIONNAME)(solver, cable_system)
```
"""
function setup_paths(cable_system::LineCableSystem, formulation::FEMFormulation)

    opts = formulation.options
    # Create base output directory if it doesn't exist
    if !isdir(opts.save_path)
        mkpath(opts.save_path)
        @info "Created base output directory: $(display_path(opts.save_path))"
    end

    # Set up case-specific paths
    case_id = cable_system.system_id
    case_dir = joinpath(opts.save_path, case_id)

    # Create case directory if needed
    if !isdir(case_dir) && (opts.force_remesh || opts.mesh_only)
        mkpath(case_dir)
        @info "Created case directory: $(display_path(case_dir))"
    end

    # Create results directory path
    results_dir = joinpath(case_dir, "results")

    # Define key file paths
    mesh_file = joinpath(case_dir, "$(case_id).msh")
    geo_file = joinpath(case_dir, "$(case_id).geo_unrolled")
    # data_file = joinpath(case_dir, "$(case_id)_data.geo")

    impedance_res = lowercase(formulation.analysis_type[1].resolution_name)
    impedance_file = joinpath(case_dir, "$(case_id)_$(impedance_res).pro")

    admittance_res = lowercase(formulation.analysis_type[2].resolution_name)
    admittance_file = joinpath(case_dir, "$(case_id)_$(admittance_res).pro")

    # Return compiled dictionary of paths
    paths = Dict{Symbol,String}(
        :base_dir => opts.save_path,
        :case_dir => case_dir,
        :results_dir => results_dir,
        :mesh_file => mesh_file,
        :geo_file => geo_file,
        :impedance_file => impedance_file,
        :admittance_file => admittance_file,
    )

    @debug "Paths configured: $(join(["$(k): $(v)" for (k,v) in paths], ", "))"

    return paths
end

"""
$(TYPEDSIGNATURES)

Clean up files based on configuration flags.

# Arguments

- `paths`: Dictionary of paths for the simulation.
- `solver`: The [`FEMSolver`](@ref) containing the configuration flags.

# Returns

- Nothing. Deletes files as specified by the configuration.

# Examples

```julia
$(FUNCTIONNAME)(paths, solver)
```
"""
function cleanup_files(paths::Dict{Symbol,String}, opts::NamedTuple)
    if opts.force_remesh
        # If force_remesh is true, delete mesh-related files
        if isfile(paths[:mesh_file])
            rm(paths[:mesh_file], force=true)
            @info "Removed existing mesh file: $(display_path(paths[:mesh_file]))"
        end

        if isfile(paths[:geo_file])
            rm(paths[:geo_file], force=true)
            @info "Removed existing geometry file: $(display_path(paths[:geo_file]))"
        end
    end

    if opts.overwrite_results && opts.run_solver

        # Add cleanup for .pro files in case_dir
        for file in readdir(paths[:case_dir])
            if endswith(file, ".pro")
                filepath = joinpath(paths[:case_dir], file)
                rm(filepath, force=true)
                @info "Removed existing problem file: $(display_path(filepath))"
            end
        end

        # If overwriting results and running solver, clear the results directory
        if isdir(paths[:results_dir])
            for file in readdir(paths[:results_dir])
                filepath = joinpath(paths[:results_dir], file)
                if isfile(filepath)
                    rm(filepath, force=true)
                end
            end
            @info "Cleared existing results in: $(display_path(paths[:results_dir]))"
        end
    end
end

function read_results_file(fem_formulation::Union{AbstractImpedanceFormulation,AbstractAdmittanceFormulation}, workspace::FEMWorkspace; file::Union{String,Nothing}=nothing)

    results_path = joinpath(workspace.paths[:results_dir], lowercase(fem_formulation.resolution_name))

    if isnothing(file)
        file = fem_formulation isa AbstractImpedanceFormulation ? "Z.dat" :
               fem_formulation isa AbstractAdmittanceFormulation ? "Y.dat" :
               throw(ArgumentError("Invalid formulation type: $(typeof(fem_formulation))"))
    end

    filepath = joinpath(results_path, file)

    isfile(filepath) || Base.error("File not found: $filepath")

    # Read all lines from file
    lines = readlines(filepath)
    n_rows = sum([length(c.design_data.components) for c in workspace.problem_def.system.cables])

    # Pre-allocate result matrix
    matrix = zeros(ComplexF64, n_rows, n_rows)

    # Process each line (matrix row)
    for (i, line) in enumerate(lines)
        # Parse all numbers, dropping the initial 0
        values = parse.(Float64, split(line))[2:end]

        # Fill matrix row with complex values
        for j in 1:n_rows
            idx = 2j - 1  # Index for real part
            matrix[i, j] = Complex(values[idx], values[idx+1])
        end
    end

    return matrix
end


# Verbosity Levels in GetDP
# Level	Output Description
# 0	     Silent (no output)
# 1	     Errors only
# 2	     Errors + warnings
# 3	     Errors + warnings + basic info
# 4	     Detailed debugging
# 5	     Full internal tracing
function map_verbosity_to_getdp(verbosity::Int)
    if is_headless()       # Prevent huge logs in CI/CD deploys
        @info "Running in headless mode, suppressing GetDP output"
        return 0            # Gmsh Silent level
    elseif verbosity >= 2   # Debug
        return 4            # GetDP Debug level
    elseif verbosity == 1   # Info
        return 3            # GetDP Info level
    else                    # Warn
        return 1            # GetDP Errors level
    end
end

# Verbosity Levels in Gmsh
# Level  Output Description
# 0      Silent (no output)
# 1      Errors only
# 2      Warnings
# 3      Direct/Important info
# 4      Information
# 5      Status messages
# 99     Debug
function map_verbosity_to_gmsh(verbosity::Int)
    if is_headless()       # Prevent huge logs in CI/CD deploys
        @info "Running in headless mode, suppressing Gmsh output"
        return 0            # Gmsh Silent level
    elseif verbosity >= 2   # Debug
        return 99           # Gmsh Debug level
    elseif verbosity == 1   # Info
        return 4            # Gmsh Information level
    else                    # Warn
        return 1            # Gmsh Errors level
    end
end

function calc_domain_size(earth_params::EarthModel, f::Vector{<:Float64}; min_radius=5.0, max_radius=5000.0)
    # Find the earth layer with the highest resistivity to determine the domain size
    if isempty(earth_params.layers)
        error("EarthModel has no layers defined.")
    end

    # Find the index of the layer with the maximum resistivity at the first frequency
    max_rho_idx = argmax([layer.rho_g[1] for layer in earth_params.layers])
    target_layer = earth_params.layers[max_rho_idx]

    rho_g = target_layer.rho_g[1]
    mu_g = target_layer.mu_g[1]
    freq = first(f) # Use the first frequency for the calculation
    skin_depth_earth = abs(sqrt(rho_g / (1im * 2 * pi * freq * mu_g)))
    return clamp(skin_depth_earth, min_radius, max_radius)
end

function archive_frequency_results(workspace::FEMWorkspace, frequency::Float64)
    try
        results_dir = workspace.paths[:results_dir]
        freq_dir = joinpath(dirname(results_dir), "results_f=$(round(frequency, sigdigits=6))")

        if isdir(results_dir)
            mv(results_dir, freq_dir, force=true)
            @debug "Archived results for f=$frequency Hz"
        end

        # Move solver files
        for ext in [".res", ".pre"]
            case_files = filter(f -> endswith(f, ext),
                readdir(workspace.paths[:case_dir], join=true))
            for f in case_files
                mv(f, joinpath(freq_dir, basename(f)), force=true)
            end
        end
    catch e
        @warn "Failed to archive results for frequency $frequency Hz" exception = e
    end
end

# Internal helper to find the executable path 
function _resolve_getdp_path(opts::NamedTuple)
    user_path = get(opts, :getdp_executable, nothing)
    if user_path isa String && isfile(user_path)
        return user_path
    end
    fallback_path = GetDP.get_getdp_executable()
    if isfile(fallback_path)
        return fallback_path
    end
    error("GetDP executable not found.")
end