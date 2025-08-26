"""
$(TYPEDEF)

FEMWorkspace - The central workspace for FEM simulations.
This is the main container that maintains all state during the simulation process.

$(TYPEDFIELDS)
"""
mutable struct FEMWorkspace
    "Line parameters problem definition."
    problem_def::LineParametersProblem
    "Formulation parameters."
    formulation::FEMFormulation
    "Computation options."
    opts::FEMOptions

    "Path information."
    paths::Dict{Symbol,String}

    "Conductor surfaces within cables."
    conductors::Vector{GmshObject{<:AbstractEntityData}}
    "Insulator surfaces within cables."
    insulators::Vector{GmshObject{<:AbstractEntityData}}
    "Domain-space physical surfaces (air and earth layers)."
    space_regions::Vector{GmshObject{<:AbstractEntityData}}
    "Domain boundary curves."
    boundaries::Vector{GmshObject{<:AbstractEntityData}}
    "Container for all pre-fragmentation entities."
    unassigned_entities::Dict{Vector{Float64},AbstractEntityData}
    "Container for all material names used in the model."
    material_registry::Dict{String,Int}
    "Container for unique physical groups."
    physical_groups::Dict{Int,Material}

    """
    $(TYPEDSIGNATURES)

    Constructs a [`FEMWorkspace`](@ref) instance.

    # Arguments

    - `cable_system`: Cable system being simulated.
    - `formulation`: Problem definition parameters.
    - `solver`: Solver parameters.
    - `frequency`: Simulation frequency \\[Hz\\]. Default: 50.0.

    # Returns

    - A [`FEMWorkspace`](@ref) instance with the specified parameters.

    # Examples

    ```julia
    # Create a workspace
    workspace = $(FUNCTIONNAME)(cable_system, formulation, solver)
    ```
    """
    function FEMWorkspace(problem::LineParametersProblem, formulation::FEMFormulation)

        # Initialize empty workspace
        opts = formulation.options
        workspace = new(
            problem, formulation, opts,
            Dict{Symbol,String}(), # Path information.
            Vector{GmshObject{<:AbstractEntityData}}(), #conductors
            Vector{GmshObject{<:AbstractEntityData}}(), #insulators
            Vector{GmshObject{<:AbstractEntityData}}(), #space_regions
            Vector{GmshObject{<:AbstractEntityData}}(), #boundaries
            Dict{Vector{Float64},AbstractEntityData}(), #unassigned_entities
            Dict{String,Int}(),  # Initialize empty material registry
            Dict{Int,Material}(), # Maps physical group tags to materials
        )

        # Set up paths
        workspace.paths = setup_paths(problem.system, formulation)

        return workspace
    end
end

function init_workspace(problem, formulation, workspace)
    if isnothing(workspace)
        @debug "Creating new workspace"
        workspace = FEMWorkspace(problem, formulation)
    else
        @debug "Reusing existing workspace"
    end

    opts = formulation.options

    setup_logging!(opts.verbosity, opts.logfile)

    # Handle existing results - check both current and archived
    results_dir = workspace.paths[:results_dir]
    base_dir = dirname(results_dir)

    # Check current results directory
    current_results_exist = isdir(results_dir) && !isempty(readdir(results_dir))

    # Check for archived frequency results (results_f* pattern)
    archived_results_exist = false
    if isdir(base_dir)
        archived_dirs = filter(d -> startswith(d, "results_f") && isdir(joinpath(base_dir, d)),
            readdir(base_dir))
        archived_results_exist = !isempty(archived_dirs)
    end

    # Handle existing results if any are found
    if current_results_exist || archived_results_exist
        if opts.force_overwrite
            # Remove both current and archived results
            if current_results_exist
                rm(results_dir, recursive=true, force=true)
            end
            if archived_results_exist
                for archived_dir in archived_dirs
                    rm(joinpath(base_dir, archived_dir), recursive=true, force=true)
                end
                @debug "Removed $(length(archived_dirs)) archived result directories"
            end
        else
            # Build informative error message
            error_msg = "Existing results found:\n"
            if current_results_exist
                error_msg *= "  - Current results: $results_dir\n"
            end
            if archived_results_exist
                error_msg *= "  - Archived results: $(length(archived_dirs)) frequency directories\n"
            end
            error_msg *= "Set force_overwrite=true to automatically delete existing results."

            Base.error(error_msg)
        end
    end

    return workspace
end