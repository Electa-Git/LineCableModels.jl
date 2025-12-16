# """
# Utility functions for the FEMTools.jl module.
# These functions provide various utilities for file management, logging, etc.
# """


# """
# $(TYPEDSIGNATURES)

# Set up directory structure and file paths for a FEM simulation.

# # Arguments

# - `solver`: The [`FEMSolver`](@ref) containing the base path.
# - `cable_system`: The [`LineCableSystem`](@ref) containing the case ID.

# # Returns

# - A dictionary of paths for the simulation.

# # Examples

# ```julia
# paths = $(FUNCTIONNAME)(solver, cable_system)
# ```
# """
# function setup_paths(cable_system::LineCableSystem, formulation::ElectroThermalFEMFormulation)

#     opts = formulation.options
#     # Create base output directory if it doesn't exist
#     if !isdir(opts.save_path)
#         mkpath(opts.save_path)
#         @info "Created base output directory: $(display_path(opts.save_path))"
#     end

#     # Set up case-specific paths
#     case_id = cable_system.system_id
#     case_dir = joinpath(opts.save_path, case_id)

#     # Create case directory if needed
#     if !isdir(case_dir) && (opts.force_remesh || opts.mesh_only)
#         mkpath(case_dir)
#         @info "Created case directory: $(display_path(case_dir))"
#     end

#     # Create results directory path
#     results_dir = joinpath(case_dir, "results")

#     # Define key file paths
#     mesh_file = joinpath(case_dir, "$(case_id).msh")
#     geo_file = joinpath(case_dir, "$(case_id).geo_unrolled")
#     # data_file = joinpath(case_dir, "$(case_id)_data.geo")
#     multiphysics_file = String[]
#     for prob in formulation.analysis_type
#         push!(multiphysics_file, joinpath(case_dir, "$(case_id)_$(prob.resolution_name).pro"))
#     end


#     # Return compiled dictionary of paths
#     paths = Dict{Symbol,Any}(
#         :base_dir => opts.save_path,
#         :case_dir => case_dir,
#         :results_dir => results_dir,
#         :mesh_file => mesh_file,
#         :geo_file => geo_file,
#         :multiphysics_file => multiphysics_file
#     )

#     @debug "Paths configured: $(join(["$(k): $(v)" for (k,v) in paths], ", "))"

#     return paths
# end


# # function read_results_file(
# #     fem_formulation::Union{AbstractImpedanceFormulation,AbstractAdmittanceFormulation},
# #     workspace::FEMWorkspace;
# #     file::Union{String,Nothing}=nothing,
# # )

# #     results_path =
# #         joinpath(workspace.paths[:results_dir], lowercase(fem_formulation.resolution_name))

# #     if isnothing(file)
# #         file =
# #             fem_formulation isa AbstractImpedanceFormulation ? "Z.dat" :
# #             fem_formulation isa AbstractAdmittanceFormulation ? "Y.dat" :
# #             throw(ArgumentError("Invalid formulation type: $(typeof(fem_formulation))"))
# #     end

# #     filepath = joinpath(results_path, file)

# #     isfile(filepath) || Base.error("File not found: $filepath")

# #     # Read all lines from file
# #     lines = readlines(filepath)
# #     n_rows =
# #         sum([length(c.design_data.components) for c in workspace.core.system.cables])

# #     # Pre-allocate result matrix
# #     matrix = zeros(ComplexF64, n_rows, n_rows)

# #     # Process each line (matrix row)
# #     for (i, line) in enumerate(lines)
# #         # Parse all numbers, dropping the initial 0
# #         values = parse.(Float64, split(line))[2:end]

# #         # Fill matrix row with complex values
# #         for j in 1:n_rows
# #             idx = 2j - 1  # Index for real part
# #             matrix[i, j] = Complex(values[idx], values[idx+1])
# #         end
# #     end

# #     return matrix
# # end

# function archive_frequency_results(workspace::ElectroThermalFEMWorkspace, frequency::Float64)
#     try
#         results_dir = workspace.paths[:results_dir]
#         freq_dir =
#             joinpath(dirname(results_dir), "results_f=$(round(frequency, sigdigits=6))")

#         if isdir(results_dir)
#             mv(results_dir, freq_dir, force=true)
#             @debug "Archived results for f=$frequency Hz"
#         end

#         # Move solver files
#         for ext in [".res", ".pre"]
#             case_files = filter(f -> endswith(f, ext),
#                 readdir(workspace.paths[:case_dir], join=true))
#             for f in case_files
#                 mv(f, joinpath(freq_dir, basename(f)), force=true)
#             end
#         end
#     catch e
#         @warn "Failed to archive results for frequency $frequency Hz" exception = e
#     end
# end
