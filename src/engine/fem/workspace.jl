import ..Engine: _get_earth_data

"""
$(TYPEDEF)

FEMWorkspace - The central workspace for FEM simulations.
This is the main container that maintains all state during the simulation process.

$(TYPEDFIELDS)
"""
struct FEMWorkspace{T <: AbstractFloat}
	"Line parameters problem definition."
	problem_def::LineParametersProblem
	"Formulation parameters."
	formulation::FEMFormulation
	"Computation options."
	opts::FEMOptions

	"Path information."
	paths::Dict{Symbol, String}

	"Conductor surfaces within cables."
	conductors::Vector{GmshObject{<:AbstractEntityData}}
	"Insulator surfaces within cables."
	insulators::Vector{GmshObject{<:AbstractEntityData}}
	"Domain-space physical surfaces (air and earth layers)."
	space_regions::Vector{GmshObject{<:AbstractEntityData}}
	"Domain boundary curves."
	boundaries::Vector{GmshObject{<:AbstractEntityData}}
	"Container for all pre-fragmentation entities."
	unassigned_entities::Dict{Vector{Float64}, AbstractEntityData}
	"Container for all material names used in the model."
	material_registry::Dict{String, Int}
	"Container for unique physical groups."
	physical_groups::Dict{Int, Material}

	"Vector of frequency values [Hz]."
	freq::Vector{T}
	"Vector of horizontal positions [m]."
	horz::Vector{T}
	"Vector of vertical positions [m]."
	vert::Vector{T}
	"Vector of internal conductor radii [m]."
	r_in::Vector{T}
	"Vector of external conductor radii [m]."
	r_ext::Vector{T}
	"Vector of internal insulator radii [m]."
	r_ins_in::Vector{T}
	"Vector of external insulator radii [m]."
	r_ins_ext::Vector{T}
	"Vector of conductor resistivities [Ω·m]."
	rho_cond::Vector{T}
	"Vector of conductor temperature coefficients [1/°C]."
	alpha_cond::Vector{T}
	"Vector of conductor relative permeabilities."
	mu_cond::Vector{T}
	"Vector of conductor relative permittivities."
	eps_cond::Vector{T}
	"Vector of insulator resistivities [Ω·m]."
	rho_ins::Vector{T}
	"Vector of insulator relative permeabilities."
	mu_ins::Vector{T}
	"Vector of insulator relative permittivities."
	eps_ins::Vector{T}
	"Vector of insulator loss tangents."
	tan_ins::Vector{T}
	"Vector of phase mapping indices."
	phase_map::Vector{Int}
	"Vector of cable mapping indices."
	cable_map::Vector{Int}
	"Effective earth resistivity (layers × freq)."
	rho_g::Matrix{T}
	"Effective earth permittivity (layers × freq)."
	eps_g::Matrix{T}
	"Effective earth permeability (layers × freq)."
	mu_g::Matrix{T}
	"Operating temperature [°C]."
	temp::T
	"Number of frequency samples."
	n_frequencies::Int
	"Number of phases in the system."
	n_phases::Int
	"Number of cables in the system."
	n_cables::Int
	"Full component-based Z matrix (before bundling/reduction)."
	Z::Array{Complex{T}, 3}
	"Full component-based Y matrix (before bundling/reduction)."
	Y::Array{Complex{T}, 3}


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
	function FEMWorkspace(
		problem::LineParametersProblem{U},
		formulation::FEMFormulation,
	) where {U <: REALSCALAR}

		# Initialize empty workspace
		opts = formulation.options

		system = problem.system
		n_frequencies = length(problem.frequencies)
		n_phases = sum(length(cable.design_data.components) for cable in system.cables)

		# Pre-allocate 1D arrays
		T = BASE_FLOAT
		freq = Vector{T}(undef, n_frequencies)
		horz = Vector{T}(undef, n_phases)
		vert = Vector{T}(undef, n_phases)
		r_in = Vector{T}(undef, n_phases)
		r_ext = Vector{T}(undef, n_phases)
		r_ins_in = Vector{T}(undef, n_phases)
		r_ins_ext = Vector{T}(undef, n_phases)
		rho_cond = Vector{T}(undef, n_phases)
		alpha_cond = Vector{T}(undef, n_phases)
		mu_cond = Vector{T}(undef, n_phases)
		eps_cond = Vector{T}(undef, n_phases)
		rho_ins = Vector{T}(undef, n_phases)
		mu_ins = Vector{T}(undef, n_phases)
		eps_ins = Vector{T}(undef, n_phases)
		tan_ins = Vector{T}(undef, n_phases)   # Loss tangent for insulator
		phase_map = Vector{Int}(undef, n_phases)
		cable_map = Vector{Int}(undef, n_phases)
		Z = zeros(Complex{T}, n_phases, n_phases, n_frequencies)
		Y = zeros(Complex{T}, n_phases, n_phases, n_frequencies)

		# Fill arrays, ensuring type promotion
		freq .= to_nominal.(problem.frequencies)

		idx = 0
		for (cable_idx, cable) in enumerate(system.cables)
			for (comp_idx, component) in enumerate(cable.design_data.components)
				idx += 1
				# Geometric properties
				horz[idx] = to_nominal(cable.horz)
				vert[idx] = to_nominal(cable.vert)
				r_in[idx] = to_nominal(component.conductor_group.r_in)
				r_ext[idx] = to_nominal(component.conductor_group.r_ex)
				r_ins_in[idx] = to_nominal(component.insulator_group.r_in)
				r_ins_ext[idx] = to_nominal(component.insulator_group.r_ex)

				# Material properties
				rho_cond[idx] = to_nominal(component.conductor_props.rho)
				alpha_cond[idx] = to_nominal(component.conductor_props.alpha)
				mu_cond[idx] = to_nominal(component.conductor_props.mu_r)
				eps_cond[idx] = to_nominal(component.conductor_props.eps_r)
				rho_ins[idx] = to_nominal(component.insulator_props.rho)
				mu_ins[idx] = to_nominal(component.insulator_props.mu_r)
				eps_ins[idx] = to_nominal(component.insulator_props.eps_r)

				# Calculate loss factor from resistivity
				ω = 2 * π * f₀  # Using default frequency
				C_eq = to_nominal(component.insulator_group.shunt_capacitance)
				G_eq = to_nominal(component.insulator_group.shunt_conductance)
				tan_ins[idx] = G_eq / (ω * C_eq)

				# Mapping
				phase_map[idx] = cable.conn[comp_idx]
				cable_map[idx] = cable_idx
			end
		end

		(rho_g, eps_g, mu_g) = _get_earth_data(
			nothing,
			problem.earth_props,
			freq,
			T,
		)


		temp = to_nominal(problem.temperature)


		workspace = new{T}(
			problem, formulation, opts,
			setup_paths(problem.system, formulation),
			# Dict{Symbol,String}(), # Path information.
			Vector{GmshObject{<:AbstractEntityData}}(), #conductors
			Vector{GmshObject{<:AbstractEntityData}}(), #insulators
			Vector{GmshObject{<:AbstractEntityData}}(), #space_regions
			Vector{GmshObject{<:AbstractEntityData}}(), #boundaries
			Dict{Vector{Float64}, AbstractEntityData}(), #unassigned_entities
			Dict{String, Int}(),  # Initialize empty material registry
			Dict{Int, Material}(), # Maps physical group tags to materials,
			freq,
			horz, vert,
			r_in, r_ext,
			r_ins_in, r_ins_ext,
			rho_cond, alpha_cond, mu_cond, eps_cond,
			rho_ins, mu_ins, eps_ins, tan_ins,
			phase_map, cable_map, rho_g,
			eps_g, mu_g,
			temp, n_frequencies, n_phases,
			system.num_cables, Z, Y,
		)

		# Set up paths
		# workspace.paths = setup_paths(problem.system, formulation)

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

	# set_verbosity!(opts.verbosity, opts.logfile)

	# Handle existing results - check both current and archived
	results_dir = workspace.paths[:results_dir]
	base_dir = dirname(results_dir)

	# Check current results directory
	current_results_exist = isdir(results_dir) && !isempty(readdir(results_dir))

	# Check for archived frequency results (results_f* pattern)
	archived_results_exist = false
	if isdir(base_dir)
		archived_dirs =
			filter(d -> startswith(d, "results_f") && isdir(joinpath(base_dir, d)),
				readdir(base_dir))
		archived_results_exist = !isempty(archived_dirs)
	end

	# Handle existing results if any are found
	if current_results_exist || archived_results_exist
		if opts.force_overwrite
			# Remove both current and archived results
			if current_results_exist
				rm(results_dir, recursive = true, force = true)
			end
			if archived_results_exist
				for archived_dir in archived_dirs
					rm(joinpath(base_dir, archived_dir), recursive = true, force = true)
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
