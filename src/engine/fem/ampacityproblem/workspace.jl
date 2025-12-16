# """
# $(TYPEDEF)

# A container for the flattened, type-stable data arrays derived from a
# [`AmpacityProblem`](@ref). This struct serves as the primary data source
# for all subsequent computational steps.

# # Fields
# $(TYPEDFIELDS)
# """
# @kwdef struct ElectroThermalFEMWorkspace{T <: REALSCALAR}
# 	"Electro Thermal problem definition."
# 	problem_def::AmpacityProblem{T}
# 	"Formulation parameters."
# 	formulation::ElectroThermalFEMFormulation
# 	"Computation options."
# 	opts::FEMOptions

# 	"Path information."
# 	paths::Dict{Symbol, Any}

# 	"Conductor surfaces within cables."
# 	conductors::Vector{GmshObject{<:AbstractEntityData}}
# 	"Insulator surfaces within cables."
# 	insulators::Vector{GmshObject{<:AbstractEntityData}}
# 	"Domain-space physical surfaces (air and earth layers)."
# 	space_regions::Vector{GmshObject{<:AbstractEntityData}}
# 	"Domain boundary curves."
# 	boundaries::Vector{GmshObject{<:AbstractEntityData}}
# 	"Container for all pre-fragmentation entities."
# 	unassigned_entities::Dict{Vector{Float64}, AbstractEntityData}
# 	"Container for all material names used in the model."
# 	material_registry::Dict{String, Int}
# 	"Container for unique physical groups."
# 	physical_groups::Dict{Int, Material}

# 	"Vector of frequency values [Hz]."
# 	frequencies::Vector{T}
# 	"Vector of horizontal positions [m]."
# 	horz::Vector{T}
# 	"Vector of vertical positions [m]."
# 	vert::Vector{T}
# 	"Vector of internal conductor radii [m]."
# 	r_in::Vector{T}
# 	"Vector of external conductor radii [m]."
# 	r_ext::Vector{T}
# 	"Vector of internal insulator radii [m]."
# 	r_ins_in::Vector{T}
# 	"Vector of external insulator radii [m]."
# 	r_ins_ext::Vector{T}
# 	"Vector of conductor resistivities [Ω·m]."
# 	rho_cond::Vector{T}
# 	"Vector of conductor temperature coefficients [1/°C]."
# 	alpha_cond::Vector{T}
# 	"Vector of conductor relative permeabilities."
# 	mu_cond::Vector{T}
# 	"Vector of conductor relative permittivities."
# 	eps_cond::Vector{T}
# 	"Vector of insulator resistivities [Ω·m]."
# 	rho_ins::Vector{T}
# 	"Vector of insulator relative permeabilities."
# 	mu_ins::Vector{T}
# 	"Vector of insulator relative permittivities."
# 	eps_ins::Vector{T}
# 	"Vector of insulator loss tangents."
# 	tan_ins::Vector{T}
# 	"Effective earth resistivity (layers × frequencies)."
# 	rho_g::Matrix{T}
# 	"Effective earth permittivity (layers × frequencies)."
# 	eps_g::Matrix{T}
# 	"Effective earth permeability (layers × frequencies)."
# 	mu_g::Matrix{T}
# 	"Effective earth conductivity (layers × frequencies)."
# 	kappa_g::Matrix{T}
# 	"Operating temperature [°C]."
# 	temp::T
# 	"Number of frequency samples."
# 	n_frequencies::Int
# 	"Number of phases in the system."
# 	n_phases::Int
# 	"Number of cables in the system."
# 	n_cables::Int

# end



# """
# $(TYPEDSIGNATURES)

# Initializes and populates the [`ElectroThermalFEMWorkspace`](@ref) by normalizing a
# [`AmpacityProblem`](@ref) into flat, type-stable arrays.
# """
# function init_workspace(
# 	problem::AmpacityProblem{U},
# 	formulation::ElectroThermalFEMFormulation,
# ) where {U <: REALSCALAR}

# 	opts = formulation.options

# 	system = problem.system
# 	n_frequencies = length(problem.frequencies)
# 	n_phases = sum(length(cable.design_data.components) for cable in system.cables)
# 	n_cables = system.num_cables

# 	# Pre-allocate 1D arrays
# 	T = BASE_FLOAT
# 	frequencies = Vector{T}(undef, n_frequencies)
# 	horz = Vector{T}(undef, n_phases)
# 	vert = Vector{T}(undef, n_phases)
# 	r_in = Vector{T}(undef, n_phases)
# 	r_ext = Vector{T}(undef, n_phases)
# 	r_ins_in = Vector{T}(undef, n_phases)
# 	r_ins_ext = Vector{T}(undef, n_phases)
# 	rho_cond = Vector{T}(undef, n_phases)
# 	alpha_cond = Vector{T}(undef, n_phases)
# 	mu_cond = Vector{T}(undef, n_phases)
# 	eps_cond = Vector{T}(undef, n_phases)
# 	rho_ins = Vector{T}(undef, n_phases)
# 	mu_ins = Vector{T}(undef, n_phases)
# 	eps_ins = Vector{T}(undef, n_phases)
# 	tan_ins = Vector{T}(undef, n_phases)


# 	# Fill arrays, ensuring type promotion
# 	frequencies .= problem.frequencies

# 	idx = 0
# 	for (cable_idx, cable) in enumerate(system.cables)
# 		for (comp_idx, component) in enumerate(cable.design_data.components)
# 			idx += 1
# 			# Geometric properties
# 			horz[idx] = T(cable.horz)
# 			vert[idx] = T(cable.vert)
# 			r_in[idx] = T(component.conductor_group.radius_in)
# 			r_ext[idx] = T(component.conductor_group.radius_ext)
# 			r_ins_in[idx] = T(component.insulator_group.radius_in)
# 			r_ins_ext[idx] = T(component.insulator_group.radius_ext)

# 			# Material properties
# 			rho_cond[idx] = T(component.conductor_props.rho)
# 			alpha_cond[idx] = T(component.conductor_props.alpha)
# 			mu_cond[idx] = T(component.conductor_props.mu_r)
# 			eps_cond[idx] = T(component.conductor_props.eps_r)
# 			rho_ins[idx] = T(component.insulator_props.rho)
# 			mu_ins[idx] = T(component.insulator_props.mu_r)
# 			eps_ins[idx] = T(component.insulator_props.eps_r)

# 			# Calculate loss factor from resistivity
# 			ω = 2 * π * f₀  # Using default frequency
# 			C_eq = T(component.insulator_group.shunt_capacitance)
# 			G_eq = T(component.insulator_group.shunt_conductance)
# 			tan_ins[idx] = G_eq / (ω * C_eq)

# 		end
# 	end

# 	(rho_g, eps_g, mu_g, kappa_g) = _get_earth_data(
# 		nothing,
# 		problem.earth_props,
# 		frequencies,
# 		T,
# 	)

# 	temp = T(problem.temperature)

# 	# Construct and return the ElectroThermalFEMWorkspace struct
# 	return ElectroThermalFEMWorkspace{T}(
# 			problem_def = problem, formulation = formulation, opts = opts,
# 			paths = setup_paths(problem.system, formulation),
# 			conductors = Vector{GmshObject{<:AbstractEntityData}}(),
# 			insulators = Vector{GmshObject{<:AbstractEntityData}}(),
# 			space_regions = Vector{GmshObject{<:AbstractEntityData}}(),
# 			boundaries = Vector{GmshObject{<:AbstractEntityData}}(),
# 			unassigned_entities = Dict{Vector{Float64}, AbstractEntityData}(),
# 			material_registry = Dict{String, Int}(), # Maps physical group tags to materials,
# 			physical_groups = Dict{Int, Material}(),  # Initialize empty material registry
# 			frequencies = frequencies,
# 			horz = horz,
# 			vert = vert,
# 			r_in = r_in,
# 			r_ext = r_ext,
# 			r_ins_in = r_ins_in,
# 			r_ins_ext = r_ins_ext,
# 			rho_cond = rho_cond,
# 			alpha_cond = alpha_cond,
# 			mu_cond = mu_cond,
# 			eps_cond = eps_cond,
# 			rho_ins = rho_ins,
# 			mu_ins = mu_ins,
# 			eps_ins = eps_ins,
# 			tan_ins = tan_ins,
# 			rho_g = rho_g,
# 			eps_g = eps_g,
# 			mu_g = mu_g,
# 			kappa_g = kappa_g,
# 			temp = temp,
# 			n_frequencies = n_frequencies,
# 			n_phases = n_phases,
# 			n_cables = n_cables,
# 		)

# end
