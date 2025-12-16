# @kwdef struct FEMOptions <: AbstractFormulationOptions
# 	"Build mesh only and preview (no solving)"
# 	mesh_only::Bool = false
# 	"Force mesh regeneration even if file exists"
# 	force_remesh::Bool = false
# 	"Skip user confirmation for overwriting results"
# 	force_overwrite::Bool = false
# 	"Generate field visualization outputs"
# 	plot_field_maps::Bool = true
# 	"Archive temporary files after each frequency run"
# 	keep_run_files::Bool = false
# 	"Reduce bundle conductors to equivalent single conductor"
# 	reduce_bundle::Bool = true
# 	"Eliminate grounded conductors from the system (Kron reduction)"
# 	kron_reduction::Bool = true
# 	"Enforce ideal transposition transposition/snaking"
# 	ideal_transposition::Bool = true
# 	"Temperature correction"
# 	temperature_correction::Bool = true
# 	"Base path for output files"
# 	save_path::String = joinpath(".", "fem_output")
# 	"Path to GetDP executable"
# 	getdp_executable::Union{String, Nothing} = nothing
# 	"Verbosity level"
# 	verbosity::Int = 0
# 	"Log file path"
# 	logfile::Union{String, Nothing} = nothing
# end

# # The one-line constructor to "promote" a NamedTuple
# FEMOptions(opts::NamedTuple) = FEMOptions(; opts...)



"""
$(TYPEDEF)

Abstract problem definition type for FEM simulation parameters.
This contains the physics-related parameters of the simulation.

$(TYPEDFIELDS)
"""
@kwdef struct FEMFormulation <: AbstractFormulationSet
	"Radius of the physical domain \\[m\\]."
	domain_radius::Float64
	"Outermost radius to apply the infinity transform \\[m\\]."
	domain_radius_inf::Float64
	"Elements per characteristic length for conductors \\[dimensionless\\]."
	elements_per_length_conductor::Int
	"Elements per characteristic length for insulators \\[dimensionless\\]."
	elements_per_length_insulator::Int
	"Elements per characteristic length for semiconductors \\[dimensionless\\]."
	elements_per_length_semicon::Int
	"Elements per characteristic length for interfaces \\[dimensionless\\]."
	elements_per_length_interfaces::Int
	"Points per circumference length (2π radians) \\[dimensionless\\]."
	points_per_circumference::Int
	"Analysis types to perform \\[dimensionless\\]."
	analysis_type::Union{Tuple{AbstractImpedanceFormulation, AbstractAdmittanceFormulation}, AmpacityFormulation}
	"Minimum mesh size \\[m\\]."
	mesh_size_min::Float64
	"Maximum mesh size \\[m\\]."
	mesh_size_max::Float64
	"Default mesh size \\[m\\]."
	mesh_size_default::Float64
	"Mesh transition regions for improved mesh quality"
	mesh_transitions::Vector{MeshTransition}
	"Mesh algorithm to use \\[dimensionless\\]."
	mesh_algorithm::Int
	"Maximum meshing retries and number of recursive subdivisions \\[dimensionless\\]."
	mesh_max_retries::Int
	"Materials database."
	materials::MaterialsLibrary
	"Solver options for FEM simulations."
	options::FEMOptions
end

"""
$(TYPEDSIGNATURES)

Inner constructor for `FEMFormulation` for impedance/admittance analysis.

This constructor is intended for internal use, called by wrapper functions like `FormulationSet`.
It accepts all parameters positionally and performs no default value setting or path validation.

# Arguments

- impedance: Impedance formulation to use`.
- admittance: Admittance formulation to use`.
- domain_radius: Domain radius for the simulation \\[m\\]. Default: 5.0.
- elements_per_length_conductor: Elements per scale length for conductors \\[dimensionless\\]. Default: 3.0.
- elements_per_length_insulator: Elements per scale length for insulators \\[dimensionless\\]. Default: 2.0.
- elements_per_length_semicon: Elements per scale length for semiconductors \\[dimensionless\\]. Default: 4.0.
- elements_per_length_interfaces: Elements per scale length for interfaces \\[dimensionless\\]. Default: 0.1.
- analysis_type: 
- mesh_size_min: Minimum mesh size \\[m\\]. Default: 1e-4.
- mesh_size_max: Maximum mesh size \\[m\\]. Default: 1.0.
- mesh_size_default: Default mesh size \\[m\\]. Default: domain_radius/10.
- mesh_algorithm: Mesh algorithm to use \\[dimensionless\\]. Default: 6.
- materials: Materials database. Default: MaterialsLibrary().

# Returns

- A [FEMFormulation](@ref) instance with the specified parameters.

# Examples

julia
# Create a problem definition with default parameters
formulation = $(FUNCTIONNAME)()

# Create a problem definition with custom parameters
formulation = $(FUNCTIONNAME)(
	domain_radius=10.0,
	elements_per_length_conductor=5.0,
	mesh_algorithm=2
)

"""
function FEMFormulation(impedance::AbstractImpedanceFormulation,
	admittance::AbstractAdmittanceFormulation,
	domain_radius::Float64,
	domain_radius_inf::Float64,
	elements_per_length_conductor::Int,
	elements_per_length_insulator::Int,
	elements_per_length_semicon::Int,
	elements_per_length_interfaces::Int,
	points_per_circumference::Int,
	mesh_size_min::Float64,
	mesh_size_max::Float64,
	mesh_size_default::Float64,
	mesh_transitions::Vector{MeshTransition},
	mesh_algorithm::Int,
	mesh_max_retries::Int,
	materials::MaterialsLibrary,
	options::FEMOptions,
)

	return new(
		domain_radius, domain_radius_inf,
		elements_per_length_conductor, elements_per_length_insulator,
		elements_per_length_semicon, elements_per_length_interfaces,
		points_per_circumference, (impedance, admittance),
		mesh_size_min, mesh_size_max, mesh_size_default,
		mesh_transitions, mesh_algorithm, mesh_max_retries, materials,
		options,
	)
end
"""
$(TYPEDSIGNATURES)

Inner constructor for `FEMFormulation` for ampacity analysis.

This constructor is intended for internal use, called by wrapper functions like `FormulationSet`.
It accepts all parameters positionally, resolves the GetDP solver path, and builds the final `FEMOptions`.

# Arguments
- `analysis_spec`: Ampacity formulation (e.g., `MagnetoThermal()`).
- `domain_radius`: Domain radius for the simulation [m].
- `domain_radius_inf`: Domain radius for the infinite element transformation [m].
- `elements_per_length_conductor`: Elements per scale length for conductors.
- `elements_per_length_insulator`: Elements per scale length for insulators.
- `elements_per_length_semicon`: Elements per scale length for semiconductors.
- `elements_per_length_interfaces`: Elements per scale length for interfaces.
- `points_per_circumference`: Number of points to discretize circular boundaries.
- `mesh_size_min`: Minimum mesh size [m].
- `mesh_size_max`: Maximum mesh size [m].
- `mesh_size_default`: Default mesh size [m].
- `mesh_transitions`: Vector of `MeshTransition` objects for mesh refinement.
- `mesh_algorithm`: Mesh algorithm to use.
- `mesh_max_retries`: Maximum number of mesh generation retries.
- `materials`: `MaterialsLibrary` instance.
- `options`: `FEMOptions` instance.

# Returns
- A new `FEMFormulation` instance.
"""
function FEMFormulation(analysis_spec::AmpacityFormulation, # analysis specification
	domain_radius::Float64,
	domain_radius_inf::Float64,
	elements_per_length_conductor::Int,
	elements_per_length_insulator::Int,
	elements_per_length_semicon::Int,
	elements_per_length_interfaces::Int,
	points_per_circumference::Int,
	mesh_size_min::Float64,
	mesh_size_max::Float64,
	mesh_size_default::Float64,
	mesh_transitions::Vector{MeshTransition},
	mesh_algorithm::Int,
	mesh_max_retries::Int,
	materials::MaterialsLibrary,
	options::FEMOptions,
)

	return new(
		analysis_type = analysis_spec,
		domain_radius, domain_radius_inf,
		elements_per_length_conductor, elements_per_length_insulator,
		elements_per_length_semicon, elements_per_length_interfaces,
		points_per_circumference, analysis_spec,
		mesh_size_min, mesh_size_max, mesh_size_default,
		mesh_transitions, mesh_algorithm, mesh_max_retries, materials,
		options = fem_opts,
	)
end


"""
$(TYPEDSIGNATURES)

Constructs a `FEMFormulation` for impedance and admittance analysis.

This function acts as a user-facing constructor, providing default values for
common simulation parameters.

# Arguments
- `impedance`: Impedance formulation to use. Default: `Darwin()`.
- `admittance`: Admittance formulation to use. Default: `Electrodynamics()`.
- `domain_radius`: Domain radius for the simulation [m]. Default: `5.0`.
- `domain_radius_inf`: Domain radius for the infinite element transformation [m]. Default: `6.25`.
- `elements_per_length_conductor`: Elements per scale length for conductors. Default: `3`.
- `elements_per_length_insulator`: Elements per scale length for insulators. Default: `2`.
- `elements_per_length_semicon`: Elements per scale length for semiconductors. Default: `4`.
- `elements_per_length_interfaces`: Elements per scale length for interfaces. Default: `3`.
- `points_per_circumference`: Number of points to discretize circular boundaries. Default: `16`.
- `mesh_size_min`: Minimum mesh size [m]. Default: `1e-4`.
- `mesh_size_max`: Maximum mesh size [m]. Default: `1.0`.
- `mesh_size_default`: Default mesh size [m]. Default: `domain_radius / 10`.
- `mesh_transitions`: Vector of `MeshTransition` objects for mesh refinement. Default: `MeshTransition[]`.
- `mesh_algorithm`: Mesh algorithm to use. Default: `5`.
- `mesh_max_retries`: Maximum number of mesh generation retries. Default: `20`.
- `materials`: Materials database. Default: `MaterialsLibrary()`.
- `options`: `NamedTuple` of options to be converted to `FEMOptions`. Default: `(; )`.

# Returns
- A `FEMFormulation` instance configured for FEM analysis.

# Examples
```julia
# Create a FEM formulation with default parameters
formulation_fem = $(FUNCTIONNAME)(Val(:FEM))

# Create a FEM formulation with custom parameters
formulation_fem_custom = $(FUNCTIONNAME)(
    Val(:FEM);
    domain_radius=10.0,
    elements_per_length_conductor=5,
    mesh_algorithm=2
)
```
"""
function FormulationSet(::Val{:FEM}; impedance::AbstractImpedanceFormulation = Darwin(),
	admittance::AbstractAdmittanceFormulation = Electrodynamics(),
	domain_radius::Float64 = 5.0,
	domain_radius_inf::Float64 = 6.25,
	elements_per_length_conductor::Int = 3,
	elements_per_length_insulator::Int = 2,
	elements_per_length_semicon::Int = 4,
	elements_per_length_interfaces::Int = 3,
	points_per_circumference::Int = 16,
	mesh_size_min::Float64 = 1e-4,
	mesh_size_max::Float64 = 1.0,
	mesh_size_default::Float64 = domain_radius / 10,
	mesh_transitions::Vector{MeshTransition} = MeshTransition[],
	mesh_algorithm::Int = 5,
	mesh_max_retries::Int = 20,
	materials::MaterialsLibrary = MaterialsLibrary(),
	options = (;),
)
	# Resolve solver path
	validated_path = _resolve_getdp_path(options)

	# Create a new NamedTuple with the validated path overwriting any user value
	final_opts = merge(options, (getdp_executable = validated_path,))
	fem_opts = build_options(FEMOptions, final_opts; strict = true)

	return FEMFormulation(
		analysis_type = (impedance, admittance),
		domain_radius = domain_radius,
		domain_radius_inf = domain_radius_inf,
		elements_per_length_conductor = elements_per_length_conductor,
		elements_per_length_insulator = elements_per_length_insulator,
		elements_per_length_semicon = elements_per_length_semicon,
		elements_per_length_interfaces = elements_per_length_interfaces,
		points_per_circumference = points_per_circumference,
		mesh_size_min = mesh_size_min,
		mesh_size_max = mesh_size_max,
		mesh_size_default = mesh_size_default,
		mesh_transitions = mesh_transitions,
		mesh_algorithm = mesh_algorithm,
		mesh_max_retries = mesh_max_retries,
		materials = materials,
		options = fem_opts,
	)
end


"""
$(TYPEDSIGNATURES)

Constructs a `FEMFormulation` for ampacity analysis.

This function acts as a user-facing constructor, providing default values for
common simulation parameters.

# Arguments
- `analysis_type`: Ampacity formulation to use. Default: `MagnetoThermal()`.
- `domain_radius`: Domain radius for the simulation [m]. Default: `5.0`.
- `domain_radius_inf`: Domain radius for the infinite element transformation [m]. Default: `6.25`.
- `elements_per_length_conductor`: Elements per scale length for conductors. Default: `3`.
- `elements_per_length_insulator`: Elements per scale length for insulators. Default: `2`.
- `elements_per_length_semicon`: Elements per scale length for semiconductors. Default: `4`.
- `elements_per_length_interfaces`: Elements per scale length for interfaces. Default: `3`.
- `points_per_circumference`: Number of points to discretize circular boundaries. Default: `16`.
- `mesh_size_min`: Minimum mesh size [m]. Default: `1e-4`.
- `mesh_size_max`: Maximum mesh size [m]. Default: `1.0`.
- `mesh_size_default`: Default mesh size [m]. Default: `domain_radius / 10`.
- `mesh_transitions`: Vector of `MeshTransition` objects for mesh refinement. Default: `MeshTransition[]`.
- `mesh_algorithm`: Mesh algorithm to use. Default: `5`.
- `mesh_max_retries`: Maximum number of mesh generation retries. Default: `20`.
- `materials`: Materials database. Default: `MaterialsLibrary()`.
- `options`: `NamedTuple` of options to be converted to `FEMOptions`. Default: `(; )`.

# Returns
- A `FEMFormulation` instance configured for Ampacity analysis.

# Examples
```julia
# Create an Ampacity formulation with default parameters
formulation_amp = $(FUNCTIONNAME)(Val(:Ampacity))

# Create an Ampacity formulation with custom parameters
formulation_amp_custom = $(FUNCTIONNAME)(
    Val(:Ampacity);
    domain_radius=10.0,
    elements_per_length_conductor=5,
    analysis_type=MagnetoThermal()
)
```
"""
function FormulationSet(::Val{:Ampacity};
	analysis_type::AmpacityFormulation = MagnetoThermal(),
	domain_radius::Float64 = 5.0,
	domain_radius_inf::Float64 = 6.25,
	elements_per_length_conductor::Int = 3,
	elements_per_length_insulator::Int = 2,
	elements_per_length_semicon::Int = 4,
	elements_per_length_interfaces::Int = 3,
	points_per_circumference::Int = 16,
	mesh_size_min::Float64 = 1e-4,
	mesh_size_max::Float64 = 1.0,
	mesh_size_default::Float64 = domain_radius / 10,
	mesh_transitions::Vector{MeshTransition} = MeshTransition[],
	mesh_algorithm::Int = 5,
	mesh_max_retries::Int = 20,
	materials::MaterialsLibrary = MaterialsLibrary(),
	options = (;),
)
	# Resolve solver path
	validated_path = _resolve_getdp_path(options)

	# Create a new NamedTuple with the validated path overwriting any user value
	final_opts = merge(options, (getdp_executable = validated_path,))
	fem_opts = build_options(FEMOptions, final_opts; strict = true)

	return FEMFormulation(analysis_type = analysis_type,
		domain_radius = domain_radius,
		domain_radius_inf = domain_radius_inf,
		elements_per_length_conductor = elements_per_length_conductor,
		elements_per_length_insulator = elements_per_length_insulator,
		elements_per_length_semicon = elements_per_length_semicon,
		elements_per_length_interfaces = elements_per_length_interfaces,
		points_per_circumference = points_per_circumference,
		mesh_size_min = mesh_size_min,
		mesh_size_max = mesh_size_max,
		mesh_size_default = mesh_size_default,
		mesh_transitions = mesh_transitions,
		mesh_algorithm = mesh_algorithm,
		mesh_max_retries = mesh_max_retries,
		materials = materials,
		options = fem_opts,
	)
end

