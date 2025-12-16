import ..Engine: _get_earth_data

"""
$(TYPEDEF)

FEMWorkspaceCore - The core workspace for FEM simulations.

$(TYPEDFIELDS)
"""
@kwdef struct FEMWorkspaceCore{T <: REALSCALAR}
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
	"Thermal conductivity of earth (layers × freq)."
	kappa_g::Matrix{T}
	"Operating temperature [°C]."
	temp::T
	"Number of frequency samples."
	n_frequencies::Int
	"Number of phases in the system."
	n_phases::Int
	"Number of cables in the system."
	n_cables::Int
	"The physical cable system to analyze."
	system::LineCableSystem{T}
	"Earth properties model."
	earth_props::EarthModel{T}
end

"""
$(TYPEDEF)

FEMWorkspaceAmpacity - The workspace for ampacity simulations.

$(TYPEDFIELDS)
"""
@kwdef  struct FEMWorkspaceAmpacity{T <: AbstractFloat}
	"Common parameters definitions."
	core::FEMWorkspaceCore
    "Vector of energizing currents [A]. Index corresponds to phase number."
    energizations::Vector{Complex{T}}
	"Velocity of the ambient wind [m/s]."
	wind_velocity::T


end

"""
$(TYPEDEF)

FEMWorkspaceLineParameters - The workspace for line parameters simulations.

$(TYPEDFIELDS)
"""
@kwdef struct FEMWorkspaceLineParameters{T <: AbstractFloat}
	"Common parameters definitions."
	core::FEMWorkspaceCore
	"Full component-based Z matrix (before bundling/reduction)."
	Z::Array{Complex{T}, 3}
	"Full component-based Y matrix (before bundling/reduction)."
	Y::Array{Complex{T}, 3}
end

function Base.show(io::IO, ::MIME"text/plain", ws::FEMWorkspaceLineParameters)
	print(io, "FEMWorkspaceLineParameters:\n")
	# iterar em todos os campos de ws.core e imprimir
	print(io, "  Z: ", ws.Z, "\n")
	print(io, "  Y: ", ws.Y, "\n")
end

"""
    FEMWorkspace

A type alias (`Union`) that can refer to an instance of either a
[`FEMWorkspaceAmpacity`](@ref) or a [`FEMWorkspaceLineParameters`](@ref).
"""
const FEMWorkspace{T} = Union{FEMWorkspaceAmpacity{T}, FEMWorkspaceLineParameters{T}}
"""
$(TYPEDSIGNATURES)

Internal function to construct the common [`FEMWorkspaceCore`](@ref).
This populates all common fields from a given problem.
"""
function build_workspace(
    problem::Union{LineParametersProblem{U}, AmpacityProblem{U}},
    formulation::FEMFormulation,
) where {U <: REALSCALAR}

    opts = formulation.options
    system = problem.system
    n_frequencies = length(problem.frequencies)
    n_phases = sum(length(cable.design_data.components) for cable in system.cables)
    
    T = BASE_FLOAT
    
    # Pre-allocate 1D arrays
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
    tan_ins = Vector{T}(undef, n_phases)
    phase_map = Vector{Int}(undef, n_phases)
    cable_map = Vector{Int}(undef, n_phases)

    # Fill arrays
    freq .= problem.frequencies

    idx = 0
    for (cable_idx, cable) in enumerate(system.cables)
        for (comp_idx, component) in enumerate(cable.design_data.components)
            idx += 1
            # Geometric properties
            horz[idx] = T(cable.horz)
            vert[idx] = T(cable.vert)
            r_in[idx] = T(component.conductor_group.radius_in)
            r_ext[idx] = T(component.conductor_group.radius_ext)
            r_ins_in[idx] = T(component.insulator_group.radius_in)
            r_ins_ext[idx] = T(component.insulator_group.radius_ext)

            # Material properties
            rho_cond[idx] = T(component.conductor_props.rho)
            alpha_cond[idx] = T(component.conductor_props.alpha)
            mu_cond[idx] = T(component.conductor_props.mu_r)
            eps_cond[idx] = T(component.conductor_props.eps_r)
            rho_ins[idx] = T(component.insulator_props.rho)
            mu_ins[idx] = T(component.insulator_props.mu_r)
            eps_ins[idx] = T(component.insulator_props.eps_r)

            # Calculate loss factor
            ω = 2 * π * f₀  # Using default frequency
            C_eq = T(component.insulator_group.shunt_capacitance)
            G_eq = T(component.insulator_group.shunt_conductance)
            tan_ins[idx] = G_eq / (ω * C_eq)

            # Mapping
            phase_map[idx] = cable.conn[comp_idx]
            cable_map[idx] = cable_idx
        end
    end

    (rho_g, eps_g, mu_g, kappa_g) = _get_earth_data(
        nothing,
        problem.earth_props,
        freq,
        T,
    )

    temp = T(problem.temperature)

    # Create and return the core workspace
    core = FEMWorkspaceCore{T}(
        formulation = formulation,
        opts = opts,
        paths = setup_paths(problem.system, formulation),
        conductors = Vector{GmshObject{<:AbstractEntityData}}(),
        insulators = Vector{GmshObject{<:AbstractEntityData}}(),
        space_regions = Vector{GmshObject{<:AbstractEntityData}}(),
        boundaries = Vector{GmshObject{<:AbstractEntityData}}(),
        unassigned_entities = Dict{Vector{Float64}, AbstractEntityData}(),
        material_registry = Dict{String, Int}(),
        physical_groups = Dict{Int, Material}(),
        freq = freq,
        horz = horz,
        vert = vert,
        r_in = r_in,
        r_ext = r_ext,
        r_ins_in = r_ins_in,
        r_ins_ext = r_ins_ext,
        rho_cond = rho_cond,
        alpha_cond = alpha_cond,
        mu_cond = mu_cond,
        eps_cond = eps_cond,
        rho_ins = rho_ins,
        mu_ins = mu_ins,
        eps_ins = eps_ins,
        tan_ins = tan_ins,
        phase_map = phase_map,
        cable_map = cable_map,
        rho_g = rho_g,
        eps_g = eps_g,
        mu_g = mu_g,
        kappa_g = kappa_g,
        temp = temp,
        n_frequencies = n_frequencies,
        n_phases = n_phases,
        n_cables = system.num_cables,
        system = system,
        earth_props = problem.earth_props,
    )
	return FEMWorkspaceSpecific(problem, core)
	
end

function FEMWorkspaceSpecific(
    problem::LineParametersProblem{T},
    core::FEMWorkspaceCore{T},
) where {T}

    # 1. Add LineParameters-specific fields
    n_phases = core.n_phases
    n_frequencies = core.n_frequencies
    Z = zeros(Complex{T}, n_phases, n_phases, n_frequencies)
    Y = zeros(Complex{T}, n_phases, n_phases, n_frequencies)
    
    # 2. Return the complete workspace
    return FEMWorkspaceLineParameters{T}(
        core = core,
        Z = Z,
        Y = Y
    )
end

function FEMWorkspaceSpecific(
    problem::AmpacityProblem{T},
    core::FEMWorkspaceCore{T},
) where {T}

    # 1. Return the complete workspace
    return FEMWorkspaceAmpacity{T}(
        core = core,
        wind_velocity = problem.wind_velocity,
        energizations = problem.energizations,
    )
end

function init_workspace(
    problem::Union{LineParametersProblem, AmpacityProblem},
    formulation::FEMFormulation, 
    workspace::Union{FEMWorkspace, Nothing}
)
    if isnothing(workspace)
        @debug "Creating new workspace"
        workspace = build_workspace(problem, formulation)
    else
        @debug "Reusing existing workspace"
    end

    opts = formulation.options
    
    # Access paths via the core
    results_dir = workspace.core.paths[:results_dir]
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
