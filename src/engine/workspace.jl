"""
$(TYPEDEF)

A container for the flattened, type-stable data arrays derived from a
[`LineParametersProblem`](@ref). This struct serves as the primary data source
for all subsequent computational steps.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct CoaxialWorkspace{T<:REALSCALAR,U<:Vector{<:REALSCALAR}}
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
    "Effective earth parameters as a vector of NamedTuples."
    earth::Vector{NamedTuple{(:rho_g, :eps_g, :mu_g),Tuple{U,U,U}}}
    "Operating temperature [°C]."
    temp::T
    "Number of frequency samples."
    n_frequencies::Int
    "Number of phases in the system."
    n_phases::Int
    "Number of cables in the system."
    n_cables::Int
end

"""
$(TYPEDSIGNATURES)

Initializes and populates the [`CoaxialWorkspace`](@ref) by normalizing a
[`LineParametersProblem`](@ref) into flat, type-stable arrays.
"""
function init_workspace(problem::LineParametersProblem, formulation::CoaxialFormulation)

    opts = formulation.options
    setup_logging!(opts.verbosity, opts.logfile)

    system = problem.system
    n_frequencies = length(problem.frequencies)
    n_components = sum(length(cable.design_data.components) for cable in system.cables)

    # Determine the common numeric type for all calculations
    T = _find_common_type(problem)

    # Pre-allocate 1D arrays
    freq = Vector{T}(undef, n_frequencies)
    horz = Vector{T}(undef, n_components)
    vert = Vector{T}(undef, n_components)
    r_in = Vector{T}(undef, n_components)
    r_ext = Vector{T}(undef, n_components)
    r_ins_in = Vector{T}(undef, n_components)
    r_ins_ext = Vector{T}(undef, n_components)
    rho_cond = Vector{T}(undef, n_components)
    mu_cond = Vector{T}(undef, n_components)
    eps_cond = Vector{T}(undef, n_components)
    rho_ins = Vector{T}(undef, n_components)
    mu_ins = Vector{T}(undef, n_components)
    eps_ins = Vector{T}(undef, n_components)
    tan_ins = Vector{T}(undef, n_components)   # Loss tangent for insulator
    phase_map = Vector{Int}(undef, n_components)
    cable_map = Vector{Int}(undef, n_components)

    # Fill arrays, ensuring type promotion
    freq .= T.(problem.frequencies)

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
            mu_cond[idx] = T(component.conductor_props.mu_r)
            eps_cond[idx] = T(component.conductor_props.eps_r)
            rho_ins[idx] = T(component.insulator_props.rho)
            mu_ins[idx] = T(component.insulator_props.mu_r)
            eps_ins[idx] = T(component.insulator_props.eps_r)

            # Calculate loss factor from resistivity
            ω = 2 * π * f₀  # Using default frequency
            C_eq = T(component.insulator_group.shunt_capacitance)
            G_eq = T(component.insulator_group.shunt_conductance)
            tan_ins[idx] = G_eq / (ω * C_eq)

            # Mapping
            phase_map[idx] = cable.conn[comp_idx]
            cable_map[idx] = cable_idx
        end
    end

    earth = _get_earth_data(
        formulation.equivalent_earth,
        problem.earth_props,
        problem.frequencies,
        T
    )

    temp = T(problem.temperature)

    # Construct and return the CoaxialWorkspace struct
    return CoaxialWorkspace(
        freq=freq,
        horz=horz, vert=vert,
        r_in=r_in, r_ext=r_ext,
        r_ins_in=r_ins_in, r_ins_ext=r_ins_ext,
        rho_cond=rho_cond, mu_cond=mu_cond, eps_cond=eps_cond,
        rho_ins=rho_ins, mu_ins=mu_ins, eps_ins=eps_ins, tan_ins=tan_ins,
        phase_map=phase_map, cable_map=cable_map, earth=earth,
        temp=temp, n_frequencies=n_frequencies, n_phases=n_components, n_cables=system.num_cables,
    )
end