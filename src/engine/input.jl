"""
$(TYPEDEF)

Immutable flattened solver input derived from one validated problem.

The conductor resistivities stored here are base-state values. Operating
temperature correction is calculated into a local array by [`compute`](@ref).

$(TYPEDFIELDS)
"""
struct AnalyticalInput{T <: Real, E <: EarthModel{T}}
    freq::Vector{T}
    jω::Vector{Complex{T}}
    horz::Vector{T}
    horz_sep::Matrix{T}
    vert::Vector{T}
    r_in::Vector{T}
    r_ext::Vector{T}
    r_ins_in::Vector{T}
    r_ins_ext::Vector{T}
    rho0_cond::Vector{T}
    T0_cond::Vector{T}
    alpha_cond::Vector{T}
    mu_cond::Vector{T}
    eps_cond::Vector{T}
    rho_ins::Vector{T}
    mu_ins::Vector{T}
    eps_ins::Vector{T}
    tan_ins::Vector{T}
    insulator_layer_ranges::Vector{UnitRange{Int}}
    r_ins_layer_in::Vector{T}
    r_ins_layer_ext::Vector{T}
    rho_ins_layer::Vector{T}
    eps_ins_layer::Vector{T}
    phase_map::Vector{Int}
    cable_map::Vector{Int}
    earth::E
    line_length::T
    n_frequencies::Int
    n_phases::Int
    n_cables::Int
end

Base.eltype(::AnalyticalInput{T}) where {T} = T
Base.eltype(::Type{<:AnalyticalInput{T}}) where {T} = T

function _check_analytical_input(input::AnalyticalInput)
    n = input.n_phases
    input.n_frequencies == length(input.freq) || throw(DimensionMismatch(
        "frequency count differs from the frequency vector",
    ))
    input.n_cables == maximum(input.cable_map) || throw(DimensionMismatch(
        "cable count differs from the cable map",
    ))
    for values in (
        input.horz, input.vert, input.r_in, input.r_ext, input.r_ins_in,
        input.r_ins_ext, input.rho0_cond, input.T0_cond, input.alpha_cond,
        input.mu_cond, input.eps_cond, input.rho_ins, input.mu_ins,
        input.eps_ins, input.tan_ins, input.phase_map, input.cable_map,
        input.insulator_layer_ranges
    )
        length(values) == n || throw(DimensionMismatch(
            "flattened analytical input arrays must have $n component entries",
        ))
    end
    size(input.horz_sep) == (n, n) || throw(DimensionMismatch(
        "horizontal separation matrix must be $n×$n",
    ))
    return nothing
end

function Validation.rules(::Type{<:AnalyticalInput})
    (Validation.OwnerRule(:analytical_input_dimensions, _check_analytical_input),)
end

function AnalyticalInput(problem::LineParametersProblem{T}) where {T <: Real}
    system = problem.system
    n_frequencies = length(problem.frequencies)
    n_phases = sum(
        length(cable.design_data.components) for cable in system.cables
    )
    n_layers = sum(
        length(component.insulator_group.layers)
    for cable in system.cables
    for component in cable.design_data.components
    )
    n_cables = ncables(system)

    freq = copy(problem.frequencies)
    jω = Complex{T}.(im .* (2 * (one(first(freq)) * π) .* freq))
    horz = Vector{T}(undef, n_phases)
    horz_sep = Matrix{T}(undef, n_phases, n_phases)
    vert = Vector{T}(undef, n_phases)
    r_in = Vector{T}(undef, n_phases)
    r_ext = Vector{T}(undef, n_phases)
    r_ins_in = Vector{T}(undef, n_phases)
    r_ins_ext = Vector{T}(undef, n_phases)
    rho0_cond = Vector{T}(undef, n_phases)
    T0_cond = Vector{T}(undef, n_phases)
    alpha_cond = Vector{T}(undef, n_phases)
    mu_cond = Vector{T}(undef, n_phases)
    eps_cond = Vector{T}(undef, n_phases)
    rho_ins = Vector{T}(undef, n_phases)
    mu_ins = Vector{T}(undef, n_phases)
    eps_ins = Vector{T}(undef, n_phases)
    tan_ins = Vector{T}(undef, n_phases)
    insulator_layer_ranges = Vector{UnitRange{Int}}(undef, n_phases)
    r_ins_layer_in = Vector{T}(undef, n_layers)
    r_ins_layer_ext = Vector{T}(undef, n_layers)
    rho_ins_layer = Vector{T}(undef, n_layers)
    eps_ins_layer = Vector{T}(undef, n_layers)
    phase_map = Vector{Int}(undef, n_phases)
    cable_map = Vector{Int}(undef, n_phases)

    component_index = 0
    layer_index = 0
    for (cable_index, cable) in enumerate(system.cables)
        for (local_index, component) in enumerate(cable.design_data.components)
            component_index += 1
            horz[component_index] = cable.horz
            vert[component_index] = cable.vert
            r_in[component_index] = component.conductor_group.r_in
            r_ext[component_index] = component.conductor_group.r_ex
            r_ins_in[component_index] = component.insulator_group.r_in
            r_ins_ext[component_index] = component.insulator_group.r_ex
            rho0_cond[component_index] = component.conductor_props.rho
            T0_cond[component_index] = component.conductor_props.T0
            alpha_cond[component_index] = component.conductor_props.alpha
            mu_cond[component_index] = component.conductor_props.mu_r
            eps_cond[component_index] = component.conductor_props.eps_r
            rho_ins[component_index] = component.insulator_props.rho
            mu_ins[component_index] = component.insulator_props.mu_r
            eps_ins[component_index] = component.insulator_props.eps_r
            reference_ω = 2 * (one(first(freq)) * π) *
                          component.insulator_group.reference_frequency
            tan_ins[component_index] = DataModel.loss_tangent(
                component.insulator_group.shunt_conductance,
                component.insulator_group.shunt_capacitance,
                reference_ω
            )
            first_layer = layer_index + 1
            for layer in component.insulator_group.layers
                layer_index += 1
                r_ins_layer_in[layer_index] = layer.r_in
                r_ins_layer_ext[layer_index] = layer.r_ex
                rho_ins_layer[layer_index] = layer.material_props.rho
                eps_ins_layer[layer_index] = layer.material_props.eps_r
            end
            insulator_layer_ranges[component_index] = first_layer:layer_index
            phase_map[component_index] = cable.conn[local_index]
            cable_map[component_index] = cable_index
        end
    end
    horizontal_separation!(horz_sep, horz, r_ext, r_ins_ext, cable_map)
    return validate(AnalyticalInput(
        freq, jω, horz, horz_sep, vert, r_in, r_ext, r_ins_in, r_ins_ext,
        rho0_cond, T0_cond, alpha_cond, mu_cond, eps_cond, rho_ins, mu_ins,
        eps_ins, tan_ins, insulator_layer_ranges, r_ins_layer_in,
        r_ins_layer_ext, rho_ins_layer, eps_ins_layer, phase_map, cable_map,
        problem.earth_props, system.line_length, n_frequencies, n_phases,
        n_cables
    ))
end
@inline function _outer_radii(cable_map, r_ext, r_ins_ext)
    length(cable_map) == length(r_ext) == length(r_ins_ext) ||
        throw(DimensionMismatch("cable maps and radius vectors must align"))
    outer = fill(zero(eltype(r_ext)), maximum(cable_map))
    @inbounds for index in eachindex(cable_map)
        cable = cable_map[index]
        outer[cable] = max(outer[cable], r_ext[index], r_ins_ext[index])
    end
    return outer
end

function horizontal_separation!(destination, horizontal, r_ext, r_ins_ext, cable_map)
    n = length(horizontal)
    size(destination) == (n, n) || throw(DimensionMismatch(
        "horizontal separation matrix must be $n×$n",
    ))
    outer = _outer_radii(cable_map, r_ext, r_ins_ext)
    @inbounds for column in 1:n, row in 1:n

        destination[row, column] = cable_map[row] == cable_map[column] ?
                                   outer[cable_map[row]] :
                                   abs(horizontal[row] - horizontal[column])
    end
    return destination
end
