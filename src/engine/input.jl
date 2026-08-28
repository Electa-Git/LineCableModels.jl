"""
$(TYPEDEF)

Own the normalized data and reusable storage for one native line-parameter
calculation.

The constructor adapts a completed physical system once, validates the aligned
numerical representation, precomputes cable and reduction indices, and
allocates every matrix used by the frequency loop. Each `compute` call owns one
workspace; no mutable state is shared between calculations.

$(TYPEDFIELDS)
"""
struct LineParametersWorkspace{
        T <: Real,
        N <: NamedTuple,
        P <: NamedTuple,
        B <: NamedTuple,
        C
}
    "Immutable numerical data normalized from the problem and formulation."
    normalized::N
    "Precomputed values and index maps reused by every frequency."
    prepared::P
    "Mutable numerical storage allocated once for the calculation."
    buffers::B
    "Optional retained diagnostic arrays, or `nothing`."
    capture::C
end

Base.eltype(::LineParametersWorkspace{T}) where {T} = T
Base.eltype(::Type{<:LineParametersWorkspace{T}}) where {T} = T

function _capture_buffers(
        ::Type,
        normalized::NamedTuple,
        ::Val{false}
)
    return nothing
end

function _capture_buffers(
        ::Type{T},
        normalized::NamedTuple,
        ::Val{true}
) where {T <: Real}
    n = normalized.n_phases
    nc = normalized.n_cables
    nf = normalized.n_frequencies
    return (
        Zin = Array{Complex{T}, 3}(undef, n, n, nf),
        Pin = Array{Complex{T}, 3}(undef, n, n, nf),
        Zg = Array{Complex{T}, 3}(undef, nc, nc, nf),
        Pg = Array{Complex{T}, 3}(undef, nc, nc, nf),
        Z = Array{Complex{T}, 3}(undef, n, n, nf),
        P = Array{Complex{T}, 3}(undef, n, n, nf)
    )
end

function _check_line_parameters_workspace(workspace::LineParametersWorkspace)
    input = workspace.normalized
    n = input.n_phases
    input.n_frequencies == length(input.freq) || throw(DimensionMismatch(
        "frequency count differs from the frequency vector"
    ))
    input.n_cables == maximum(input.cable_map) || throw(DimensionMismatch(
        "cable count differs from the cable map"
    ))
    for values in (
        input.horz, input.vert, input.r_in, input.r_ext, input.r_ins_in,
        input.r_ins_ext, input.rho0_cond, input.T0_cond, input.alpha_cond,
        input.mu_cond, input.eps_cond, input.rho_ins, input.mu_ins,
        input.eps_ins, input.tan_ins, input.phase_map, input.cable_map,
        input.insulator_layer_ranges
    )
        length(values) == n || throw(DimensionMismatch(
            "normalized engine arrays must have $n component entries"
        ))
    end
    size(input.horz_sep) == (n, n) || throw(DimensionMismatch(
        "horizontal separation matrix must be $n×$n"
    ))
    length(workspace.prepared.cable_indices) == input.n_cables || throw(
        DimensionMismatch("cable indices must align with the cable count")
    )
    all(!isempty, workspace.prepared.cable_indices) || throw(ArgumentError(
        "every cable must contain one retained primitive conductor"
    ))
    size(workspace.buffers.Zprimitive) == (n, n) || throw(DimensionMismatch(
        "primitive impedance storage must be $n×$n"
    ))
    size(workspace.buffers.Pprimitive) == (n, n) || throw(DimensionMismatch(
        "primitive potential-coefficient storage must be $n×$n"
    ))
    return nothing
end

function Validation.rules(::Type{<:LineParametersWorkspace})
    (Validation.OwnerRule(
        :line_parameters_workspace_dimensions,
        _check_line_parameters_workspace
    ),)
end

"""
$(TYPEDSIGNATURES)

Construct the native workspace for one validated line-parameter problem.

All physical adaptation, temperature correction, earth-property evaluation,
index preparation, and numerical allocation occurs before the frequency loop.

# Arguments

- `engine`: Native backend identity.
- `problem`: Completed line-parameter problem.
- `formulation`: Selected physical-method bundle.
- `execution`: Normalized native computation options.

# Returns

- A validated [`LineParametersWorkspace`](@ref).
"""
function LineParametersWorkspace(
        ::LineCableModelsEngine,
        problem::LineParametersProblem{T},
        formulation::LineParametersFormulation,
        execution::NamedTuple
) where {T <: Real}
    system = problem.system
    dielectric_frequency = oftype(first(problem.frequencies), 50)
    terminals_by_design = [
        homogeneous_components(formulation, design, dielectric_frequency, T)
        for design in system.designs
    ]
    n_frequencies = length(problem.frequencies)
    n_phases = sum(length, terminals_by_design)
    n_layers = sum(
        length(terminal.dielectric.layers)
        for terminals in terminals_by_design
        for terminal in terminals
    )
    n_cables = ncables(system)

    freq = copy(problem.frequencies)
    jω = Complex{T}.(im .* (2 * (one(first(freq)) * π) .* freq))
    horz = Vector{T}(undef, n_phases)
    horz_sep = Matrix{T}(undef, n_phases, n_phases)
    vert = Vector{T}(undef, n_phases)
    r_in_values = Vector{T}(undef, n_phases)
    r_ext_values = Vector{T}(undef, n_phases)
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
    for (cable_index, (terminals, position, connections)) in enumerate(zip(
            terminals_by_design,
            system.positions,
            system.connections
    ))
        for (local_index, terminal) in enumerate(terminals)
            component_index += 1
            conductor = terminal.conductor
            dielectric = terminal.dielectric
            horz[component_index] = position.x
            vert[component_index] = position.y
            r_in_values[component_index] = conductor.r_in
            r_ext_values[component_index] = conductor.r_ex
            r_ins_in[component_index] = dielectric.r_in
            r_ins_ext[component_index] = dielectric.r_ex
            rho0_cond[component_index] = conductor.material.rho
            T0_cond[component_index] = conductor.material.T0
            alpha_cond[component_index] = conductor.material.alpha
            mu_cond[component_index] = conductor.material.mu_r
            eps_cond[component_index] = conductor.material.eps_r
            rho_ins[component_index] = dielectric.material.rho
            mu_ins[component_index] = dielectric.material.mu_r
            eps_ins[component_index] = dielectric.material.eps_r
            reference_ω = 2 * (one(first(freq)) * π) * dielectric.frequency
            tan_ins[component_index] = isempty(dielectric.layers) ? zero(T) :
                                       DataModel.loss_tangent(
                dielectric.shunt_conductance,
                dielectric.shunt_capacitance,
                reference_ω
            )
            first_layer = layer_index + 1
            for layer in dielectric.layers
                layer_index += 1
                r_ins_layer_in[layer_index] = layer.r_in
                r_ins_layer_ext[layer_index] = layer.r_ex
                rho_ins_layer[layer_index] = layer.material.rho
                eps_ins_layer[layer_index] = layer.material.eps_r
            end
            insulator_layer_ranges[component_index] = first_layer:layer_index
            phase_map[component_index] = connections[local_index]
            cable_map[component_index] = cable_index
        end
    end
    horizontal_separation!(
        horz_sep,
        horz,
        r_ext_values,
        r_ins_ext,
        cable_map
    )
    normalized = (
        freq,
        jω,
        horz,
        horz_sep,
        vert,
        r_in = r_in_values,
        r_ext = r_ext_values,
        r_ins_in,
        r_ins_ext,
        rho0_cond,
        T0_cond,
        alpha_cond,
        mu_cond,
        eps_cond,
        rho_ins,
        mu_ins,
        eps_ins,
        tan_ins,
        insulator_layer_ranges,
        r_ins_layer_in,
        r_ins_layer_ext,
        rho_ins_layer,
        eps_ins_layer,
        phase_map,
        cable_map,
        earth = problem.earth_props,
        line_length = system.line_length,
        n_frequencies,
        n_phases,
        n_cables
    )

    rho_cond = copy(rho0_cond)
    if formulation.options.temperature_correction
        @inbounds for index in eachindex(rho_cond)
            rho_cond[index] *= DataModel.temperature_factor(
                alpha_cond[index], problem.temperature, T0_cond[index]
            )
        end
    end
    earth = _earth_data(formulation, normalized)
    cable_indices = [Int[] for _ in 1:n_cables]
    @inbounds for index in 1:n_phases
        push!(cable_indices[cable_map[index]], index)
    end
    cable_representatives = first.(cable_indices)
    permutation, reordered_map, kron_map = _reduction_map(phase_map, formulation)
    nkeep = kron_map === nothing ? n_phases : count(!=(0), kron_map)
    Prepared = NamedTuple{
        (:rho_cond, :earth, :cable_indices, :cable_representatives,
         :permutation, :reordered_map, :kron_map, :nkeep),
        Tuple{
            Vector{T}, typeof(earth), Vector{Vector{Int}}, Vector{Int},
            Vector{Int}, Vector{Int}, Union{Nothing, Vector{Int}}, Int
        }
    }
    prepared = Prepared((
        rho_cond,
        earth,
        cable_indices,
        cable_representatives,
        permutation,
        reordered_map,
        kron_map,
        nkeep
    ))

    Zbuffer = Matrix{Complex{T}}(undef, n_phases, n_phases)
    Pbuffer = similar(Zbuffer)
    Zprimitive = similar(Zbuffer)
    Pprimitive = similar(Zbuffer)
    Pinverse = similar(Zbuffer)
    reduced = Matrix{Complex{T}}(undef, nkeep, nkeep)
    reduced_inverse = similar(reduced)
    identity_full = Matrix{Complex{T}}(I, n_phases, n_phases)
    identity_reduced = Matrix{Complex{T}}(I, nkeep, nkeep)
    Zout = Array{Complex{T}, 3}(undef, nkeep, nkeep, n_frequencies)
    Yout = similar(Zout)
    earth_matrix = Matrix{Complex{T}}(undef, n_cables, n_cables)
    largest_cable = maximum(length, cable_indices)
    coefficients = Vector{Complex{T}}(undef, largest_cable)
    tails = similar(coefficients)
    buffers = (;
        Zbuffer,
        Pbuffer,
        Zprimitive,
        Pprimitive,
        Pinverse,
        reduced,
        reduced_inverse,
        identity_full,
        identity_reduced,
        Zout,
        Yout,
        earth_matrix,
        coefficients,
        tails
    )
    capture = _capture_buffers(T, normalized, execution.trace)
    workspace = LineParametersWorkspace{
        T,
        typeof(normalized),
        typeof(prepared),
        typeof(buffers),
        typeof(capture)
    }(normalized, prepared, buffers, capture)
    validate(workspace)
    return workspace
end

@inline function _outer_radii(cable_map, r_ext_values, r_ins_ext)
    length(cable_map) == length(r_ext_values) == length(r_ins_ext) ||
        throw(DimensionMismatch("cable maps and radius vectors must align"))
    outer = fill(zero(eltype(r_ext_values)), maximum(cable_map))
    @inbounds for index in eachindex(cable_map)
        cable = cable_map[index]
        outer[cable] = max(outer[cable], r_ext_values[index], r_ins_ext[index])
    end
    return outer
end

function horizontal_separation!(
        destination,
        horizontal,
        r_ext_values,
        r_ins_ext,
        cable_map
)
    n = length(horizontal)
    size(destination) == (n, n) || throw(DimensionMismatch(
        "horizontal separation matrix must be $n×$n"
    ))
    outer = _outer_radii(cable_map, r_ext_values, r_ins_ext)
    @inbounds for column in 1:n, row in 1:n
        destination[row, column] = cable_map[row] == cable_map[column] ?
                                   outer[cable_map[row]] :
                                   abs(horizontal[row] - horizontal[column])
    end
    return destination
end
