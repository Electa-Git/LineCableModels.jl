"""
$(TYPEDEF)

Own the numerical input and reusable storage for one coaxial line-parameter
calculation.

The constructor adapts a completed physical system once, validates the aligned
numerical representation, constructs cable and reduction indices, and
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
    "Immutable numerical input derived from the problem and formulation."
    input::N
    "Physical values and index maps invariant across the frequency loop."
    invariants::P
    "Mutable numerical storage allocated once for the calculation."
    buffers::B
    "Optional retained diagnostic arrays, or `nothing`."
    capture::C
end

constitutive(
    ::LineParametersFormulation,
    ::Val{:conductor},
    material::Material
) = material

function constitutive(
        ::LineParametersFormulation,
        ::Val{Kind},
        ::Material
) where {Kind}
    throw(ArgumentError(
        "the coaxial formulation does not support material class :$Kind"
    ))
end

Base.eltype(::LineParametersWorkspace{T}) where {T} = T
Base.eltype(::Type{<:LineParametersWorkspace{T}}) where {T} = T

function _capture_buffers(
        ::Type,
        input::NamedTuple,
        ::Val{false}
)
    return nothing
end

function _capture_buffers(
        ::Type{T},
        input::NamedTuple,
        ::Val{true}
) where {T <: Real}
    n = input.n_phases
    nc = input.n_cables
    nf = input.n_frequencies
    return (
        Zin = Array{Complex{T}, 3}(undef, n, n, nf),
        Pin = Array{Complex{T}, 3}(undef, n, n, nf),
        Zg = Array{Complex{T}, 3}(undef, nc, nc, nf),
        Pg = Array{Complex{T}, 3}(undef, nc, nc, nf),
        Z = Array{Complex{T}, 3}(undef, n, n, nf),
        P = Array{Complex{T}, 3}(undef, n, n, nf)
    )
end

function _check_workspace(workspace::LineParametersWorkspace)
    input = workspace.input
    n = input.n_phases
    input.n_frequencies == length(input.freq) || throw(DimensionMismatch(
        "frequency count differs from the frequency vector"
    ))
    input.Γ !== nothing && length(input.Γ) != input.n_frequencies &&
        throw(DimensionMismatch(
            "longitudinal propagation constants must align with frequencies"
        ))
    input.n_cables == maximum(input.cable_map) || throw(DimensionMismatch(
        "cable count differs from the cable map"
    ))
    for values in (
        input.horz, input.vert, input.r_in, input.r_ext, input.r_ins_in,
        input.r_ins_ext, input.rho0_cond, input.T0_cond, input.alpha_cond,
        input.mu_cond, input.eps_cond, input.mu_ins,
        input.phase_map, input.cable_map,
        input.insulator_layer_ranges
    )
        length(values) == n || throw(DimensionMismatch(
            "engine input arrays must have $n component entries"
        ))
    end
    n_layers = length(input.dielectric_materials)
    for values in (
        input.r_ins_layer_in, input.r_ins_layer_ext, input.rho_ins_layer,
        input.eps_ins_layer, workspace.buffers.layer_coefficients
    )
        length(values) == n_layers || throw(DimensionMismatch(
            "dielectric-layer arrays must contain $n_layers entries"
        ))
    end
    sort!(vcat(
        input.insulation_layer_indices,
        input.semicon_layer_indices
    )) == collect(1:n_layers) || throw(DimensionMismatch(
        "insulation and semicon indices must partition the dielectric layers"
    ))
    size(input.horz_sep) == (n, n) || throw(DimensionMismatch(
        "horizontal separation matrix must be $n×$n"
    ))
    length(workspace.invariants.cable_indices) == input.n_cables || throw(
        DimensionMismatch("cable indices must align with the cable count")
    )
    length(workspace.invariants.earth_pairs) ==
    input.n_cables * (input.n_cables + 1) ÷ 2 || throw(DimensionMismatch(
        "earth pairs must contain the upper triangular cable interactions"
    ))
    length(workspace.invariants.homogeneous_pairs) ==
    length(workspace.invariants.earth_pairs) || throw(DimensionMismatch(
        "physical and homogeneous earth pairs must align"
    ))
    all(!isempty, workspace.invariants.cable_indices) || throw(ArgumentError(
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
    (Validation.OwnerRule(:line_parameters_workspace_dimensions, _check_workspace),)
end

"""
$(TYPEDSIGNATURES)

Construct the coaxial workspace for one validated line-parameter problem.

All physical adaptation, temperature correction, earth-property evaluation,
index construction, and numerical allocation occurs before the frequency loop.

# Arguments

- `engine`: Coaxial backend identity.
- `problem`: Completed line-parameter problem.
- `formulation`: Selected physical-method bundle.
- `execution`: Normalized coaxial computation options.

# Returns

- A validated [`LineParametersWorkspace`](@ref).
"""
function LineParametersWorkspace(
        ::LineCableModelsCoaxial,
        problem::LineParametersProblem{T},
        formulation::LineParametersFormulation,
        execution::NamedTuple
) where {T <: Real}
    system = problem.system
    # DataModel's homogeneous export view requires a positive dielectric
    # frequency. Coaxial shunt calculations ignore its equivalent dielectric
    # material and consume the canonical physical layers below.
    material_reference_frequency = first(problem.frequencies)
    terminals_by_design = [DataModel.flatten(
                               design, material_reference_frequency, T
                           )
                           for design in system.designs]
    for (cable_index, terminals) in enumerate(terminals_by_design)
        reference_position = first(terminals).conductor.position
        all(terminals) do terminal
            DataModel.same_radial_position(
                terminal.conductor.position,
                reference_position
            )
        end || throw(ArgumentError(
            "the coaxial backend requires one concentric " *
            "assembly member per cable; cable $cable_index contains " *
            "independently positioned members"
        ))
    end
    n_frequencies = length(problem.frequencies)
    n_phases = length(system.terminal_order)
    sum(length, terminals_by_design) == n_phases || throw(DimensionMismatch(
        "DataModel terminal order differs from the homogeneous component count"
    ))
    n_layers = sum(
        length(terminal.dielectric.layers)
    for terminals in terminals_by_design
    for terminal in terminals
    )
    n_cables = ncables(system)

    freq = copy(problem.frequencies)
    Γ = problem.Γ === nothing ? nothing : copy(problem.Γ)
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
    mu_ins = Vector{T}(undef, n_phases)
    insulator_layer_ranges = Vector{UnitRange{Int}}(undef, n_phases)
    r_ins_layer_in = Vector{T}(undef, n_layers)
    r_ins_layer_ext = Vector{T}(undef, n_layers)
    rho_ins_layer = Vector{T}(undef, n_layers)
    eps_ins_layer = Vector{T}(undef, n_layers)
    dielectric_materials = Vector{Material{T}}(undef, n_layers)
    insulation_layer_indices = Int[]
    semicon_layer_indices = Int[]
    sizehint!(insulation_layer_indices, n_layers)
    sizehint!(semicon_layer_indices, n_layers)
    phase_map = copy(system.connection_order)
    cable_map = Int[entry.cable for entry in system.terminal_order]

    component_index = 0
    layer_index = 0
    for (cable_index, (terminals, position)) in enumerate(zip(
        terminals_by_design,
        system.positions
    ))
        for terminal in terminals
            component_index += 1
            canonical = system.terminal_order[component_index]
            canonical == (cable = cable_index, terminal = terminal.name) ||
                throw(DimensionMismatch(
                    "DataModel terminal order is not aligned with homogeneous components"
                ))
            conductor = terminal.conductor
            dielectric = terminal.dielectric
            conductor_material = constitutive(
                formulation,
                Val(conductor.material.kind),
                conductor.material
            )
            local_x, local_y = conductor.position
            horz[component_index] = position.x +
                                    cos(position.φ) * local_x -
                                    sin(position.φ) * local_y
            vert[component_index] = position.y +
                                    sin(position.φ) * local_x +
                                    cos(position.φ) * local_y
            r_in_values[component_index] = conductor.r_in
            r_ext_values[component_index] = conductor.r_ex
            r_ins_in[component_index] = dielectric.r_in
            r_ins_ext[component_index] = dielectric.r_ex
            rho0_cond[component_index] = conductor_material.rho
            T0_cond[component_index] = conductor_material.T0
            alpha_cond[component_index] = conductor_material.alpha
            mu_cond[component_index] = conductor_material.mu_r
            eps_cond[component_index] = conductor_material.eps_r
            for layer in dielectric.layers
                layer.material.kind in (:insulator, :semicon) || throw(ArgumentError(
                    "the coaxial backend requires dielectric layers to be :insulator or :semicon; " *
                    "got :$(layer.material.kind)"
                ))
            end
            mu_ins[component_index] = isempty(dielectric.layers) ? one(T) :
                                      DataModel.equivalent_dielectric_permeability(
                dielectric.layers,
                conductor.num_turns,
                conductor.r_ex,
                dielectric.r_ex
            )
            first_layer = layer_index + 1
            for layer in dielectric.layers
                layer_index += 1
                r_ins_layer_in[layer_index] = layer.r_in
                r_ins_layer_ext[layer_index] = layer.r_ex
                rho_ins_layer[layer_index] = layer.material.rho
                eps_ins_layer[layer_index] = layer.material.eps_r
                dielectric_materials[layer_index] = layer.material
                if layer.material.kind === :insulator
                    push!(insulation_layer_indices, layer_index)
                else
                    push!(semicon_layer_indices, layer_index)
                end
            end
            insulator_layer_ranges[component_index] = first_layer:layer_index
        end
    end
    horizontal_separation!(
        horz_sep,
        horz,
        r_ext_values,
        r_ins_ext,
        cable_map
    )
    input = (
        freq,
        Γ,
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
        mu_ins,
        insulator_layer_ranges,
        r_ins_layer_in,
        r_ins_layer_ext,
        rho_ins_layer,
        eps_ins_layer,
        dielectric_materials,
        insulation_layer_indices,
        semicon_layer_indices,
        phase_map,
        cable_map,
        earth = problem.earth_props,
        temperature = problem.temperature,
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
    earth = _earth_data(formulation, input)
    cable_indices = [Int[] for _ in 1:n_cables]
    @inbounds for (index, terminal) in pairs(system.terminal_order)
        push!(cable_indices[terminal.cable], index)
    end
    cable_representatives = first.(cable_indices)
    earth_pairs = _earth_pairs(
        cable_representatives,
        horz,
        vert,
        horz_sep,
        problem.earth_props
    )
    homogeneous_pairs = _homogeneous_pairs(earth_pairs)
    permutation, reordered_map, kron_map = _reduction_map(phase_map, formulation)
    bundle_pairs = bundle_operations(reordered_map)
    keep_indices = kron_map === nothing ? Int[] : findall(!=(0), kron_map)
    eliminate_indices = kron_map === nothing ? Int[] : findall(==(0), kron_map)
    nkeep = kron_map === nothing ? n_phases : count(!=(0), kron_map)
    Invariants = NamedTuple{
        (:rho_cond, :earth, :cable_indices, :cable_representatives, :earth_pairs,
            :homogeneous_pairs,
            :permutation, :reordered_map, :bundle_pairs, :kron_map,
            :keep_indices, :eliminate_indices, :nkeep),
        Tuple{
            Vector{T}, typeof(earth), Vector{Vector{Int}}, Vector{Int},
            Vector{EarthPair{T}}, Vector{EarthPair{T}},
            Vector{Int}, Vector{Int}, Vector{Tuple{Int, Int}},
            Union{Nothing, Vector{Int}}, Vector{Int}, Vector{Int}, Int
        }
    }
    invariants = Invariants((
        rho_cond,
        earth,
        cable_indices,
        cable_representatives,
        earth_pairs,
        homogeneous_pairs,
        permutation,
        reordered_map,
        bundle_pairs,
        kron_map,
        keep_indices,
        eliminate_indices,
        nkeep
    ))

    Zbuffer = Matrix{Complex{T}}(undef, n_phases, n_phases)
    Pbuffer = similar(Zbuffer)
    Zprimitive = similar(Zbuffer)
    Pprimitive = similar(Zbuffer)
    Pinverse = similar(Zbuffer)
    reduced = Matrix{Complex{T}}(undef, nkeep, nkeep)
    reduced_inverse = similar(reduced)
    neliminate = length(eliminate_indices)
    kron_factor = Matrix{Complex{T}}(undef, neliminate, neliminate)
    kron_coupling = Matrix{Complex{T}}(undef, nkeep, neliminate)
    kron_rhs = Matrix{Complex{T}}(undef, neliminate, nkeep)
    identity_full = Matrix{Complex{T}}(I, n_phases, n_phases)
    identity_reduced = Matrix{Complex{T}}(I, nkeep, nkeep)
    Zout = Array{Complex{T}, 3}(undef, nkeep, nkeep, n_frequencies)
    Yout = similar(Zout)
    earth_matrix = Matrix{Complex{T}}(undef, n_cables, n_cables)
    pair_count = length(earth_pairs)
    earth_media = (
        rho = Matrix{T}(undef, 2, pair_count),
        epsilon = Matrix{T}(undef, 2, pair_count),
        mu = Matrix{T}(undef, 2, pair_count)
    )
    n_earth_layers = length(problem.earth_props.layers)
    earth_layers = (
        rho = Matrix{T}(undef, n_earth_layers, pair_count),
        epsilon = Matrix{T}(undef, n_earth_layers, pair_count),
        mu = Matrix{T}(undef, n_earth_layers, pair_count),
        thickness = Vector{T}(undef, n_earth_layers)
    )
    integration_type = typeof(float(nominal(one(T))))
    earth_impedance_segments = alloc_segbuf(
        integration_type,
        Complex{T},
        integration_type;
        size = 128
    )
    earth_admittance_segments = alloc_segbuf(
        integration_type,
        Complex{T},
        integration_type;
        size = 128
    )
    largest_cable = maximum(length, cable_indices)
    coefficients = Vector{Complex{T}}(undef, largest_cable)
    tails = similar(coefficients)
    layer_coefficients = Vector{Complex{T}}(undef, n_layers)
    buffers = (;
        Zbuffer,
        Pbuffer,
        Zprimitive,
        Pprimitive,
        Pinverse,
        reduced,
        reduced_inverse,
        kron_factor,
        kron_coupling,
        kron_rhs,
        identity_full,
        identity_reduced,
        Zout,
        Yout,
        earth_matrix,
        earth_media,
        earth_layers,
        earth_impedance_segments,
        earth_admittance_segments,
        layer_coefficients,
        coefficients,
        tails
    )
    capture = _capture_buffers(T, input, execution.trace)
    workspace = LineParametersWorkspace{
        T,
        typeof(input),
        typeof(invariants),
        typeof(buffers),
        typeof(capture)
    }(input, invariants, buffers, capture)
    validate(workspace)
    return workspace
end

function _earth_layer(model::EarthModel, horizontal, vertical)
    vertical > zero(vertical) && return 1
    iszero(vertical) && throw(ArgumentError(
        "a conductor on the air-earth interface has no physical layer"
    ))
    model.vertical_layers && return 2
    depth = -vertical
    boundary = zero(depth)
    @inbounds for layer in 2:length(model.layers)
        thickness = model.layers[layer].thickness
        isinf(thickness) && return layer
        boundary += thickness
        depth <= boundary && return layer
    end
    throw(ArgumentError(
        "conductor depth $depth m is outside the earth-layer model"
    ))
end

function _earth_pairs(
        cables::AbstractVector{Int},
        horizontal,
        vertical,
        separation,
        earth::EarthModel
)
    T = eltype(vertical)
    pairs = EarthPair{T}[]
    sizehint!(pairs, length(cables) * (length(cables) + 1) ÷ 2)
    @inbounds for column in eachindex(cables), row in firstindex(cables):column

        left = cables[row]
        right = cables[column]
        layers = (
            _earth_layer(earth, horizontal[left], vertical[left]),
            _earth_layer(earth, horizontal[right], vertical[right])
        )
        push!(pairs, EarthPair(
            row,
            column,
            (vertical[left], vertical[right]),
            separation[left, right],
            layers
        ))
    end
    return pairs
end

@inline function _layout(pair::EarthPair)
    source_air = pair.layers[1] == 1
    target_air = pair.layers[2] == 1
    source_air && target_air && return Val(:overhead)
    !source_air && !target_air && return Val(:underground)
    return Val(:mixed)
end

@inline _homogeneous_layer(layer::Int) = layer == 1 ? 1 : 2

function _homogeneous_pairs(pairs::AbstractVector{<:EarthPair{T}}) where {T <: Real}
    mapped = Vector{EarthPair{T}}(undef, length(pairs))
    @inbounds for index in eachindex(pairs)
        pair = pairs[index]
        mapped[index] = EarthPair(
            pair.row,
            pair.column,
            pair.heights,
            pair.separation,
            (
                _homogeneous_layer(pair.layers[1]),
                _homogeneous_layer(pair.layers[2])
            )
        )
    end
    return mapped
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
