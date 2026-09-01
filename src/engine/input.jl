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

Base.eltype(::LineParametersWorkspace{T}) where {T} = T
Base.eltype(::Type{<:LineParametersWorkspace{T}}) where {T} = T

@inline _capture_buffers(::Type, ::Any, ::Val{false}) = nothing

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
    cable = input.cable
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
        input.horz, input.vert, input.phase_map, input.cable_map,
        input.design_map, cable.terminals, cable.positions, cable.r_in,
        cable.r_ext, cable.r_ins_in, cable.r_ins_ext, cable.rho0_cond,
        cable.T0_cond, cable.alpha_cond, cable.mu_cond, cable.mu_ins,
        cable.dielectric_ranges
    )
        length(values) == n || throw(DimensionMismatch(
            "engine input arrays must have $n component entries"
        ))
    end
    n_layers = length(cable.dielectric_materials)
    for values in (
        cable.r_layer_in, cable.r_layer_ext,
        workspace.buffers.layer_coefficients
    )
        length(values) == n_layers || throw(DimensionMismatch(
            "dielectric-layer arrays must contain $n_layers entries"
        ))
    end
    sort!(vcat(
        cable.insulation_indices,
        cable.semicon_indices
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

- `problem`: Completed line-parameter problem.
- `formulation`: Selected physical-method bundle.
- `execution`: Normalized coaxial computation options.
- `blueprints`: One frequency-independent blueprint per selected design.

# Returns

- A validated [`LineParametersWorkspace`](@ref).
"""
function LineParametersWorkspace(
        problem::LineParametersProblem{T},
        formulation::LineParametersFormulation,
        execution::NamedTuple,
        blueprints::Vector{CableBlueprint{T}}
) where {T <: Real}
    system = problem.system
    length(blueprints) == length(system.designs) || throw(DimensionMismatch(
        "line-parameter blueprints must align with the selected system designs",
    ))
    cable = LocalCableData(blueprints)
    n_frequencies = length(problem.frequencies)
    n_phases = length(system.terminal_order)
    length(cable.terminals) == n_phases || throw(DimensionMismatch(
        "DataModel terminal order differs from the cable blueprint count"
    ))
    n_layers = length(cable.dielectric_materials)
    n_cables = length(cable.assemblies)

    freq = copy(problem.frequencies)
    Γ = problem.Γ === nothing ? nothing : copy(problem.Γ)
    jω = Complex{T}.(im .* (2 * (one(first(freq)) * π) .* freq))
    horz = Vector{T}(undef, n_phases)
    horz_sep = Matrix{T}(undef, n_phases, n_phases)
    vert = Vector{T}(undef, n_phases)
    phase_map = copy(system.connection_order)
    design_map = Int[entry.cable for entry in system.terminal_order]
    cable_map = Vector{Int}(undef, n_phases)
    @inbounds for (assembly, indices) in pairs(cable.assemblies), index in indices
        cable_map[index] = assembly
    end

    @inbounds for index in eachindex(cable.terminals)
        canonical = system.terminal_order[index]
        canonical.terminal === cable.terminals[index] || throw(DimensionMismatch(
            "DataModel terminal order is not aligned with the cable blueprint"
        ))
        design_index = design_map[index]
        cable.assembly_designs[cable_map[index]] == design_index ||
            throw(DimensionMismatch(
                "blueprint assembly ownership differs from system terminal order"
            ))
        position = system.positions[design_index]
        local_x, local_y = cable.positions[index]
        horz[index] = position.x + cos(position.φ) * local_x -
                      sin(position.φ) * local_y
        vert[index] = position.y + sin(position.φ) * local_x +
                      cos(position.φ) * local_y
    end
    horizontal_separation!(
        horz_sep,
        horz,
        cable.r_ext,
        cable.r_ins_ext,
        cable_map
    )
    input = (
        freq,
        Γ,
        jω,
        horz,
        horz_sep,
        vert,
        cable,
        phase_map,
        cable_map,
        design_map,
        earth = problem.earth_props,
        temperature = problem.temperature,
        line_length = system.line_length,
        n_frequencies,
        n_phases,
        n_cables
    )

    rho_cond = copy(cable.rho0_cond)
    if formulation.options.temperature_correction
        @inbounds for index in eachindex(rho_cond)
            rho_cond[index] *= one(T) +
                               cable.alpha_cond[index] *
                               (problem.temperature - cable.T0_cond[index])
        end
    end
    earth = _earth_data(formulation, input)
    cable_indices = [collect(indices) for indices in cable.assemblies]
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
