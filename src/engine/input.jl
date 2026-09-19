"""
$(TYPEDEF)

Own the numerical input and reusable storage for one coaxial line-parameter
calculation.

The constructor adapts a completed physical system once, validates the aligned
numerical representation, constructs cable and reduction indices, and
allocates every matrix used by the frequency loop. Each `compute` call owns one
workspace; no mutable state is shared between calculations. Constant fields fix
its input and buffer bindings; reference identity avoids copying this large
record when dispatching heterogeneous equation groups.

$(TYPEDFIELDS)
"""
mutable struct LineParametersWorkspace{
    T <: Real,
    N <: NamedTuple,
    P <: NamedTuple,
    B <: NamedTuple,
    C
}
    "Immutable numerical input derived from the problem and formulation."
    const input::N
    "Physical values and index maps invariant across the frequency loop."
    const invariants::P
    "Mutable numerical storage allocated once for the calculation."
    const buffers::B
    "Optional retained diagnostic arrays, or `nothing`."
    const capture::C

    function LineParametersWorkspace{T, N, P, B, C}(
            input::N,
            invariants::P,
            buffers::B,
            capture::C
    ) where {T <: Real, N <: NamedTuple, P <: NamedTuple, B <: NamedTuple, C}
        return validate(new{T, N, P, B, C}(
            input,
            invariants,
            buffers,
            capture
        ))
    end
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
        P = Array{Complex{T}, 3}(undef, n, n, nf),
        integrals = NamedTuple[]
    )
end

function validate(workspace::LineParametersWorkspace)
    input = workspace.input
    cable = input.cable
    n = input.n_phases
    input.n_frequencies == length(input.freq) || throw(DimensionMismatch(
        "frequency count differs from the frequency vector"
    ))
    input.n_cables == maximum(input.cable_map) || throw(DimensionMismatch(
        "cable count differs from the cable map"
    ))
    for values in (
        input.horz, input.vert, input.phase_map, input.cable_map,
        input.design_map, cable.terminals, cable.positions, cable.r_in,
        cable.r_ext, cable.r_ins_in, cable.r_ins_ext, cable.conductor_materials,
        cable.mu_cond, cable.mu_ins,
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
    sort(vcat(
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
    all(!isempty, workspace.invariants.cable_indices) || throw(ArgumentError(
        "every cable must contain one retained primitive conductor"
    ))
    size(workspace.buffers.Zprimitive) == (n, n) || throw(DimensionMismatch(
        "primitive impedance storage must be $n×$n"
    ))
    size(workspace.buffers.Pprimitive) == (n, n) || throw(DimensionMismatch(
        "primitive potential-coefficient storage must be $n×$n"
    ))
    return workspace
end

"""
$(TYPEDSIGNATURES)

Construct the formulation-independent coaxial input for one validated
line-parameter problem.

The selected designs have already been flattened into frequency-independent
blueprints. This step constructs local cable arrays, physical geometry,
terminal indices, and frequency coordinates once. It does not apply
temperature correction, earth-property/EquivalentHomogeneous formulas, reduction policy, or
allocate formula-specific buffers.

# Arguments

- `problem`: Completed line-parameter problem.
- `blueprints`: One frequency-independent blueprint per selected design.

# Returns

- A read-only named tuple shared by independent formulation workspaces.
"""
function lineinput(
        problem::LineParametersProblem{T},
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
    return (
        freq,
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
end

function lineinput(::Type{T}, input::NamedTuple) where {T <: Real}
    T === eltype(input.freq) && return input
    return merge(input, (
        freq = T.(input.freq), jω = Complex{T}.(input.jω),
        horz = T.(input.horz), vert = T.(input.vert), horz_sep = T.(input.horz_sep),
        cable = convert(LocalCableData{T}, input.cable),
        earth = convert(EarthModel{T}, input.earth),
        temperature = convert(T, input.temperature),
        line_length = convert(T, input.line_length)))
end

function LineParametersWorkspace(
        problem::LineParametersProblem{T},
        formulation::LineParametersFormulation,
        execution::ComputationOptions,
        blueprints::Vector{CableBlueprint{T}}
) where {T <: Real}
    return LineParametersWorkspace(
        problem,
        formulation,
        execution,
        lineinput(problem, blueprints)
    )
end

function LineParametersWorkspace(
        problem::LineParametersProblem{T},
        formulation::LineParametersFormulation,
        execution::ComputationOptions,
        input::NamedTuple
) where {T <: Real}
    cable = input.cable
    horz = input.horz
    horz_sep = input.horz_sep
    vert = input.vert
    phase_map = input.phase_map
    cable_indices = [collect(indices) for indices in cable.assemblies]
    cable_representatives = first.(cable_indices)
    physical_pairs = earth_pairs(
        cable_representatives,
        horz,
        vert,
        horz_sep,
        problem.earth_props
    )
    homogeneous_pairs = _homogeneous_pairs(physical_pairs)
    bindings = map(formulation.methods[(:earth_impedance, :earth_admittance)]) do selected
        leaves = [Formulation(selected, Val.(layer_index(pair))...) for pair in physical_pairs]
        cases = NamedTuple[]
        for leaf in unique(leaves)
            # Validate the selected material inventory while binding it, not
            # by reconstructing selections in a later workspace preflight.
            validate(leaf, leaf.equivalent_earth === nothing ? problem.earth_props : 2)
            indices = findall(value -> value === leaf, leaves)
            push!(cases, earth_bindings(leaf, physical_pairs, homogeneous_pairs, indices))
        end
        (selection = selected, cases = cases)
    end
    for (zi, z) in pairs(bindings.earth_impedance.cases)
        for (pi, p) in pairs(bindings.earth_admittance.cases)
            p.partner == 0 || continue
            paired = earth_bindings(z.selection, p.selection, z, p)
            paired === nothing && continue
            bindings.earth_impedance.cases[zi] = merge(paired.impedance, (partner = pi,))
            bindings.earth_admittance.cases[pi] = merge(paired.admittance, (partner = zi,))
            break
        end
    end
    permutation, reordered_map, kron_map = _reduction_map(phase_map, formulation)
    bundle_pairs = bundle_operations(reordered_map)
    keep_indices = kron_map === nothing ? Int[] : findall(!=(0), kron_map)
    eliminate_indices = kron_map === nothing ? Int[] : findall(==(0), kron_map)
    Invariants = NamedTuple{
        (:cable_indices, :permutation, :reordered_map, :bundle_pairs, :kron_map,
            :keep_indices, :eliminate_indices),
        Tuple{
            Vector{Vector{Int}},
            Vector{Int}, Vector{Int}, Vector{Tuple{Int, Int}},
            Union{Nothing, Vector{Int}}, Vector{Int}, Vector{Int}
        }
    }
    invariants = Invariants((
        cable_indices,
        permutation,
        reordered_map,
        bundle_pairs,
        kron_map,
        keep_indices,
        eliminate_indices
    ))

    scalar = foldl(values(bindings); init = T) do current, bound
        bound.selection isa NamedTuple ||
            return computation_type(current, bound.selection, input.freq)
        foldl(bound.cases; init = current) do representation, call
            computation_type(representation, call.selection, input.freq)
        end
    end
    return LineParametersWorkspace{scalar}(problem, formulation, execution,
        lineinput(scalar, input), merge(invariants, (earth_bindings = bindings,)))
end

function LineParametersWorkspace{T}(
        problem::LineParametersProblem,
        formulation::LineParametersFormulation,
        execution::ComputationOptions,
        input::NamedTuple,
        invariants::NamedTuple
) where {T <: Real}
    cable = input.cable
    n_phases, n_cables, n_frequencies = input.n_phases, input.n_cables, input.n_frequencies
    n_layers = length(cable.dielectric_materials)
    cable_indices = invariants.cable_indices
    representatives = first.(cable_indices)
    eliminate_indices = invariants.eliminate_indices
    nkeep = invariants.kron_map === nothing ? n_phases : length(invariants.keep_indices)
    bindings = invariants.earth_bindings
    geometry = (horizontal = input.horz[representatives], height = input.vert[representatives],
        radius = _outer_radii(input.cable_map, cable.r_ext, cable.r_ins_ext))
    invariants = merge(invariants, (; geometry))
    rho_cond = Vector{T}(undef, length(cable.conductor_materials))
    earth = _earth_data(input, bindings)

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
    Zearth = Matrix{Complex{T}}(undef, n_cables, n_cables)
    Pearth = similar(Zearth)
    MaterialStorage = NamedTuple{(:rho, :epsilon, :mu, :thickness),
        Tuple{Matrix{T}, Matrix{T}, Matrix{T}, Union{Nothing, Vector{T}}}}
    earth_materials = map(_ -> MaterialStorage[], bindings)
    for (family, bound) in pairs(bindings)
        values = getproperty(earth_materials, family)
        for binding in bound.cases
            if family === :earth_admittance && binding.partner != 0
                push!(values, earth_materials.earth_impedance[binding.partner])
                continue
            end
            stratified = media(binding.selection) === Val(:stratified)
            count = stratified ? length(problem.earth_props.layers) : 2
            columns = length(binding.interactions)
            push!(values,
                MaterialStorage((Matrix{T}(undef, count, columns),
                    Matrix{T}(undef, count, columns), Matrix{T}(undef, count, columns),
                    stratified ? Vector{T}(undef, count) : nothing)))
        end
    end
    capture = _capture_buffers(T, input, execution.data.trace)
    observations = capture === nothing ? nothing : capture.integrals
    R = typeof(float(nominal(one(T))))
    # Concrete formulas provision numerical storage, never option-key inspection.
    quadrature = integration_workspace(R, Complex{T}; size = 0)
    largest_cable = maximum(length, cable_indices)
    coefficients = Vector{Complex{T}}(undef, largest_cable)
    tails = similar(coefficients)
    layer_coefficients = Vector{Complex{T}}(undef, n_layers)
    dielectric_admittivity = similar(layer_coefficients)
    buffers = (;
        rho_cond,
        earth,
        dielectric_admittivity,
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
        Zearth,
        Pearth,
        earth_materials,
        quadrature,
        observations,
        layer_coefficients,
        coefficients,
        tails
    )
    # Allocation consumes the same active selections as indexed execution.
    # Unused declarations impose no storage requirement.
    external = map(bindings) do bound
        bound.selection isa NamedTuple ? Tuple(call.selection for call in bound.cases) :
        bound.selection
    end
    selected_internal = formulation.methods.internal_impedance
    internal = if selected_internal isa NamedTuple
        kinds = any(>(0), cable.r_in) ? (:inner, :outer, :transfer) : (:outer,)
        Tuple(unique(Formulation(selected_internal, Val(kind)) for kind in kinds))
    else
        selected_internal
    end
    allocations = merge(formulation.methods, external, (internal_impedance = internal,))
    buffers = initialize_buffers(allocations, T, input, invariants, buffers)
    workspace = LineParametersWorkspace{
        T,
        typeof(input),
        typeof(invariants),
        typeof(buffers),
        typeof(capture)
    }(input, invariants, buffers, capture)
    return workspace
end

function earth_bindings(
        selected::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        physical::AbstractVector{<:EarthPair}, homogeneous, indices)
    pairs = selected.equivalent_earth === nothing ? physical[indices] : homogeneous[indices]
    declarations = validate(selected, pairs)
    reductions = if selected.equivalent_earth === nothing
        nothing
    else
        rule = EquivalentHomogeneous.rule(selected.equivalent_earth)
        foreach(declaration -> validate(declaration.equation, rule), declarations)
        validate(rule, physical[indices])
    end
    interactions = [(index = position, pair = pairs[position],
                        physical_pair = physical[index])
                    for (position, index) in enumerate(indices)]
    equations = [(declaration = declaration,
                     indices = findall(==(declaration), declarations))
                 for declaration in unique(declarations)]
    return (selection = selected, equations, interactions, reductions, partner = 0)
end

earth_bindings(::EarthImpedanceFormulation, ::EarthAdmittanceFormulation, z, p) = nothing
function initialize_buffers(
        ::Union{AbstractFormulation, Nothing}, ::Type, input, invariants, buffers)
    buffers
end

function initialize_buffers(
        selected::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        ::Type{T}, input, invariants, buffers) where {T}
    return initialize_buffers(selected.equivalent_earth, T, input, invariants, buffers)
end

function initialize_buffers(selected::EquivalentHomogeneous.AbstractSequence,
        ::Type{T}, input, invariants, buffers) where {T}
    return initialize_buffers(
        EquivalentHomogeneous.rule(selected), T, input, invariants, buffers)
end

function initialize_buffers(
        selections::Union{NamedTuple, Tuple}, ::Type{T}, input, invariants, buffers) where {T}
    return foldl(values(selections); init = buffers) do accumulated, selected
        initialized = initialize_buffers(selected, T, input, invariants, accumulated)
        # Extension methods may append arrays but cannot replace another owner's
        # storage. Only initially empty QuadGK capacity may be provisioned.
        provision = haskey(accumulated, :quadrature) &&
                    isempty(accumulated.quadrature.segments)
        retained = map(keys(accumulated), values(accumulated),
            values(initialized[keys(accumulated)])) do name, before, after
            (name === :quadrature && provision) || before === after
        end
        all(retained) || throw(ArgumentError(
            "buffer initialization replaced existing storage :$(keys(accumulated)[findfirst(!, retained)])"))
        initialized
    end
end

function _earth_data(input::NamedTuple, bindings::NamedTuple)
    static = (rho = collect(getproperty.(input.earth.layers, :rho)),
        eps_r = collect(getproperty.(input.earth.layers, :eps_r)),
        mu_r = collect(getproperty.(input.earth.layers, :mu_r)))
    needed = any(
        bound -> any(
            case -> !(case.selection.equivalent_earth isa EquivalentHomogeneous.BeforeFD), bound.cases),
        bindings)
    evaluated = needed ?
                map(
        _ -> Matrix{eltype(input.freq)}(undef,
            length(input.earth.layers), input.n_frequencies), static) : nothing
    Evaluated = NamedTuple{(:rho, :eps_r, :mu_r), NTuple{3, Matrix{eltype(input.freq)}}}
    State = NamedTuple{
        (:static, :evaluated), Tuple{typeof(static), Union{Nothing, Evaluated}}}
    return State((static, evaluated))
end

layer_index(problem::LineParametersProblem, horizontal, vertical) =
    layer_index(problem.earth_props, horizontal, vertical)
layer_index(pair::EarthPair) = pair.layers

function layer_index(model::EarthModel, horizontal, vertical)
    vertical > zero(vertical) && return 1
    iszero(vertical) && throw(ArgumentError(
        "a conductor on the air-earth interface has no physical layer"
    ))
    model.vertical_layers && length(model.layers) > 2 &&
        throw(ArgumentError(
            "physical source/target indexing for vertical earth interfaces is not implemented"))
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

"""
$(TYPEDSIGNATURES)

Construct every ordered external interaction from resolved conductor geometry.
Source columns and target rows retain physical earth-layer indices; air is 1.

`cables` contains representative conductor indices. `horizontal`, `vertical`,
and `separation` are aligned coordinates/distances in meters. Diagonal
separations supply the external self radius. Layer assignment uses `earth`'s
physical interfaces and rejects conductors on the air/earth interface.

Return ordered geometry payloads for indexed formula dispatch. No formula
selection or equivalent-earth reduction is performed here.
"""
function earth_pairs(
        cables::AbstractVector{Int},
        horizontal,
        vertical,
        separation,
        earth::EarthModel
)
    T = eltype(vertical)
    pairs = EarthPair{T}[]
    sizehint!(pairs, length(cables)^2)
    placed_layers = [layer_index(earth, horizontal[index], vertical[index]) for index in cables]
    @inbounds for column in eachindex(cables), row in eachindex(cables)

        source = cables[column]
        target = cables[row]
        layers = (placed_layers[column], placed_layers[row])
        push!(pairs,
            EarthPair(
                row,
                column,
                (vertical[source], vertical[target]),
                row == column ? zero(T) : separation[target, source],
                layers; radius = row == column ? separation[target, source] : nothing
            ))
    end
    return pairs
end

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
                pair.layers[1] == 1 ? 1 : 2,
                pair.layers[2] == 1 ? 1 : 2
            ); radius = pair.radius
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
