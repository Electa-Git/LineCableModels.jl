"""
$(TYPEDEF)

Own the numerical input and reusable storage for one coaxial line-parameter
calculation.

The constructor adapts a completed physical system once and validates the aligned numerical representation. It constructs cable and reduction indices, and
allocates every matrix used by the frequency loop. Each `compute` call uses an independent workspace. Constant fields fix
its input and buffer bindings. Reference identity avoids copying this large
record when dispatching heterogeneous equation groups.

Bound earth calculations and their material arrays are stored as tuples.
The complete scan specializes on their concrete types once, before frequency
traversal. Conductor layout is stored as runtime data.

$(TYPEDFIELDS)
"""
mutable struct LineParametersWorkspace{
    T <: Real,
    N <: NamedTuple,
    P <: NamedTuple,
    B <: NamedTuple,
    C
}
    # Immutable numerical input derived from the problem and formulation.
    const input::N
    # Physical values and index maps invariant across the frequency loop.
    const plan::P
    # Mutable numerical storage allocated once for the calculation.
    const buffers::B
    # Optional retained diagnostic arrays, or `nothing`.
    const trace::C

    function LineParametersWorkspace{T, N, P, B, C}(
            input::N,
            plan::P,
            buffers::B,
            trace::C
    ) where {T <: Real, N <: NamedTuple, P <: NamedTuple, B <: NamedTuple, C}
        return validate(new{T, N, P, B, C}(
            input,
            plan,
            buffers,
            trace
        ))
    end
end

Base.eltype(::LineParametersWorkspace{T}) where {T} = T
Base.eltype(::Type{<:LineParametersWorkspace{T}}) where {T} = T

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
        cable.mu_r_cond, cable.mu_r_ins,
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
    length(workspace.plan.cable_indices) == input.n_cables || throw(
        DimensionMismatch("cable indices must align with the cable count")
    )
    all(!isempty, workspace.plan.cable_indices) || throw(ArgumentError(
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
terminal indices, and frequency coordinates once.

# Arguments

- `problem`: completed line-parameter problem.
- `blueprints`: one frequency-independent blueprint per selected design.

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
        terminal = system.terminal_order[index]
        terminal.terminal === cable.terminals[index] || throw(DimensionMismatch(
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
    return merge(input,
        (
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
        # The earth is decided once for the slot. A recipe then picks each pair's formula on
        # the decided layers.
        equivalent = _equivalent_earth(selected, problem.earth_props, first(physical_pairs))
        decided = equivalent === nothing ? physical_pairs : homogeneous_pairs
        leaves = selected isa NamedTuple ?
                 [FormulaMethod(selected, pair).selection for pair in decided] :
                 fill(selected, length(decided))
        cases = NamedTuple[]
        for leaf in unique(leaves)
            indices = findall(value -> value === leaf, leaves)
            push!(cases, earth_bindings(leaf, problem.earth_props, equivalent, physical_pairs,
                homogeneous_pairs, indices))
        end
        (selection = selected, cases = cases)
    end
    options = formulation.options.data
    reduction = ReductionPlan(phase_map; options.reduce_bundle, options.kron_reduction,
        options.ideal_transposition)
    plan = (; cable_indices, reduction)

    scalar = foldl(values(bindings); init = T) do current, bound
        bound.selection isa NamedTuple ||
            return computation_type(current, bound.selection, input.freq)
        foldl(bound.cases; init = current) do representation, call
            computation_type(representation, call.selection, input.freq)
        end
    end
    geometry = (layers = [layer_index(pair)[1]
                          for pair in physical_pairs
                          if pair.row == pair.column],)
    return LineParametersWorkspace{scalar}(problem, formulation, execution,
        lineinput(scalar, input), merge(plan, (; geometry)), bindings)
end

function LineParametersWorkspace{T}(
        problem::LineParametersProblem,
        formulation::LineParametersFormulation,
        execution::ComputationOptions,
        input::NamedTuple,
        plan::NamedTuple,
        bindings::NamedTuple
) where {T <: Real}
    cable = input.cable
    n_phases, n_cables, n_frequencies = input.n_phases, input.n_cables, input.n_frequencies
    n_layers = length(cable.dielectric_materials)
    cable_indices = plan.cable_indices
    nkeep = length(plan.reduction.keep)
    geometry = (
        radius = _outer_radii(input.cable_map, cable.r_ext, cable.r_ins_ext),
        layers = plan.geometry.layers)
    # Resolve shared physical outputs once. These temporary lists are not used
    # by the frequency loop. Each completed calculation has its concrete type.
    calculations = NamedTuple[]
    remaining_potential = copy(bindings.earth_admittance.cases)
    for impedance in bindings.earth_impedance.cases
        potential_indices = Int[]
        for (index, potential) in pairs(remaining_potential)
            shared = earth_bindings(impedance.selection, potential.selection,
                impedance, potential)
            shared === nothing && continue
            impedance = shared.impedance
            potential_indices = shared.admittance.output_indices
            deleteat!(remaining_potential, index)
            break
        end
        push!(calculations, merge(impedance[filter(!=(:output_indices), keys(impedance))],
            (impedance_indices = impedance.output_indices, potential_indices)))
    end
    for potential in remaining_potential
        push!(calculations, merge(potential[filter(!=(:output_indices), keys(potential))],
            (impedance_indices = Int[], potential_indices = potential.output_indices)))
    end
    earth_calculations = map(Tuple(calculations)) do calculation
        earth_binding = earth_bindings(calculation.selection, calculation, geometry)
        previous = zeros(Int, length(earth_binding.interactions))
        for group in earth_binding.equations
            for (ordinal, position) in pairs(group.indices)
                for earlier in (ordinal - 1):-1:1
                    previous_interaction = group.indices[earlier]
                    if same_physical_state(earth_binding.reuse_inputs[position],
                        earth_binding.reuse_inputs[previous_interaction])
                        previous[position] = previous_interaction
                        break
                    end
                end
            end
        end
        layers = calculation.earth isa EarthModel ? geometry.layers :
                 [layer == 1 ? 1 : 2 for layer in geometry.layers]
        merge(earth_binding, (equations = Tuple(earth_binding.equations), previous, layers))
    end
    # The workspace layout is independent of the required layer signatures.
    # _solve! specializes once on these concrete tuples before its frequency loop.
    plan = merge(plan,
        NamedTuple{(:earth_calculations, :geometry), Tuple{Tuple, typeof(geometry)}}(
            (earth_calculations, geometry)))
    rho_cond = Vector{T}(undef, length(cable.conductor_materials))
    earth = _earth_data(input, bindings)

    Zprimitive = Matrix{Complex{T}}(undef, n_phases, n_phases)
    Pprimitive = similar(Zprimitive)
    Zout = Array{Complex{T}, 3}(undef, nkeep, nkeep, n_frequencies)
    Yout = similar(Zout)
    Zearth = Matrix{Complex{T}}(undef, n_cables, n_cables)
    Pearth = similar(Zearth)
    # The layered earth keeps every physical layer, and its thicknesses when it has
    # interior layers. A reduction keeps air and one equivalent medium.
    earth_materials = map(earth_calculations) do calculation
        count = calculation.earth isa EarthModel ? length(input.earth.layers) : 2
        columns = length(calculation.interactions)
        thickness = count > 2 ? T[layer.thickness for layer in input.earth.layers] : nothing
        (rho = Matrix{T}(undef, count, columns), epsilon = Matrix{T}(undef, count, columns),
            mu = Matrix{T}(undef, count, columns), thickness)
    end
    trace = if execution.data.trace isa Val{true}
        (Zin = Array{Complex{T}, 3}(undef, n_phases, n_phases, n_frequencies),
            Pin = Array{Complex{T}, 3}(undef, n_phases, n_phases, n_frequencies),
            Zg = Array{Complex{T}, 3}(undef, n_cables, n_cables, n_frequencies),
            Pg = Array{Complex{T}, 3}(undef, n_cables, n_cables, n_frequencies),
            Z = Array{Complex{T}, 3}(undef, n_phases, n_phases, n_frequencies),
            P = Array{Complex{T}, 3}(undef, n_phases, n_phases, n_frequencies),
            integrals = NamedTuple[])
    else
        nothing
    end
    observations = trace === nothing ? nothing : trace.integrals
    largest_cable = maximum(length, cable_indices)
    coefficients = Vector{Complex{T}}(undef, largest_cable)
    tails = similar(coefficients)
    layer_coefficients = Vector{Complex{T}}(undef, n_layers)
    dielectric_admittivity = similar(layer_coefficients)
    buffers = (;
        rho_cond,
        earth,
        dielectric_admittivity,
        Zprimitive,
        Pprimitive,
        Zout,
        Yout,
        Zearth,
        Pearth,
        observations,
        layer_coefficients,
        coefficients,
        tails
    )
    buffers = initialize_buffers(plan.reduction, Complex{T}, input, plan, buffers)
    buffers = initialize_buffers(earth!, T, input, plan, buffers)
    buffers = merge(buffers,
        NamedTuple{(:earth_materials,), Tuple{Tuple}}((earth_materials,)))
    # Allocation consumes the same active selections as indexed execution.
    # Unused declarations do not impose a storage requirement.
    external = map(bindings) do bound
        bound.selection isa NamedTuple ? Tuple(call.selection for call in bound.cases) :
        bound.selection
    end
    # The internal formula has an expression for each surface impedance that the geometry
    # needs.
    for kind in (any(>(0), cable.r_in) ? (:inner, :outer, :transfer) : (:outer,))
        validate(FormulaMethod(formulation.methods.internal_impedance,
            InternalImpedance.internal_impedance, Val(kind)))
    end
    # Concrete formulas provision numerical storage, never option-key inspection.
    allocations = merge(formulation.methods, external)
    buffers = initialize_buffers(allocations, T, input, plan, buffers)
    # Earth traversal clears this shared warning scratch before and after each use.
    buffers = haskey(buffers, :quadrature) ?
              merge(buffers, (quadrature = merge(buffers.quadrature,
        (warnings = buffers.earth_interactions.warnings,)),)) : buffers
    workspace = LineParametersWorkspace{
        T,
        typeof(input),
        typeof(plan),
        typeof(buffers),
        typeof(trace)
    }(input, plan, buffers, trace)
    return workspace
end

"""
$(TYPEDSIGNATURES)

Bind the interactions `physical[indices]` of an earth `model` to the earth formula
`selected`, on the earth that its slot decided: `reduction`, or the layered `model` when it
is `nothing`. The slot decides once, from the signatures of its formulas' expressions. An
explicit `equivalent_earth` reduction always applies. Without one, a formula that does not
admit any layer from 3 to N consumes the `:default` reduction of a model with N > 2 layers,
and every other formula sees the layered earth. The record's `earth` holds the decision:
the `EarthModel` or the reduction. Each reduced interaction stores the binding of the
reduction's equation.
"""
function earth_bindings(
        selected::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        model::EarthModel, reduction, physical::AbstractVector{<:EarthPair}, homogeneous,
        indices)
    reduction === nothing && validate(model, selected)
    pairs = reduction === nothing ? physical[indices] : homogeneous[indices]
    declarations = bindings(selected, pairs)
    distinct = unique(declarations)
    for declaration in distinct
        validate(declaration.equation, model)
    end
    # A reduced interaction also stores the binding of the reduction's own equation.
    interactions = if reduction === nothing
        [(index = position, pair = pairs[position], physical_pair = physical[index])
         for (position, index) in enumerate(indices)]
    else
        rule = EquivalentHomogeneous.rule(reduction)
        foreach(declaration -> validate(rule, declaration.equation), declarations)
        reductions = bindings(rule, physical[indices])
        [(index = position, pair = pairs[position], physical_pair = physical[index],
             reduction = reductions[position])
         for (position, index) in enumerate(indices)]
    end
    equations = [(declaration = declaration,
                     indices = findall(==(declaration), declarations))
                 for declaration in distinct]
    return (selection = selected, earth = reduction === nothing ? model : reduction,
        equations, interactions, output_indices = collect(eachindex(interactions)))
end

earth_bindings(::EarthImpedanceFormulation, ::EarthAdmittanceFormulation, z, p) = nothing

# Indices remain arithmetic inputs unless the selected equation declares otherwise.
function earth_bindings(::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        binding::NamedTuple, geometry::NamedTuple)
    inputs = [(interaction.pair.row, interaction.pair.column)
              for interaction in binding.interactions]
    return merge(binding, (reuse_inputs = inputs,))
end
function initialize_buffers(
        selected::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        ::Type{T}, input, plan, buffers) where {T}
    return initialize_buffers(selected.equivalent_earth, T, input, plan, buffers)
end

function _earth_data(input::NamedTuple, bindings::NamedTuple)
    static = (rho = collect(getproperty.(input.earth.layers, :rho)),
        eps_r = collect(getproperty.(input.earth.layers, :eps_r)),
        mu_r = collect(getproperty.(input.earth.layers, :mu_r)))
    needed = any(
        bound -> any(case -> !(case.earth isa EquivalentHomogeneous.BeforeFD), bound.cases),
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

function layer_index(problem::LineParametersProblem, horizontal, vertical)
    layer_index(problem.earth_props, horizontal, vertical)
end
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
Source columns and target rows retain physical earth-layer indices. Air is 1.

`cables` contains representative conductor indices. `horizontal`, `vertical`,
and `separation` are aligned coordinates and distances in meters. Diagonal
separations supply the external self radius. Layer assignment uses `earth`'s
physical interfaces and rejects conductors on the interface between air and earth.

Return ordered geometry payloads for indexed formula dispatch. Formula
selection and equivalent-earth reduction occur separately.
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
    placed_layers = [layer_index(earth, horizontal[index], vertical[index])
                     for index in cables]
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
