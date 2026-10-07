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
    geometry = (radius = _outer_radii(input.cable_map, cable.r_ext, cable.r_ins_ext),
        layers = [layer_index(pair)[1] for pair in physical_pairs if pair.row == pair.column])
    # Each slot decides its earth once. A recipe then picks each pair's formula on the decided
    # layers. Impedance and admittance calculations that solve the same system merge.
    earth = EarthPlan(
        EarthPlan(formulation.methods.earth_impedance, problem.earth_props, physical_pairs,
            geometry),
        EarthPlan(formulation.methods.earth_admittance, problem.earth_props, physical_pairs,
            geometry))
    options = formulation.options.data
    reduction = ReductionPlan(phase_map; options.reduce_bundle, options.kron_reduction,
        options.ideal_transposition)
    plan = (; cable_indices, reduction, earth, geometry)
    # The earth formulas that the computation uses: each slot's formula, or the formulas of a
    # recipe that a calculation uses, the others `nothing`. They widen the scalar type and
    # provision their arrays.
    formulas = map(formulation.methods[(:earth_impedance, :earth_admittance)]) do selected
        selected isa NamedTuple || return selected
        map(selected) do leaf
            any(calculation -> any(entry -> entry !== nothing && entry.formula === leaf,
                (calculation.impedance, calculation.admittance)), earth.calculations) ?
            leaf : nothing
        end
    end
    scalar = computation_type(T, formulas, input.freq)
    return LineParametersWorkspace{scalar}(problem, formulation, execution,
        lineinput(scalar, input), plan, formulas)
end

function LineParametersWorkspace{T}(
        problem::LineParametersProblem,
        formulation::LineParametersFormulation,
        execution::ComputationOptions,
        input::NamedTuple,
        plan::NamedTuple,
        formulas::NamedTuple
) where {T <: Real}
    cable = input.cable
    n_phases, n_cables, n_frequencies = input.n_phases, input.n_cables, input.n_frequencies
    n_layers = length(cable.dielectric_materials)
    cable_indices = plan.cable_indices
    nkeep = length(plan.reduction.keep)
    # The plan keeps the conductor radii in the computation's scalar type.
    plan = merge(plan, (geometry = (radius = convert(Vector{T}, plan.geometry.radius),
        layers = plan.geometry.layers),))
    rho_cond = Vector{T}(undef, length(cable.conductor_materials))

    Zprimitive = Matrix{Complex{T}}(undef, n_phases, n_phases)
    Pprimitive = similar(Zprimitive)
    Zout = Array{Complex{T}, 3}(undef, nkeep, nkeep, n_frequencies)
    Yout = similar(Zout)
    Zearth = Matrix{Complex{T}}(undef, n_cables, n_cables)
    Pearth = similar(Zearth)
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
    # The internal formula has an expression for each surface impedance that the geometry
    # needs.
    for kind in (any(>(0), cable.r_in) ? (:inner, :outer, :transfer) : (:outer,))
        validate(Expression(formulation.methods.internal_impedance,
            InternalImpedance.internal_impedance, Val(kind)))
    end
    # Concrete formulas provision numerical storage, never option-key inspection. A recipe
    # formula provisions arrays only when a calculation uses it. The earth plan then
    # provisions its calculations.
    buffers = initialize_buffers(merge(formulation.methods, formulas), T, input, plan, buffers)
    buffers = initialize_buffers(plan.earth, T, input, plan, buffers)
    workspace = LineParametersWorkspace{
        T,
        typeof(input),
        typeof(plan),
        typeof(buffers),
        typeof(trace)
    }(input, plan, buffers, trace)
    return workspace
end

function initialize_buffers(
        selected::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        ::Type{T}, input, plan, buffers) where {T}
    return initialize_buffers(selected.equivalent_earth, T, input, plan, buffers)
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
