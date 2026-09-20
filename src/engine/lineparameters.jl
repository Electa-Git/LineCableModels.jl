# LineParameters computation remains independent from CableConstants.

@inline function _stash!(destination, frequency::Int, source::AbstractMatrix)
    destination === nothing && return nothing
    @views copyto!(destination[:, :, frequency], source)
    return nothing
end

_capture_target(::Nothing, ::Symbol) = nothing
_capture_target(capture, name::Symbol) = getproperty(capture, name)

@inline function _reorder_into!(destination, source, permutation)
    @inbounds for column in eachindex(permutation), row in eachindex(permutation)

        destination[row, column] = source[permutation[row], permutation[column]]
    end
    return destination
end

function _reduction_map(phase_map, formulation)
    permutation = reorder_indices(phase_map)
    reordered = phase_map[permutation]
    reduced = copy(reordered)
    seen = Set{Int}()
    @inbounds for (index, phase) in pairs(reordered)
        if phase > 0 && phase in seen
            reduced[index] = 0
        elseif phase > 0
            push!(seen, phase)
        end
    end
    kron_map = if formulation.options.data.reduce_bundle
        if formulation.options.data.kron_reduction
            reduced
        else
            map(eachindex(reduced)) do index
                reordered[index] == 0 ? -1 : reduced[index]
            end
        end
    else
        formulation.options.data.kron_reduction ? reordered : nothing
    end
    return permutation, reordered, kron_map
end

function _solve!(
        workspace::LineParametersWorkspace{T},
        formulation::LineParametersFormulation
) where {T <: Real}
    input = workspace.input
    invariants = workspace.invariants
    buffers = workspace.buffers
    workspace.capture===nothing || empty!(workspace.capture.integrals)
    Zbuffer = buffers.Zbuffer
    Pbuffer = buffers.Pbuffer
    Zprimitive = buffers.Zprimitive
    Pprimitive = buffers.Pprimitive
    Pinverse = buffers.Pinverse
    reduced = buffers.reduced
    reduced_inverse = buffers.reduced_inverse
    kron_factor = buffers.kron_factor
    kron_coupling = buffers.kron_coupling
    kron_rhs = buffers.kron_rhs
    Zout = buffers.Zout
    Yout = buffers.Yout
    permutation = invariants.permutation
    bundle_pairs = invariants.bundle_pairs
    kron_map = invariants.kron_map
    keep_indices = invariants.keep_indices
    eliminate_indices = invariants.eliminate_indices

    materials!(workspace, formulation)
    @info "Starting line parameters computation"
    for frequency in 1:input.n_frequencies
        materials!(workspace, formulation, frequency)
        cable_impedance!(Zprimitive, input.cable, buffers.rho_cond,
            formulation.methods, input.jω[frequency]; workspace)
        cable_potential!(Pprimitive, input.cable, buffers.dielectric_admittivity,
            input.jω[frequency], buffers.layer_coefficients, buffers.coefficients, buffers.tails)
        _stash!(_capture_target(workspace.capture, :Zin), frequency, Zprimitive)
        _stash!(_capture_target(workspace.capture, :Pin), frequency, Pprimitive)
        earth!(workspace, frequency)
        _stash!(_capture_target(workspace.capture, :Zg), frequency, buffers.Zearth)
        _stash!(_capture_target(workspace.capture, :Pg), frequency, buffers.Pearth)
        impedance!(Zprimitive, workspace, frequency)
        admittance!(Pprimitive, workspace, frequency)
        _reorder_into!(Zbuffer, Zprimitive, permutation)
        _reorder_into!(Pbuffer, Pprimitive, permutation)
        if formulation.options.data.reduce_bundle
            merge_bundles!(Zbuffer, bundle_pairs)
            merge_bundles!(Pbuffer, bundle_pairs)
        end

        if kron_map === nothing
            formulation.options.data.ideal_transposition && ideal_transposition!(Zbuffer)
            @views Zout[:, :, frequency] .= Zbuffer

            factorization = lu!(Pbuffer)
            ldiv!(Pinverse, factorization, buffers.identity_full)
            Pinverse .*= input.jω[frequency]
            formulation.options.data.ideal_transposition && ideal_transposition!(Pinverse)
            @views Yout[:, :, frequency] .= Pinverse
        else
            kronify!(
                Zbuffer,
                keep_indices,
                eliminate_indices,
                reduced,
                kron_factor,
                kron_coupling,
                kron_rhs
            )
            formulation.options.data.ideal_transposition && ideal_transposition!(reduced)
            @views Zout[:, :, frequency] .= reduced

            kronify!(
                Pbuffer,
                keep_indices,
                eliminate_indices,
                reduced,
                kron_factor,
                kron_coupling,
                kron_rhs
            )
            factorization = lu!(reduced)
            ldiv!(reduced_inverse, factorization, buffers.identity_reduced)
            reduced_inverse .*= input.jω[frequency]
            formulation.options.data.ideal_transposition &&
                ideal_transposition!(reduced_inverse)
            @views Yout[:, :, frequency] .= reduced_inverse
        end
    end

    return workspace
end

function _retained_details(workspace::LineParametersWorkspace{
        <:Real, <:NamedTuple, <:NamedTuple,
        <:NamedTuple, Nothing})
    # Local model diagnostics vary with geometry and selection, not the result type.
    ComputationDetails(NamedTuple{(:shunt_model,), Tuple{NamedTuple}}((workspace.input.cable.shunt_details,)))
end

function _retained_details(workspace::LineParametersWorkspace)
    capture = workspace.capture
    shunt = NamedTuple{(:shunt_model,), Tuple{NamedTuple}}((workspace.input.cable.shunt_details,))
    capture === nothing && return ComputationDetails(shunt)
    input = workspace.input
    return ComputationDetails(merge(shunt,
        (
            trace = (
            phase_map = copy(input.phase_map),
            cable_map = copy(input.cable_map),
            Zin = copy(capture.Zin),
            Pin = copy(capture.Pin),
            Zg = copy(capture.Zg),
            Pg = copy(capture.Pg),
            Z = copy(capture.Z),
            P = copy(capture.P),
            integrals = copy(capture.integrals)
        ),
        )))
end

function _finish(
        workspace::LineParametersWorkspace,
        problem::LineParametersProblem,
        formulation::LineParametersFormulation,
        ::Val{Basis}
) where {Basis}
    impedance = copy(workspace.buffers.Zout)
    admittance = copy(workspace.buffers.Yout)
    if Basis === :total
        impedance .*= workspace.input.line_length
        admittance .*= workspace.input.line_length
    end
    all(isfinite, impedance) || throw(DomainError(impedance,
        "completed series impedance must contain only finite entries"))
    all(isfinite, admittance) || throw(DomainError(admittance,
        "completed shunt admittance must contain only finite entries"))
    retained = _retained_details(workspace)
    names=["cable:$(terminal.cable):$(terminal.terminal)"
           for terminal in problem.system.terminal_order]
    permutation=workspace.invariants.permutation
    indices=workspace.invariants.kron_map === nothing ? permutation :
            permutation[workspace.invariants.keep_indices]
    coordinates=map(indices) do index
        phase=problem.system.connection_order[index]
        members=findall(==(phase), problem.system.connection_order)
        formulation.options.data.reduce_bundle && phase > 0 && length(members) > 1 ?
        "bundle:[" * join(names[members], ",") * "]" : names[index]
    end
    result = LineParameters(PhaseDomain,
        SeriesImpedance{eltype(impedance), Basis}(impedance),
        ShuntAdmittance{eltype(admittance), Basis}(admittance),
        workspace.input.freq,
        ComputationDetails(merge(retained.data,
            (; formulations = NamedTuple(formulation),
                coordinates))))
    return result
end

function _compute(
        engine::LineCableModelsCoaxial,
        problem::LineParametersProblem,
        formulation::LineParametersFormulation,
        execution::ComputationOptions,
        input::NamedTuple
)
    workspace = LineParametersWorkspace(problem, formulation, execution, input)
    _solve!(workspace, formulation)
    return _finish(workspace, problem, formulation, execution.data.output_basis)
end

function _compute(
        engine::LineCableModelsCoaxial,
        problem::LineParametersProblem,
        formulation::LineParametersFormulation,
        execution::ComputationOptions
)
    values = _compute(
        engine,
        problem,
        typeof(formulation)[formulation],
        execution
    )
    return first(values)
end

function _compute(
        engine::LineCableModelsCoaxial,
        problem::LineParametersProblem,
        formulations::AbstractVector{<:LineParametersFormulation},
        execution::ComputationOptions
)
    isempty(formulations) && throw(ArgumentError(
        "line-parameter formulation collections cannot be empty",
    ))
    validate(problem)
    for design in problem.system.designs, formulation in formulations

        Formulation(engine, formulation.methods.pipe_impedance, design)
    end
    maximum(problem.frequencies) > oftype(first(problem.frequencies), 1e8) &&
        @warn("Frequencies above 100 MHz exceed the quasi-TEM validity range.",
            max_frequency=maximum(problem.frequencies),)
    T = eltype(problem)
    blueprints = flatten(engine, problem.system.designs, T, formulations)
    inputs = [lineinput(problem, first(blueprints))]
    for index in 2:length(blueprints)
        previous = findfirst(other -> other === blueprints[index], blueprints)
        push!(inputs, previous < index ? inputs[previous] :
                      lineinput(problem, blueprints[index]))
    end
    return map(formulations, inputs, eachindex(formulations)) do formulation, input, index
        value = _compute(
            engine,
            problem,
            formulation,
            execution,
            input
        )
        execution.data.on_result === nothing ||
            execution.data.on_result(problem, index, value)
        @info "Line parameters computation completed successfully"
        value
    end
end

"""
$(TYPEDSIGNATURES)

Compute line parameters with the coaxial backend and default formulation.

# Arguments

- `problem`: Completed line-parameter problem.

# Keywords

- `options`: Coaxial-backend computation options.

# Returns

- One [`LineParameters`](@ref) result.
"""
function compute(
        problem::LineParametersProblem;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions()
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    return compute(LineCableModelsCoaxial(), problem, Formulation(); options)
end

"""
$(TYPEDSIGNATURES)

Compute frequency-dependent line parameters with the coaxial backend.

The completed data model supplies the equivalent concentric representation
used for series impedance and ordinary radial dielectric intervals. Eligible
open wire/tape domains retain their physical geometry for the explicitly selected
`shunt_model=:boundary` calculation; the default uses annular geometry. The
physical system is normalized once into a backend-owned
workspace, and all reusable numerical storage is allocated before the frequency
loop. `trace=true` retains completed
intermediate matrices under `details(result).data.trace`; it does not change the
result type.

# Arguments

- `problem`: Completed line-parameter problem.
- `formulation`: Selected line-parameter physical methods.

# Keywords

- `options`: Named tuple containing `verbosity`, `output_basis`, `trace`, and
  `on_result`. The optional callable `on_result(problem, index, result)` runs
  synchronously after each completed formulation and
  before the next calculation. `index` is local to the formulation collection
  (`1` for a scalar call). Its return value is ignored; exceptions propagate.
  The callback must not mutate the problem or result. The default is `nothing`.
  The selected local shunt coefficients are constructed in the cable blueprints;
  no separate preparation call or execution option is required.

# Returns

- One [`LineParameters`](@ref) result.
"""
function compute(
        problem::LineParametersProblem,
        formulation::LineParametersFormulation;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions()
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    return compute(LineCableModelsCoaxial(), problem, formulation; options)
end

function compute(
        problem::LineParametersProblem,
        formulations::AbstractVector{<:LineParametersFormulation};
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions()
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    return compute(LineCableModelsCoaxial(), problem, formulations; options)
end

"""
$(TYPEDSIGNATURES)

Compute line parameters through an explicit coaxial backend tag.

The tag owns execution dispatch while `formulation.methods` retains the
selected physical recipes. Ordinary callers can omit the tag and use the
two-argument `compute` method.
"""
function compute(
        engine::LineCableModelsCoaxial,
        problem::LineParametersProblem,
        formulation::LineParametersFormulation = Formulation();
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions()
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    execution = computation_options(LineCableModelsCoaxial, options)
    console = ConsoleLogger(stderr, Logging.Debug)
    logger = ConsoleVerbosityLogger(console, execution.data.verbosity)
    return with_logger(logger) do
        _compute(engine, problem, formulation, execution)
    end
end

function compute(
        engine::LineCableModelsCoaxial,
        problem::LineParametersProblem,
        formulations::AbstractVector{<:LineParametersFormulation};
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions()
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    execution = computation_options(LineCableModelsCoaxial, options)
    console = ConsoleLogger(stderr, Logging.Debug)
    logger = ConsoleVerbosityLogger(console, execution.data.verbosity)
    return with_logger(logger) do
        _compute(engine, problem, formulations, execution)
    end
end

function computation_details(
        ::Type{<:LineParametersFormulation},
        result::LineParameters
)::ComputationDetails
    return details(result)
end

function computation_details(
        ::Type{<:LineCableModelsFEM},
        result::LineParameters
)::ComputationDetails
    return details(result)
end

function materials!(
        destination::NamedTuple,
        relation,
        model::EarthModel{T},
        frequencies::AbstractVector{T}; workspace = nothing
) where {T <: Real}
    rho, eps_r, mu_r = destination.rho, destination.eps_r, destination.mu_r
    @inbounds for row in eachindex(model.layers)
        static = EarthMaterial(model.layers[row])
        for column in eachindex(frequencies)
            material = row == firstindex(model.layers) ?
                       static :
                       constitutive(relation, static, frequencies[column]; workspace)
            rho[row, column] = material.rho
            eps_r[row, column] = material.eps_r
            mu_r[row, column] = material.mu_r
        end
    end
    return destination
end

function materials!(workspace::LineParametersWorkspace, formulation::LineParametersFormulation)
    input, buffers = workspace.input, workspace.buffers
    for (index, material) in pairs(input.cable.conductor_materials)
        buffers.rho_cond[index] = constitutive(formulation.methods.temperature_dependence,
            material, input.temperature; workspace)
    end
    buffers.earth.evaluated === nothing || materials!(buffers.earth.evaluated,
        formulation.methods.earth_properties, input.earth, input.freq; workspace)
    return workspace
end

function materials!(
        workspace::LineParametersWorkspace, formulation::LineParametersFormulation,
        frequency::Int)
    homogenize!(workspace, frequency, formulation)
    dielectric!(workspace.buffers.dielectric_admittivity, workspace.input.cable,
        formulation.methods, workspace.input.freq[frequency], workspace.input.temperature; workspace)
    return workspace
end

@inline function _media!(destination, column::Int, air, earth)
    unit = one(earth.rho)
    epsilon0 = unit * 88541878128 * (unit * 10)^(-22)
    mu0 = unit * 4 * (unit * π) * (unit * 10)^(-7)
    destination.rho[1, column] = air.rho
    destination.rho[2, column] = earth.rho
    destination.epsilon[1, column] = epsilon0 * air.eps_r
    destination.epsilon[2, column] = epsilon0 * earth.eps_r
    destination.mu[1, column] = mu0 * air.mu_r
    destination.mu[2, column] = mu0 * earth.mu_r
    return nothing
end

function homogenize!(
        workspace::LineParametersWorkspace,
        frequency::Int,
        formulation::LineParametersFormulation
)
    for (family, cases) in
        pairs(map(bound -> bound.cases, workspace.invariants.earth_bindings))
        destinations = getproperty(workspace.buffers.earth_materials, family)
        for (index, binding) in pairs(cases)
            family === :earth_admittance && binding.partner != 0 && continue
            homogenize!(destinations[index], binding, workspace, frequency,
                formulation.methods.earth_properties)
        end
    end
    return workspace.buffers.earth_materials
end

function homogenize!(destination,
        binding::NamedTuple,
        workspace::LineParametersWorkspace, frequency_index::Int, relation)
    state = workspace.buffers.earth
    model = workspace.input.earth
    frequency = workspace.input.freq[frequency_index]
    selected = binding.selection
    if media(selected) === Val(:stratified)
        layers!(destination, state.evaluated, model, frequency_index, binding.interactions)
    else
        data = selected.equivalent_earth isa EquivalentHomogeneous.BeforeFD ?
               state.static : state.evaluated
        homogenize!(destination, selected.equivalent_earth, relation, data, model,
            binding.interactions, frequency, frequency_index, binding.reductions; workspace)
    end
    return destination
end

# Reuse layerwise FrequencyDependent values already evaluated when preparing the calculation.
function layers!(destination, evaluated::NamedTuple,
        model::EarthModel, frequency::Integer, interactions)
    unit=one(eltype(destination.rho))
    epsilon0=unit*88541878128*(unit*10)^(-22)
    mu0=unit*4*(unit*π)*(unit*10)^(-7)
    for row in eachindex(model.layers)
        destination.thickness[row]=model.layers[row].thickness
        for interaction in interactions
            column = interaction.index
            destination.rho[row, column]=evaluated.rho[row, frequency]
            destination.epsilon[row, column]=epsilon0*evaluated.eps_r[row, frequency]
            destination.mu[row, column]=mu0*evaluated.mu_r[row, frequency]
        end
    end
    return destination
end

function homogenize!(
        destination,
        sequence::EquivalentHomogeneous.AfterFD,
        relation,
        evaluated,
        model::EarthModel,
        interactions,
        frequency,
        frequency_index::Int,
        bindings; workspace = nothing
)
    rho = @view evaluated.rho[:, frequency_index]
    eps_r = @view evaluated.eps_r[:, frequency_index]
    mu_r = @view evaluated.mu_r[:, frequency_index]
    air = EarthMaterial(rho[1], eps_r[1], mu_r[1])
    @inbounds for (index, interaction) in enumerate(interactions)
        pair = interaction.physical_pair
        earth = sequence.rule(
            rho, eps_r, mu_r, model, pair, frequency; binding = bindings[index], workspace
        )
        _media!(destination, interaction.index, air, earth)
    end
    return destination
end

function homogenize!(
        destination,
        sequence::EquivalentHomogeneous.BeforeFD,
        relation,
        static,
        model::EarthModel,
        interactions,
        frequency,
        frequency_index::Int,
        bindings; workspace = nothing
)
    air = EarthMaterial(static.rho[1], static.eps_r[1], static.mu_r[1])
    @inbounds for (index, interaction) in enumerate(interactions)
        pair = interaction.physical_pair
        reconstructed = sequence.rule(
            static.rho, static.eps_r, static.mu_r,
            model, pair, frequency; binding = bindings[index], workspace
        )
        earth = constitutive(relation, reconstructed, frequency; workspace)
        _media!(destination, interaction.index, air, earth)
    end
    return destination
end

function homogenize!(
        destination,
        ::Nothing,
        relation,
        evaluated,
        model::EarthModel,
        interactions,
        frequency,
        frequency_index::Int,
        bindings; workspace = nothing
)
    rho = @view evaluated.rho[:, frequency_index]
    eps_r = @view evaluated.eps_r[:, frequency_index]
    mu_r = @view evaluated.mu_r[:, frequency_index]
    air = EarthMaterial(rho[1], eps_r[1], mu_r[1])
    earth = EarthMaterial(rho[2], eps_r[2], mu_r[2])
    @inbounds for interaction in interactions
        _media!(destination, interaction.index, air, earth)
    end
    return destination
end
