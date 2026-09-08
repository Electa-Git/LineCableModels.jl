# LineParameters computation remains independent from CableConstants.
_basis_result(parameters, ::LineParametersWorkspace, ::Val{:pul}) = parameters

function _basis_result(
        parameters::LineParameters{T, U, D},
        workspace::LineParametersWorkspace,
        ::Val{:total}
) where {T, U, D}
    line_length = workspace.input.line_length
    impedance = parameters.Z.values .* line_length
    admittance = parameters.Y.values .* line_length
    return LineParameters(
        parameters.domain,
        SeriesImpedance{eltype(impedance), :total}(impedance),
        ShuntAdmittance{eltype(admittance), :total}(admittance),
        parameters.f
    )
end

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
    kron_map = if formulation.options.reduce_bundle
        if formulation.options.kron_reduction
            reduced
        else
            map(eachindex(reduced)) do index
                reordered[index] == 0 ? -1 : reduced[index]
            end
        end
    else
        formulation.options.kron_reduction ? reordered : nothing
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

    @info "Starting line parameters computation"
    for frequency in 1:input.n_frequencies
        homogenize!(workspace, frequency, formulation)
        impedance!(
            Zprimitive,
            workspace,
            frequency,
            formulation
        )
        admittance!(
            Pprimitive,
            workspace,
            frequency,
            formulation
        )
        _reorder_into!(Zbuffer, Zprimitive, permutation)
        _reorder_into!(Pbuffer, Pprimitive, permutation)
        if formulation.options.reduce_bundle
            merge_bundles!(Zbuffer, bundle_pairs)
            merge_bundles!(Pbuffer, bundle_pairs)
        end

        if kron_map === nothing
            formulation.options.ideal_transposition && ideal_transposition!(Zbuffer)
            @views Zout[:, :, frequency] .= Zbuffer

            factorization = lu!(Pbuffer)
            ldiv!(Pinverse, factorization, buffers.identity_full)
            Pinverse .*= input.jω[frequency]
            formulation.options.ideal_transposition && ideal_transposition!(Pinverse)
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
            formulation.options.ideal_transposition && ideal_transposition!(reduced)
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
            formulation.options.ideal_transposition &&
                ideal_transposition!(reduced_inverse)
            @views Yout[:, :, frequency] .= reduced_inverse
        end
    end

    return LineParameters(
        PhaseDomain,
        SeriesImpedance{Complex{T}, :pul}(Zout),
        ShuntAdmittance{Complex{T}, :pul}(Yout),
        input.freq
    )
end

function _retained_details(::LineParametersWorkspace{<:Real, <:NamedTuple, <:NamedTuple,
        <:NamedTuple, Nothing})
    (;)
end

function _retained_details(workspace::LineParametersWorkspace)
    capture = workspace.capture
    capture === nothing && return (;)
    input = workspace.input
    return (
        trace = (
        phase_map = input.phase_map,
        cable_map = input.cable_map,
        Zin = capture.Zin,
        Pin = capture.Pin,
        Zg = capture.Zg,
        Pg = capture.Pg,
        Z = capture.Z,
        P = capture.P
    ),
    )
end

function _finish(
        parameters::LineParameters,
        workspace::LineParametersWorkspace,
        formulation::LineParametersFormulation,
        execution::NamedTuple
)
    parameters = _basis_result(parameters, workspace, execution.output_basis)
    retained = _retained_details(workspace)
    @info "Line parameters computation completed successfully"
    return LineParameters(
        parameters.domain,
        parameters.Z,
        parameters.Y,
        parameters.f,
        retained
    )
end

function _compute(
        engine::LineCableModelsCoaxial,
        problem::LineParametersProblem,
        formulation::LineParametersFormulation,
        execution::NamedTuple,
        input::NamedTuple
)
    requested = formulation
    formulation = Formulation(engine, problem, requested)
    workspace = LineParametersWorkspace(
        problem,
        formulation,
        execution,
        input
    )
    validate(workspace, formulation)
    parameters = _solve!(workspace, formulation)
    result = _finish(parameters, workspace, formulation, execution)
    fields = keys(formulation.methods)
    # Optional FrequencyDependent/EquivalentHomogeneous selections keep a fixed metadata type across grids.
    Identifiers = NamedTuple{fields, NTuple{length(fields), Union{Nothing, Symbol}}}
    requested_ids = map(requested.definitions[fields]) do value
        value === nothing && return nothing
        value isa Symbol && return value
        value isa EquivalentHomogeneous.AbstractSequence &&
            return formula_id(EquivalentHomogeneous.rule(value))
        formula_id(value)
    end
    effective_ids = map(formulation.methods) do value
        value === nothing && return nothing
        value isa EquivalentHomogeneous.AbstractSequence &&
            return formula_id(EquivalentHomogeneous.rule(value))
        formula_id(value)
    end
    selections = (
        requested = Identifiers(requested_ids), effective = Identifiers(effective_ids))
    modified_fields = filter(!=(:pipe_impedance), fields)
    Modifications = NamedTuple{modified_fields, NTuple{length(modified_fields), Bool}}
    modifications::Modifications = Modifications(map(modified_fields) do name
        selected = getproperty(formulation.methods, name)
        selected === nothing && return false
        !isempty(selected.hooks) || !isempty(selected.parameters)
    end)
    Record = NamedTuple{
        (:kind, :source, :target, :options), Tuple{Symbol, Int, Int, NamedTuple}}
    Numerical = NamedTuple{
        (:earth_impedance, :earth_admittance), Tuple{Vector{Record}, Vector{Record}}}
    external_numerical::Numerical = Numerical(map(workspace.invariants.earth_bindings) do bound
        Record[(case.declaration.kind,
                   first(case.interactions).pair.layers...,
                   case.declaration.options) for case in bound.cases]
    end)
    local_fields = (:internal_impedance, :insulation_impedance,
        :insulation_admittance, :semicon_admittance, :earth_properties)
    LocalNumerical = NamedTuple{local_fields, NTuple{length(local_fields), NamedTuple}}
    local_numerical::LocalNumerical = LocalNumerical(map(local_fields) do name
        selected = getproperty(formulation.methods, name)
        selected === nothing && return (;)
        if name === :internal_impedance
            kinds = any(indices -> length(indices) > 1, workspace.invariants.cable_indices) ?
                    (:inner, :outer, :mutual) : (:outer,)
            return selected.options[kinds]
        end
        selected.options
    end)
    FormulaNumerical = NamedTuple{(local_fields..., :earth_impedance, :earth_admittance),
        Tuple{NamedTuple, NamedTuple, NamedTuple, NamedTuple, NamedTuple,
            Vector{Record}, Vector{Record}}}
    numerical::FormulaNumerical = FormulaNumerical((values(local_numerical)...,
        values(external_numerical)...))
    Equivalent = NamedTuple{(:identifier, :order, :parameters, :numerical, :modified),
        Tuple{Symbol, Symbol, NamedTuple, Vector{Record}, Bool}}
    Equivalents = NamedTuple{(:earth_impedance, :earth_admittance),
        Tuple{Union{Nothing, Equivalent}, Union{Nothing, Equivalent}}}
    equivalents::Equivalents = Equivalents(
        map(workspace.invariants.earth_bindings) do bound
        sequence = bound.selection.equivalent_earth
        sequence === nothing && return nothing
        rule = EquivalentHomogeneous.rule(sequence)
        records = Record[(case.equation.arguments[1] === Val(:self) ? :self : :mutual,
                             pair.layers..., case.options)
                         for (case, pair) in
                             zip(bound.reductions, workspace.invariants.earth_pairs)]
        Equivalent((formula_id(rule),
            sequence isa EquivalentHomogeneous.BeforeFD ? :before : :after,
            rule.parameters, unique(records), !isempty(rule.hooks) ||
                !isempty(rule.parameters)))
    end)
    Provenance = NamedTuple{
        (:requested, :effective, :modified, :numerical, :equivalent_earth),
        Tuple{Identifiers, Identifiers, Modifications, FormulaNumerical, Equivalents}}
    provenance::Provenance = Provenance((
        Identifiers(requested_ids), Identifiers(effective_ids),
        modifications, numerical, equivalents))
    return LineParameters(result.domain, result.Z, result.Y, result.f,
        merge(details(result), (; formulations = provenance)))
end

"""
$(TYPEDSIGNATURES)

Bind the selected source equations to the coaxial backend before numerical work.

Earth defaults retain their `:default` identity and prepare the physical
assumptions for the overhead or underground expressions. Existing
explicit choices and FrequencyDependent/EquivalentHomogeneous ordering are retained. Mixed placement requires
an explicitly supported formula. Resolution takes place before workspace
initialization and frequency evaluation.
"""
function Formulation(
        ::LineCableModelsCoaxial,
        problem::LineParametersProblem,
        requested::LineParametersFormulation
)
    methods = merge(requested.methods,
        (
            earth_impedance = Formulation(LineCableModelsCoaxial(), requested.methods.earth_impedance),
            earth_admittance = Formulation(LineCableModelsCoaxial(), requested.methods.earth_admittance)
        ))
    for selected in (methods.earth_impedance, methods.earth_admittance)
        if media(selected) === Val(:homogeneous) && selected.equivalent_earth !== nothing
            validate(problem.earth_props)
            validate(selected, 2)
        else
            validate(selected, problem.earth_props)
        end
        problem.Γ !== nothing && haskey(selected.hooks, :Γ) &&
            throw(ArgumentError(
                "an explicit problem Γ conflicts with the explicit :$(formula_id(selected)) Γ hook"))
        if problem.Γ !== nothing && selected.assumptions.longitudinal === :zero
            all(iszero, problem.Γ) ||
                throw(ArgumentError("formula :$(formula_id(selected)) fixes Γ to zero"))
        end
    end
    return LineParametersFormulation(methods, requested.options, requested.definitions)
end

function _compute(
        engine::LineCableModelsCoaxial,
        problem::LineParametersProblem,
        formulation::LineParametersFormulation,
        execution::NamedTuple
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
        execution::NamedTuple
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
    blueprints = CableBlueprint{T}[flatten(engine, design, T)
                                   for design in problem.system.designs]
    input = lineinput(problem, blueprints)
    first_result = _compute(
        engine,
        problem,
        first(formulations),
        execution,
        input
    )
    values = Vector{typeof(first_result)}(undef, length(formulations))
    values[1] = first_result
    execution.on_result === nothing || execution.on_result(problem, 1, first_result)
    for index in 2:length(formulations)
        value = _compute(
            engine,
            problem,
            formulations[index],
            execution,
            input
        )
        typeof(value) === eltype(values) || throw(ArgumentError(
            "line-parameter formulations produced inconsistent result types",
        ))
        values[index] = value
        execution.on_result === nothing || execution.on_result(problem, index, value)
    end
    return values
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
        options::NamedTuple = (;)
)
    return compute(LineCableModelsCoaxial(), problem, Formulation(); options)
end

"""
$(TYPEDSIGNATURES)

Compute frequency-dependent line parameters with the coaxial backend.

Every nonconcentric cable part must already have an equivalent concentric
representation in the completed data model. The physical system is normalized
once into a backend-owned workspace, and all reusable numerical storage is
allocated before the frequency loop. `trace=true` retains completed
intermediate matrices under `details(result).trace`; it does not change the
result type.

# Arguments

- `problem`: Completed line-parameter problem.
- `formulation`: Selected line-parameter physical methods.

# Keywords

- `options`: Named tuple containing `verbosity`, `output_basis`, `trace`, and
  `on_result`. The optional callable `on_result(problem, index, result)` runs
  synchronously after each completed formulation, including reused results, and
  before the next calculation. `index` is local to the formulation collection
  (`1` for a scalar call). Its return value is ignored; exceptions propagate.
  The callback must not mutate the problem or result. The default is `nothing`.

# Returns

- One [`LineParameters`](@ref) result.
"""
function compute(
        problem::LineParametersProblem,
        formulation::LineParametersFormulation;
        options::NamedTuple = (;)
)
    return compute(LineCableModelsCoaxial(), problem, formulation; options)
end

function compute(
        problem::LineParametersProblem,
        formulations::AbstractVector{<:LineParametersFormulation};
        options::NamedTuple = (;)
)
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
        options::NamedTuple = (;)
)
    execution = computation_options(LineCableModelsCoaxial, options)
    console = ConsoleLogger(stderr, Logging.Debug)
    logger = ConsoleVerbosityLogger(console, execution.verbosity)
    return with_logger(logger) do
        _compute(engine, problem, formulation, execution)
    end
end

function compute(
        engine::LineCableModelsCoaxial,
        problem::LineParametersProblem,
        formulations::AbstractVector{<:LineParametersFormulation};
        options::NamedTuple = (;)
)
    execution = computation_options(LineCableModelsCoaxial, options)
    console = ConsoleLogger(stderr, Logging.Debug)
    logger = ConsoleVerbosityLogger(console, execution.verbosity)
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

function _earth_data(
        formulation::LineParametersFormulation,
        input::NamedTuple
)
    methods = formulation.methods
    consumers = (methods.earth_impedance, methods.earth_admittance)
    static = (rho = collect(getproperty.(input.earth.layers, :rho)),
        eps_r = collect(getproperty.(input.earth.layers, :eps_r)),
        mu_r = collect(getproperty.(input.earth.layers, :mu_r)))
    needed = any(method -> !(method.equivalent_earth isa EquivalentHomogeneous.BeforeFD), consumers)
    evaluated = needed ?
                _earth_data(nothing, methods.earth_properties, input.earth, input.freq) :
                nothing
    return (; static, evaluated)
end

function _earth_data(
        sequence::Union{Nothing, EquivalentHomogeneous.AfterFD},
        relation,
        model::EarthModel{T},
        frequencies::AbstractVector{T}
) where {T <: Real}
    nlayer = length(model.layers)
    nfrequency = length(frequencies)
    rho = Matrix{T}(undef, nlayer, nfrequency)
    eps_r = Matrix{T}(undef, nlayer, nfrequency)
    mu_r = Matrix{T}(undef, nlayer, nfrequency)
    @inbounds for row in eachindex(model.layers)
        static = EarthMaterial(model.layers[row])
        for column in eachindex(frequencies)
            material = row == firstindex(model.layers) ?
                       static : constitutive(relation, static, frequencies[column])
            rho[row, column] = material.rho
            eps_r[row, column] = material.eps_r
            mu_r[row, column] = material.mu_r
        end
    end
    return (; rho, eps_r, mu_r)
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
    input = workspace.input
    for name in (:earth_impedance, :earth_admittance)
        selected = getproperty(formulation.methods, name)
        bindings = getproperty(workspace.invariants.earth_bindings, name)
        destination = getproperty(workspace.buffers.earth_materials, name)
        if media(selected) === Val(:stratified)
            layers!(destination, workspace.invariants.earth.evaluated, input.earth, frequency)
        else
            data = selected.equivalent_earth isa EquivalentHomogeneous.BeforeFD ?
                   workspace.invariants.earth.static : workspace.invariants.earth.evaluated
            homogenize!(destination, selected.equivalent_earth,
                formulation.methods.earth_properties, data, input.earth,
                workspace.invariants.earth_pairs, input.freq[frequency], frequency,
                bindings.reductions)
        end
    end
    return workspace.buffers.earth_materials
end

# Reuse layerwise FrequencyDependent values already evaluated when preparing the calculation.
function layers!(destination, evaluated::NamedTuple, model::EarthModel, frequency::Integer)
    unit=one(eltype(destination.rho))
    epsilon0=unit*88541878128*(unit*10)^(-22)
    mu0=unit*4*(unit*π)*(unit*10)^(-7)
    for row in eachindex(model.layers)
        destination.thickness[row]=model.layers[row].thickness
        for column in axes(destination.rho, 2)
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
        pairs,
        frequency,
        frequency_index::Int,
        bindings
)
    rho = @view evaluated.rho[:, frequency_index]
    eps_r = @view evaluated.eps_r[:, frequency_index]
    mu_r = @view evaluated.mu_r[:, frequency_index]
    air = EarthMaterial(rho[1], eps_r[1], mu_r[1])
    @inbounds for index in eachindex(pairs)
        pair = pairs[index]
        earth = sequence.rule(
            rho, eps_r, mu_r, model, pair, frequency; binding = bindings[index]
        )
        _media!(destination, index, air, earth)
    end
    return destination
end

function homogenize!(
        destination,
        sequence::EquivalentHomogeneous.BeforeFD,
        relation,
        static,
        model::EarthModel,
        pairs,
        frequency,
        frequency_index::Int,
        bindings
)
    air = EarthMaterial(static.rho[1], static.eps_r[1], static.mu_r[1])
    @inbounds for index in eachindex(pairs)
        pair = pairs[index]
        reconstructed = sequence.rule(
            static.rho, static.eps_r, static.mu_r,
            model, pair, frequency; binding = bindings[index]
        )
        earth = constitutive(relation, reconstructed, frequency)
        _media!(destination, index, air, earth)
    end
    return destination
end

function homogenize!(
        destination,
        ::Nothing,
        relation,
        evaluated,
        model::EarthModel,
        pairs,
        frequency,
        frequency_index::Int,
        bindings
)
    rho = @view evaluated.rho[:, frequency_index]
    eps_r = @view evaluated.eps_r[:, frequency_index]
    mu_r = @view evaluated.mu_r[:, frequency_index]
    air = EarthMaterial(rho[1], eps_r[1], mu_r[1])
    earth = EarthMaterial(rho[2], eps_r[2], mu_r[2])
    @inbounds for index in eachindex(pairs)
        _media!(destination, index, air, earth)
    end
    return destination
end
