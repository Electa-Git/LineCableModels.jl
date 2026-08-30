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
            reciprocity!(Zbuffer)
            formulation.options.ideal_transposition && ideal_transposition!(Zbuffer)
            @views Zout[:, :, frequency] .= Zbuffer

            factorization = lu!(Pbuffer)
            ldiv!(Pinverse, factorization, buffers.identity_full)
            Pinverse .*= input.jω[frequency]
            reciprocity!(Pinverse)
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
            reciprocity!(reduced)
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
            reciprocity!(reduced_inverse)
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
        execution::NamedTuple
)
    validate(problem)
    workspace = LineParametersWorkspace(engine, problem, formulation, execution)
    parameters = _solve!(workspace, formulation)
    return _finish(parameters, workspace, formulation, execution)
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

- `options`: Named tuple containing `verbosity`, `output_basis`, and `trace`.

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
    execution = computation_options(Val(LineCableModelsCoaxial), options)
    console = ConsoleLogger(stderr, Logging.Debug)
    logger = ConsoleVerbosityLogger(console, execution.verbosity)
    return with_logger(logger) do
        _compute(engine, problem, formulation, execution)
    end
end

computation_owner(::LineParametersFormulation) = LineCableModelsCoaxial
computation_owner(::LineCableModelsFEM) = LineCableModelsFEM

function computation_details(
        ::Val{LineCableModelsCoaxial},
        result::LineParameters
)::ComputationDetails
    return details(result)
end

computation_details(
    ::Val{LineCableModelsCoaxial},
    ::DataModel.CableConstants
)::ComputationDetails = (;)

function computation_details(
        ::Val{LineCableModelsFEM},
        result::LineParameters
)::ComputationDetails
    return details(result)
end

"""
$(TYPEDSIGNATURES)

Calculate the per-length cable constants of one retained cable terminal through
the ordinary line-parameter engine.

The selected terminal is assigned to phase 1; every other retained terminal is
assigned to phase 0 and removed by Kron reduction. The default problem uses a
single cable at `(0, -1)` m, a 1 m line, 100 Ω·m earth, 20 °C, and 1 mHz.

# Arguments

- `design`: Completed physical cable design.

# Keywords

- `core`: Retained terminal used as phase 1.
- `formulation`: Line-parameter formulation used by `compute`.
- `options`: Computation options forwarded to `compute`.
- `position`: Cable placement in metres.
- `line_length`: Physical line length \\[m\\].
- `earth_props`: Static earth model.
- `temperature`: Operating temperature \\[°C\\].
- `frequency`: Analysis frequency \\[Hz\\].

# Returns

- A [`CableConstants`](@ref) value observed from the line parameters.
"""
function DataModel.CableConstants(
        design::CableDesign;
        core::Symbol = :core,
        formulation::AbstractFormulation = Formulation(),
        options::NamedTuple = (;),
        position = DataModel.Pose2(0, -1, 0),
        line_length::Real = 1,
        earth_props::EarthModel = EarthModel(100),
        temperature::Real = 20,
        frequency::Real = 1e-3
)
    matches = findall(==(core), design.terminal_order)
    length(matches) == 1 || throw(ArgumentError(
        "cable design must contain exactly one retained terminal named :$core"
    ))
    phases = ntuple(
        index -> index == only(matches) ? 1 : 0,
        length(design.terminal_order)
    )
    connections = NamedTuple{Tuple(design.terminal_order)}(phases)
    problem = LineParametersProblem(
        design,
        position;
        connections,
        line_length,
        temperature,
        earth_props,
        frequencies = [frequency]
    )
    parameters = compute(problem, formulation; options)
    return DataModel.CableConstants(
        @observe(parameters, R[1, 1, 1]),
        @observe(parameters, L[1, 1, 1]),
        @observe(parameters, C[1, 1, 1])
    )
end

function _earth_data(
        formulation::LineParametersFormulation,
        input::NamedTuple
)
    sequence = formulation.methods.equivalent_earth
    if sequence === nothing && length(input.earth.layers) != 2 &&
       (media(formulation.methods.earth_impedance) === Val(:homogeneous) ||
        media(formulation.methods.earth_admittance) === Val(:homogeneous))
        throw(ArgumentError(
            "a homogeneous earth-return formulation requires an EHEM rule for a multilayer earth"
        ))
    end
    return _earth_data(
        sequence,
        formulation.methods.earth_properties,
        input.earth,
        input.freq
    )
end

function _earth_data(
        sequence::Union{Nothing, EHEM.AfterFD},
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

function _earth_data(
        ::EHEM.BeforeFD,
        relation,
        model::EarthModel{T},
        frequencies::AbstractVector{T}
) where {T <: Real}
    nlayer = length(model.layers)
    rho = Vector{T}(undef, nlayer)
    eps_r = Vector{T}(undef, nlayer)
    mu_r = Vector{T}(undef, nlayer)
    @inbounds for row in eachindex(model.layers)
        material = EarthMaterial(model.layers[row])
        rho[row] = material.rho
        eps_r[row] = material.eps_r
        mu_r[row] = material.mu_r
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
    layers!(
        workspace.buffers.earth_layers,
        formulation.methods.earth_properties,
        workspace.input.earth,
        workspace.input.freq[frequency]
    )
    return homogenize!(
        workspace.buffers.earth_media,
        formulation.methods.equivalent_earth,
        formulation.methods.earth_properties,
        workspace.invariants.earth,
        workspace.input.earth,
        workspace.invariants.earth_pairs,
        workspace.input.freq[frequency],
        frequency
    )
end

function layers!(destination, relation, model::EarthModel, frequency)
    unit = one(frequency)
    epsilon0 = unit * 88541878128 * (unit * 10)^(-22)
    mu0 = unit * 4 * (unit * π) * (unit * 10)^(-7)
    @inbounds for row in eachindex(model.layers)
        static = EarthMaterial(model.layers[row])
        material = row == firstindex(model.layers) ?
                   static : constitutive(relation, static, frequency)
        destination.thickness[row] = model.layers[row].thickness
        for column in axes(destination.rho, 2)
            destination.rho[row, column] = material.rho
            destination.epsilon[row, column] = epsilon0 * material.eps_r
            destination.mu[row, column] = mu0 * material.mu_r
        end
    end
    return destination
end

function homogenize!(
        destination,
        sequence::EHEM.AfterFD,
        relation,
        evaluated,
        model::EarthModel,
        pairs,
        frequency,
        frequency_index::Int
)
    rho = @view evaluated.rho[:, frequency_index]
    eps_r = @view evaluated.eps_r[:, frequency_index]
    mu_r = @view evaluated.mu_r[:, frequency_index]
    air = EarthMaterial(rho[1], eps_r[1], mu_r[1])
    @inbounds for index in eachindex(pairs)
        pair = pairs[index]
        earth = sequence.rule(
            _layout(pair), rho, eps_r, mu_r, model, pair, frequency
        )
        _media!(destination, index, air, earth)
    end
    return destination
end

function homogenize!(
        destination,
        sequence::EHEM.BeforeFD,
        relation,
        static,
        model::EarthModel,
        pairs,
        frequency,
        frequency_index::Int
)
    air = EarthMaterial(static.rho[1], static.eps_r[1], static.mu_r[1])
    @inbounds for index in eachindex(pairs)
        pair = pairs[index]
        reconstructed = sequence.rule(
            _layout(pair), static.rho, static.eps_r, static.mu_r,
            model, pair, frequency
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
        frequency_index::Int
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
