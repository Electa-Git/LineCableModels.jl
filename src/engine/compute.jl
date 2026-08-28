_modal_result(parameters, ::Nothing) = parameters

function _modal_result(parameters, transform::AbstractTransformFormulation)
    _, transformed = transform(parameters)
    return transformed
end

_basis_result(parameters, ::LineParametersWorkspace, ::Val{:pul}) = parameters

function _basis_result(
        parameters::LineParameters{T, U, D},
        workspace::LineParametersWorkspace,
        ::Val{:total}
) where {T, U, D}
    line_length = workspace.normalized.line_length
    impedance = parameters.Z.values .* line_length
    admittance = parameters.Y.values .* line_length
    return LineParameters(
        D,
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
    input = workspace.normalized
    prepared = workspace.prepared
    buffers = workspace.buffers
    Zbuffer = buffers.Zbuffer
    Pbuffer = buffers.Pbuffer
    Zprimitive = buffers.Zprimitive
    Pprimitive = buffers.Pprimitive
    Pinverse = buffers.Pinverse
    reduced = buffers.reduced
    reduced_inverse = buffers.reduced_inverse
    Zout = buffers.Zout
    Yout = buffers.Yout
    permutation = prepared.permutation
    reordered_map = prepared.reordered_map
    kron_map = prepared.kron_map

    @info "Starting line parameters computation"
    for frequency in 1:input.n_frequencies
        compute_impedance_matrix!(
            Zprimitive,
            workspace,
            frequency,
            formulation
        )
        compute_admittance_matrix!(
            Pprimitive,
            workspace,
            frequency,
            formulation
        )
        _reorder_into!(Zbuffer, Zprimitive, permutation)
        _reorder_into!(Pbuffer, Pprimitive, permutation)
        if formulation.options.reduce_bundle
            merge_bundles!(Zbuffer, reordered_map)
            merge_bundles!(Pbuffer, reordered_map)
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
            kronify!(Zbuffer, kron_map, reduced)
            reciprocity!(reduced)
            formulation.options.ideal_transposition && ideal_transposition!(reduced)
            @views Zout[:, :, frequency] .= reduced

            kronify!(Pbuffer, kron_map, reduced)
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

function _transform(
        parameters::LineParameters,
        formulation::LineParametersFormulation
)
    return _modal_result(parameters, formulation.methods.modal_transform)
end

_retained_details(::LineParametersWorkspace{<:Real, <:NamedTuple, <:NamedTuple,
                                            <:NamedTuple, Nothing}) = (;)

function _retained_details(workspace::LineParametersWorkspace)
    capture = workspace.capture
    capture === nothing && return (;)
    input = workspace.normalized
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
        domain(parameters),
        parameters.Z,
        parameters.Y,
        parameters.f,
        retained
    )
end

function _compute_line_parameters(
        engine::LineCableModelsEngine,
        problem::LineParametersProblem,
        formulation::LineParametersFormulation,
        execution::NamedTuple
)
    validate(problem)
    workspace = LineParametersWorkspace(engine, problem, formulation, execution)
    parameters = _solve!(workspace, formulation)
    parameters = _transform(parameters, formulation)
    return _finish(parameters, workspace, formulation, execution)
end

"""
$(TYPEDSIGNATURES)

Compute line parameters with the native engine and default formulation.

# Arguments

- `problem`: Completed line-parameter problem.

# Keywords

- `options`: Native computation options.

# Returns

- One [`LineParameters`](@ref) result.
"""
function compute(
        problem::LineParametersProblem;
        options::NamedTuple = (;)
)
    return compute(LineCableModelsEngine(), problem, Formulation(); options)
end

"""
$(TYPEDSIGNATURES)

Compute frequency-dependent line parameters with the native engine.

The completed physical system is normalized once into a backend-owned
workspace. All reusable numerical storage is allocated before the frequency
loop. `trace=true` retains completed intermediate matrices under
`details(result).trace`; it does not change the result type.

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
    return compute(LineCableModelsEngine(), problem, formulation; options)
end

function compute(
        engine::LineCableModelsEngine,
        problem::LineParametersProblem,
        formulation::LineParametersFormulation = Formulation();
        options::NamedTuple = (;)
)
    execution = computation_options(Val(LineCableModelsEngine), options)
    console = ConsoleLogger(stderr, Logging.Debug)
    logger = ConsoleVerbosityLogger(console, execution.verbosity)
    return with_logger(logger) do
        _compute_line_parameters(engine, problem, formulation, execution)
    end
end

computation_owner(::LineParametersFormulation) = LineCableModelsEngine

function computation_details(
        ::Val{LineCableModelsEngine},
        result::LineParameters
)::ComputationDetails
    return details(result)
end

computation_details(
    ::Val{LineCableModelsEngine},
    ::DataModel.CableConstants
)::ComputationDetails = (;)

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
    evaluated = EarthProperties.evaluate(
        formulation.methods.earth_properties,
        input.earth,
        input.freq
    )
    formulation.methods.equivalent_earth === nothing && return evaluated
    return formulation.methods.equivalent_earth(evaluated, input.earth)
end
