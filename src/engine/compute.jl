_trace_buffers(::Type, ::AnalyticalInput, ::Val{:parameters}) = nothing

function _trace_buffers(::Type{T}, input::AnalyticalInput{T}, ::Val{:trace}) where {T}
    n, nc, nf = input.n_phases, input.n_cables, input.n_frequencies
    return (
        Zin = Array{Complex{T}, 3}(undef, n, n, nf),
        Pin = Array{Complex{T}, 3}(undef, n, n, nf),
        Zg = Array{Complex{T}, 3}(undef, nc, nc, nf),
        Pg = Array{Complex{T}, 3}(undef, nc, nc, nf),
        Z = Array{Complex{T}, 3}(undef, n, n, nf),
        P = Array{Complex{T}, 3}(undef, n, n, nf)
    )
end

_solve_result(parameters, ::AnalyticalInput, ::Nothing, ::Val{:parameters}) = parameters

function _solve_result(parameters, input::AnalyticalInput, trace, ::Val{:trace})
    return LineParametersTrace(
        parameters, copy(input.freq), copy(input.phase_map), copy(input.cable_map),
        trace.Zin, trace.Pin, trace.Zg, trace.Pg, trace.Z, trace.P
    )
end

_modal_result(parameters, ::Nothing) = parameters

function _modal_result(parameters, transform::AbstractTransformFormulation)
    _, transformed = transform(parameters)
    return transformed
end

function _basis_result(parameters, ::AnalyticalInput, ::Val{:pul})
    parameters
end

function _basis_result(
        parameters::LineParameters{T, U, D},
        input::AnalyticalInput,
        ::Val{:total}
) where {T, U, D}
    impedance = parameters.Z.values .* input.line_length
    admittance = parameters.Y.values .* input.line_length
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

_trace_target(::Nothing, ::Symbol) = nothing
_trace_target(trace, name::Symbol) = getproperty(trace, name)

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

function _operating_resistivity(
        input::AnalyticalInput{T},
        problem::LineParametersProblem{T},
        formulation::AnalyticalFormulation
) where {T <: Real}
    rho = copy(input.rho0_cond)
    formulation.options.temperature_correction || return rho
    @inbounds for index in eachindex(rho)
        rho[index] *= DataModel.temperature_factor(
            input.alpha_cond[index], problem.temperature, input.T0_cond[index]
        )
    end
    return rho
end

_input(problem::LineParametersProblem, formulation::AnalyticalFormulation) =
    AnalyticalInput(problem, formulation)

function _prepare(
        problem::LineParametersProblem,
        input::AnalyticalInput,
        formulation::AnalyticalFormulation
)
    return _operating_resistivity(input, problem, formulation),
    _earth_data(formulation, input)
end

function _solve(
        input::AnalyticalInput{T},
        rho_cond::AbstractVector{T},
        earth,
        formulation::AnalyticalFormulation
) where {T <: Real}
    n, nf = input.n_phases, input.n_frequencies

    Zbuffer = Matrix{Complex{T}}(undef, n, n)
    Pbuffer = Matrix{Complex{T}}(undef, n, n)
    Zprimitive = similar(Zbuffer)
    Pprimitive = similar(Pbuffer)
    Pinverse = similar(Pbuffer)
    output = formulation.options.output
    trace = _trace_buffers(T, input, output)

    permutation, reordered_map, kron_map = _reduction_map(input.phase_map, formulation)
    nkeep = kron_map === nothing ? n : count(!=(0), kron_map)
    Zout = Array{Complex{T}, 3}(undef, nkeep, nkeep, nf)
    Yout = Array{Complex{T}, 3}(undef, nkeep, nkeep, nf)
    reduced = Matrix{Complex{T}}(undef, nkeep, nkeep)
    reduced_inverse = similar(reduced)
    identity_full = Matrix{Complex{T}}(I, n, n)
    identity_reduced = Matrix{Complex{T}}(I, nkeep, nkeep)

    @info "Starting line parameters computation"
    for frequency in 1:nf
        compute_impedance_matrix!(
            Zprimitive, input, rho_cond, earth, frequency, formulation, trace
        )
        compute_admittance_matrix!(
            Pprimitive, input, earth, frequency, formulation, trace
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
            ldiv!(Pinverse, factorization, identity_full)
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
            ldiv!(reduced_inverse, factorization, identity_reduced)
            reduced_inverse .*= input.jω[frequency]
            reciprocity!(reduced_inverse)
            formulation.options.ideal_transposition && ideal_transposition!(reduced_inverse)
            @views Yout[:, :, frequency] .= reduced_inverse
        end
    end

    parameters = LineParameters(
        PhaseDomain,
        SeriesImpedance{Complex{T}, :pul}(Zout),
        ShuntAdmittance{Complex{T}, :pul}(Yout),
        input.freq
    )
    return parameters, trace
end

function _transform(parameters::LineParameters, formulation::AnalyticalFormulation)
    return _modal_result(parameters, formulation.modal_transform)
end

function _finish(
        parameters::LineParameters,
        input::AnalyticalInput,
        trace,
        formulation::AnalyticalFormulation,
        execution::NamedTuple
)
    parameters = _basis_result(parameters, input, execution.output_basis)
    @info "Line parameters computation completed successfully"
    return _solve_result(parameters, input, trace, formulation.options.output)
end

function _compute_analytical(
        problem::LineParametersProblem,
        formulation::AnalyticalFormulation,
        execution::NamedTuple
)
    validate(problem)
    input = _input(problem, formulation)
    rho_cond, earth = _prepare(problem, input, formulation)
    parameters, trace = _solve(input, rho_cond, earth, formulation)
    parameters = _transform(parameters, formulation)
    return _finish(parameters, input, trace, formulation, execution)
end

"""
$(TYPEDSIGNATURES)

Compute line parameters with the default analytical formulation.

# Arguments

- `problem`: Completed line-parameter problem.

# Keywords

- `options`: Computation options normalized for `AnalyticalFormulation`.

# Returns

One `LineParameters` result.
"""
function compute(
        problem::LineParametersProblem;
        options::NamedTuple = (;)
)
    return compute(problem, Formulation(); options)
end

"""
$(TYPEDSIGNATURES)

Compute frequency-dependent line parameters for one materialised problem.

The problem and formulation are validated, immutable solver input is flattened once,
and operating-temperature resistivity is calculated into a local array before the
frequency loop. Neither the problem nor its cable design is mutated.

# Keywords

- `options=(verbosity=(default=0,), output_basis=:pul)`: Execution
  verbosity and output basis.

Trace output is selected by the formulation:

```julia
Formulation(:analytical; options = (output = :trace,))
```

# Returns

- [`LineParameters`](@ref) for `output=:parameters`.
- [`LineParametersTrace`](@ref) for `output=:trace`.
"""
function compute(
        problem::LineParametersProblem,
        formulation::AnalyticalFormulation;
        options::NamedTuple = (;)
)
    execution = computation_options(Val(AnalyticalFormulation), options)
    console = ConsoleLogger(stderr, Logging.Debug)
    logger = ConsoleVerbosityLogger(console, execution.verbosity)
    return with_logger(logger) do
        _compute_analytical(problem, formulation, execution)
    end
end

computation_details(
    ::Val{AnalyticalFormulation},
    ::Union{LineParameters, LineParametersTrace, DataModel.CableConstants}
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
- `formulation`: Analytical formulation used by `compute`.
- `options`: Computation options forwarded to `compute`.
- `position`: Cable placement in metres.
- `line_length`: Physical line length in metres.
- `earth_props`: Static earth model.
- `temperature`: Operating temperature in °C.
- `frequency`: Analysis frequency in Hz.

# Returns

A [`CableConstants`](@ref) value observed from the resulting line parameters.
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

function _cable_indices(input)
    indices = [Int[] for _ in 1:input.n_cables]
    @inbounds for index in 1:input.n_phases
        push!(indices[input.cable_map[index]], index)
    end
    return indices, first.(indices)
end

function _earth_data(formulation::AnalyticalFormulation, input::AnalyticalInput)
    evaluated = EarthProperties.evaluate(
        formulation.earth_properties, input.earth, input.freq
    )
    formulation.equivalent_earth === nothing && return evaluated
    return formulation.equivalent_earth(evaluated, input.earth)
end
