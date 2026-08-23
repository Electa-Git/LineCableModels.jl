struct ConsoleVerbosityLogger{L <: AbstractLogger, V <: NamedTuple} <: AbstractLogger
    console::L
    levels::V
end

Logging.min_enabled_level(::ConsoleVerbosityLogger) = Logging.Debug
Logging.catch_exceptions(::ConsoleVerbosityLogger) = false

function _log_level(level::Integer)
    level == 0 ? Logging.Warn :
    level == 1 ? Logging.Info : Logging.Debug
end

function _log_key(source::Module)
    path = Base.fullname(source)
    return isempty(path) ? :default : first(path)
end

function Logging.shouldlog(
        logger::ConsoleVerbosityLogger,
        level,
        source,
        group,
        id
)
    selected = get(logger.levels, _log_key(source), logger.levels.default)
    return level >= _log_level(selected) &&
           Logging.shouldlog(logger.console, level, source, group, id)
end

function Logging.handle_message(
        logger::ConsoleVerbosityLogger,
        level,
        message,
        source,
        group,
        id,
        file,
        line;
        kwargs...
)
    return Logging.handle_message(
        logger.console, level, message, source, group, id, file, line; kwargs...
    )
end

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

function _basis_result(parameters, ::AnalyticalInput, ::Val{:per_length})
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

_input(problem::LineParametersProblem, ::AnalyticalFormulation) = AnalyticalInput(problem)

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
        SeriesImpedance{Complex{T}, :per_length}(Zout),
        ShuntAdmittance{Complex{T}, :per_length}(Yout),
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

Compute frequency-dependent line parameters for one materialized problem.

The problem and formulation are validated, immutable solver input is flattened once,
and operating-temperature resistivity is calculated into a local array before the
frequency loop. Neither the problem nor its cable design is mutated.

# Keywords

- `options=(verbosity=(default=0,), output_basis=:per_length)`: Execution
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
        options = (;)
)
    execution = computation_options(options)
    console = ConsoleLogger(stderr, Logging.Debug)
    logger = ConsoleVerbosityLogger(console, execution.verbosity)
    return with_logger(logger) do
        _compute_analytical(problem, formulation, execution)
    end
end

function compute(
        problem::CableConstantsProblem,
        formulation::AnalyticalFormulation;
        options = (;)
)
    formulation.options.output isa Val{:parameters} || throw(ArgumentError(
        "CableConstantsProblem supports only formulation output=:parameters",
    ))
    execution = computation_options(options)
    _output_basis(execution) === :per_length || throw(ArgumentError(
        "CableConstantsProblem supports only output_basis=:per_length",
    ))
    return DataModel._compute_cable_constants(
        problem.design; S = problem.separation, rho_e = problem.earth_resistivity
    )
end

function DataModel._base_parameters(design::CableDesign)
    compute(CableConstantsProblem(design), Formulation())
end

@inline function compute_earth_return_matrix!(
        destination::AbstractMatrix{Complex{T}},
        cables::AbstractVector{Int},
        input::AnalyticalInput{T},
        earth,
        frequency::Int,
        formulation
) where {T <: Real}
    rho = @view earth.rho[:, frequency]
    epsilon = @view earth.epsilon[:, frequency]
    mu = @view earth.mu[:, frequency]
    s = input.jω[frequency]
    @inbounds for column in eachindex(cables), row in eachindex(cables)

        left = cables[row]
        right = cables[column]
        separation = input.horz_sep[left, right]
        heights = (input.vert[left], input.vert[right])
        kind = row == column ? Val(:self) : Val(:mutual)
        destination[row, column] = formulation(
            kind, heights, separation, rho, epsilon, mu, s
        )
    end
    return destination
end

function compute_impedance_matrix!(
        destination::AbstractMatrix{Complex{T}},
        input::AnalyticalInput{T},
        rho_cond::AbstractVector{T},
        earth,
        frequency::Int,
        formulation::AnalyticalFormulation,
        trace
) where {T <: Real}
    fill!(destination, zero(Complex{T}))
    indices, cables = _cable_indices(input)
    earth_matrix = Matrix{Complex{T}}(undef, input.n_cables, input.n_cables)
    compute_earth_return_matrix!(
        earth_matrix, cables, input, earth, frequency,
        formulation.earth_impedance
    )
    _stash!(_trace_target(trace, :Zg), frequency, earth_matrix)
    s = input.jω[frequency]

    @inbounds for cable in 1:input.n_cables
        conductors = indices[cable]
        count = length(conductors)
        for position in count:-1:1
            index = conductors[position]
            inner = input.r_in[index]
            outer = input.r_ext[index]
            rho = rho_cond[index]
            mu = input.mu_cond[index]
            outside = formulation.internal_impedance(:outer, inner, outer, rho, mu, s)
            inside = position < count ?
                     formulation.internal_impedance(
                :inner,
                input.r_in[conductors[position + 1]],
                input.r_ext[conductors[position + 1]],
                rho_cond[conductors[position + 1]],
                input.mu_cond[conductors[position + 1]],
                s
            ) : zero(outside)
            mutual = formulation.internal_impedance(:mutual, inner, outer, rho, mu, s)
            insulation = formulation.insulation_impedance(
                outer, input.r_ins_ext[index], input.mu_ins[index], s
            )
            loop = outside + inside + insulation
            if position > 1
                for row in 1:(position - 1), column in 1:(position - 1)

                    destination[conductors[row], conductors[column]] += loop - 2 * mutual
                end
                for row in 1:(position - 1)
                    destination[conductors[position], conductors[row]] += loop - mutual
                    destination[conductors[row], conductors[position]] += loop - mutual
                end
            end
            destination[index, index] += loop
        end
    end
    _stash!(_trace_target(trace, :Zin), frequency, destination)

    @inbounds for cable in 1:input.n_cables
        conductors = indices[cable]
        self = earth_matrix[cable, cable]
        for row in conductors, column in conductors

            destination[row, column] += self
        end
    end
    @inbounds for left in 1:(input.n_cables - 1)
        for right in (left + 1):input.n_cables
            mutual = earth_matrix[left, right]
            for row in indices[left], column in indices[right]

                destination[row, column] += mutual
                destination[column, row] += mutual
            end
        end
    end
    _stash!(_trace_target(trace, :Z), frequency, destination)
    return destination
end

function compute_admittance_matrix!(
        destination::AbstractMatrix{Complex{T}},
        input::AnalyticalInput{T},
        earth,
        frequency::Int,
        formulation::AnalyticalFormulation,
        trace
) where {T <: Real}
    fill!(destination, zero(Complex{T}))
    indices, cables = _cable_indices(input)
    earth_matrix = Matrix{Complex{T}}(undef, input.n_cables, input.n_cables)
    compute_earth_return_matrix!(
        earth_matrix, cables, input, earth, frequency,
        formulation.earth_admittance
    )
    _stash!(_trace_target(trace, :Pg), frequency, earth_matrix)
    s = input.jω[frequency]

    @inbounds for cable in 1:input.n_cables
        conductors = indices[cable]
        count = length(conductors)
        coefficients = Vector{Complex{T}}(undef, count)
        for component in 1:count
            coefficients[component] = InsulationAdmittance.potential_coefficient(
                formulation.insulation_admittance, input, conductors[component], s
            )
        end
        tails = Vector{Complex{T}}(undef, count)
        tails[count] = coefficients[count]
        for gap in (count - 1):-1:1
            tails[gap] = coefficients[gap] + tails[gap + 1]
        end
        for row in 1:count, column in 1:count

            destination[conductors[row], conductors[column]] += tails[max(row, column)]
        end
    end
    _stash!(_trace_target(trace, :Pin), frequency, destination)

    @inbounds for cable in 1:input.n_cables
        self = earth_matrix[cable, cable]
        for row in indices[cable], column in indices[cable]

            destination[row, column] += self
        end
    end
    @inbounds for left in 1:(input.n_cables - 1)
        for right in (left + 1):input.n_cables
            mutual = earth_matrix[left, right]
            for row in indices[left], column in indices[right]

                destination[row, column] += mutual
                destination[column, row] += mutual
            end
        end
    end
    _stash!(_trace_target(trace, :P), frequency, destination)
    return destination
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
