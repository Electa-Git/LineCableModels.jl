@inline function _check_operators(
        maps::ModalOperators,
        parameters::LineParameters
)
    expected = size(parameters.Z.values)
    size(maps.voltage) == expected || throw(DimensionMismatch(
        "voltage operators must have dimensions $expected"
    ))
    size(maps.current) == expected || throw(DimensionMismatch(
        "current operators must have dimensions $expected"
    ))
    all(isfinite, maps.voltage) || throw(
        DomainError(maps.voltage, "voltage operators must be finite")
    )
    all(isfinite, maps.current) || throw(
        DomainError(maps.current, "current operators must be finite")
    )
    return maps
end

function computation_options(
        ::Type{LineCableModelsModal},
        record::ComputationOptions
)::ComputationOptions
    options = record.data
    isempty(setdiff(keys(options), (:offdiagonal_tolerance,))) || throw(ArgumentError(
        "unknown modal action options: $(Tuple(keys(options)))"))
    tolerance = get(options, :offdiagonal_tolerance, 1e-6)
    tolerance isa Real && isfinite(tolerance) && tolerance >= 0 || throw(ArgumentError(
        "offdiagonal_tolerance must be finite and nonnegative"))
    return ComputationOptions(; offdiagonal_tolerance = tolerance)
end

function computation_options(
        ::Type{<:ModalTransformationProblem{P}}, options::ComputationOptions
)::ComputationOptions where {T, U, D <: ModalDomain, P <: LineParameters{T, U, D}}
    isempty(options.data) || throw(ArgumentError(
        "inverse modal transformation uses the stored operators and accepts no action options"))
    return ComputationOptions()
end

function _forward(
        parameters::LineParameters{T, U, PhaseDomain, Basis},
        formulation::ModalTransformationFormulation,
        execution::ComputationOptions
) where {T <: Complex, U <: Real, Basis}
    formula = formulation.formula
    workspace = (fallback_frequencies = Int[],)
    maps = modal_operators(formula, parameters, formula.parameters, formula.options, workspace)
    maps isa ModalOperators || throw(ArgumentError("modal_operators must return ModalOperators"))
    _check_operators(maps, parameters)
    voltage = maps.voltage
    current = maps.current
    S = promote_type(T, eltype(voltage), eltype(current))
    dimensions = size(parameters.Z.values)
    impedance = Array{S, 3}(undef, dimensions)
    admittance = Array{S, 3}(undef, dimensions)
    n = dimensions[1]
    source = Matrix{S}(undef, n, n)
    product = similar(source)
    factor = similar(source)
    tolerance = execution.data.offdiagonal_tolerance
    identifier = formula_id(formula)

    @inbounds for frequency in axes(impedance, 3)
        A = @view voltage[:, :, frequency]
        B = @view current[:, :, frequency]
        copyto!(source, @view(parameters.Z.values[:, :, frequency]))
        mul!(product, A, source)
        copyto!(factor, B)
        rdiv!(product, lu!(factor))
        copyto!(@view(impedance[:, :, frequency]), product)

        copyto!(source, @view(parameters.Y.values[:, :, frequency]))
        mul!(product, B, source)
        copyto!(factor, A)
        rdiv!(product, lu!(factor))
        copyto!(@view(admittance[:, :, frequency]), product)

        ratio = offdiagonal_ratio(@view(impedance[:, :, frequency]))
        ratio <= tolerance ||
            @warn(":$identifier transformed Z exceeds the off-diagonal tolerance",
                ratio,
                tolerance)
        ratio = offdiagonal_ratio(@view(admittance[:, :, frequency]))
        ratio <= tolerance ||
            @warn(":$identifier transformed Y exceeds the off-diagonal tolerance",
                ratio,
                tolerance)
    end

    modal = ModalDomain{typeof(maps), AbstractFormulation}(maps, formula)
    return LineParameters(
        modal,
        SeriesImpedance{eltype(impedance), Basis}(impedance),
        ShuntAdmittance{eltype(admittance), Basis}(admittance),
        parameters.f,
        ComputationDetails(merge(parameters.details.data,
            (modal = (identifier = identifier,
                requested = NamedTuple(formulation.definition),
                effective = NamedTuple(formula),
                options = formula.options.data, acceptance = execution.data,
                fallback_frequencies = copy(workspace.fallback_frequencies)),)))
    )
end

function _inverse(
        parameters::LineParameters{
        T, U, D, Basis}
) where {T <: Complex, U <: Real, D <: ModalDomain, Basis}
    maps = _check_operators(parameters.domain.operators, parameters)
    voltage = maps.voltage
    current = maps.current
    S = promote_type(T, eltype(voltage), eltype(current))
    dimensions = size(parameters.Z.values)
    impedance = Array{S, 3}(undef, dimensions)
    admittance = Array{S, 3}(undef, dimensions)
    n = dimensions[1]
    product = Matrix{S}(undef, n, n)
    solution = similar(product)
    factor = similar(product)

    @inbounds for frequency in axes(impedance, 3)
        A = @view voltage[:, :, frequency]
        B = @view current[:, :, frequency]
        Zm = @view parameters.Z.values[:, :, frequency]
        Ym = @view parameters.Y.values[:, :, frequency]
        mul!(product, Zm, B)
        copyto!(factor, A)
        ldiv!(solution, lu!(factor), product)
        copyto!(@view(impedance[:, :, frequency]), solution)

        mul!(product, Ym, A)
        copyto!(factor, B)
        ldiv!(solution, lu!(factor), product)
        copyto!(@view(admittance[:, :, frequency]), solution)
    end

    return LineParameters(
        PhaseDomain,
        SeriesImpedance{eltype(impedance), Basis}(impedance),
        ShuntAdmittance{eltype(admittance), Basis}(admittance),
        parameters.f,
        parameters.details
    )
end

function compute(
        problem::ModalTransformationProblem{P},
        formulation::ModalTransformationFormulation;
        options::Union{NamedTuple,ComputationOptions} = ComputationOptions()
) where {T, U, P <: LineParameters{T, U, PhaseDomain}}
    options = options isa NamedTuple ? ComputationOptions(options) : options
    return compute(LineCableModelsModal(), problem, formulation; options)
end

function compute(
        ::LineCableModelsModal,
        problem::ModalTransformationProblem{P},
        formulation::ModalTransformationFormulation;
        options::Union{NamedTuple,ComputationOptions} = ComputationOptions()
) where {T, U, P <: LineParameters{T, U, PhaseDomain}}
    options = options isa NamedTuple ? ComputationOptions(options) : options
    execution = computation_options(LineCableModelsModal, options)
    validate(problem)
    return _forward(problem.parameters, formulation, execution)
end

function compute(
        problem::ModalTransformationProblem{P};
        options::Union{NamedTuple,ComputationOptions} = ComputationOptions()
) where {T, U, D <: ModalDomain, P <: LineParameters{T, U, D}}
    options = options isa NamedTuple ? ComputationOptions(options) : options
    return compute(LineCableModelsModal(), problem; options)
end

function compute(
        ::LineCableModelsModal,
        problem::ModalTransformationProblem{P};
        options::Union{NamedTuple,ComputationOptions} = ComputationOptions()
) where {T, U, D <: ModalDomain, P <: LineParameters{T, U, D}}
    options = options isa NamedTuple ? ComputationOptions(options) : options
    computation_options(typeof(problem), options)
    validate(problem)
    return _inverse(problem.parameters)
end

function computation_details(
        ::Type{<:ModalTransformationFormulation},
        result::LineParameters
)::ComputationDetails
    return details(result)
end
