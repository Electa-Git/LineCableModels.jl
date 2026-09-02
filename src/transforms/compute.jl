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

@inline function _execution(options::NamedTuple)
    isempty(options) || throw(ArgumentError(
        "the modal backend has no computation options; got $(collect(keys(options)))"
    ))
    return nothing
end

function _forward(
        parameters::LineParameters{T, U, PhaseDomain, Basis},
        formulation::ModalTransformationFormulation
) where {T <: Complex, U <: Real, Basis}
    formula = formulation.formula
    maps = _check_operators(formula(parameters), parameters)
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
    tolerance = formula.assumptions.tolerance
    identifier = formula_id(formula)

    @inbounds for frequency in axes(impedance, 3)
        A = @view voltage[:, :, frequency]
        B = @view current[:, :, frequency]
        copyto!(source, @view(parameters.Z.values[:, :, frequency]))
        reciprocity!(source)
        mul!(product, A, source)
        copyto!(factor, B)
        rdiv!(product, lu!(factor))
        copyto!(@view(impedance[:, :, frequency]), product)

        copyto!(source, @view(parameters.Y.values[:, :, frequency]))
        reciprocity!(source)
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

    modal = ModalDomain{typeof(maps), Formula}(maps, formula)
    return LineParameters(
        modal,
        SeriesImpedance{eltype(impedance), Basis}(impedance),
        ShuntAdmittance{eltype(admittance), Basis}(admittance),
        parameters.f,
        parameters.details
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
        reciprocity!(solution)
        copyto!(@view(impedance[:, :, frequency]), solution)

        mul!(product, Ym, A)
        copyto!(factor, B)
        ldiv!(solution, lu!(factor), product)
        reciprocity!(solution)
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
        options::NamedTuple = (;)
) where {T, U, P <: LineParameters{T, U, PhaseDomain}}
    return compute(LineCableModelsModal(), problem, formulation; options)
end

function compute(
        ::LineCableModelsModal,
        problem::ModalTransformationProblem{P},
        formulation::ModalTransformationFormulation;
        options::NamedTuple = (;)
) where {T, U, P <: LineParameters{T, U, PhaseDomain}}
    _execution(options)
    validate(problem)
    return _forward(problem.parameters, formulation)
end

function compute(
        problem::ModalTransformationProblem{P};
        options::NamedTuple = (;)
) where {T, U, D <: ModalDomain, P <: LineParameters{T, U, D}}
    return compute(LineCableModelsModal(), problem; options)
end

function compute(
        ::LineCableModelsModal,
        problem::ModalTransformationProblem{P};
        options::NamedTuple = (;)
) where {T, U, D <: ModalDomain, P <: LineParameters{T, U, D}}
    _execution(options)
    validate(problem)
    return _inverse(problem.parameters)
end

computation_owner(::ModalTransformationFormulation) = LineCableModelsModal

function computation_details(
        ::Val{LineCableModelsModal},
        result::LineParameters
)::ComputationDetails
    return details(result)
end
