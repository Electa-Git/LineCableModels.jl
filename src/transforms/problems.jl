"""
$(TYPEDEF)

Define a bidirectional transformation of fully coupled line parameters.

The source domain determines the direction. Phase-domain input is transformed
to a modal domain by the selected formulation. Modal-domain input is returned
to the phase domain with its stored operators.

$(TYPEDFIELDS)
"""
struct ModalTransformationProblem{P <: LineParameters} <: AbstractProblemDefinition
    "Source line parameters."
    parameters::P

    function ModalTransformationProblem{P}(parameters::P) where {P <: LineParameters}
        return validate(new{P}(parameters))
    end
end

function validate(problem::ModalTransformationProblem)
    validate(problem.parameters)
    if problem.parameters.domain isa ModalDomain
        maps = problem.parameters.domain.operators
        expected = size(problem.parameters.Z.values)
        size(maps.voltage) == expected || throw(DimensionMismatch(
            "ModalTransformationProblem.parameters.domain.operators.voltage must " *
            "have size $expected; received $(size(maps.voltage))"
        ))
        size(maps.current) == expected || throw(DimensionMismatch(
            "ModalTransformationProblem.parameters.domain.operators.current must " *
            "have size $expected; received $(size(maps.current))"
        ))
        all(isfinite, maps.voltage) || throw(DomainError(
            maps.voltage,
            "ModalTransformationProblem voltage operators must be finite"
        ))
        all(isfinite, maps.current) || throw(DomainError(
            maps.current,
            "ModalTransformationProblem current operators must be finite"
        ))
    end
    return problem
end

function ModalTransformationProblem(
        parameters::LineParameters{T, U, PhaseDomain}
) where {T, U}
    return ModalTransformationProblem{typeof(parameters)}(parameters)
end

function ModalTransformationProblem(
        parameters::LineParameters{T, U, D}
) where {T, U, D <: ModalDomain}
    return ModalTransformationProblem{typeof(parameters)}(parameters)
end

Base.eltype(problem::ModalTransformationProblem) = eltype(problem.parameters)
Base.eltype(::Type{ModalTransformationProblem{P}}) where {P} = eltype(P)
