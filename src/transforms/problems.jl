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
