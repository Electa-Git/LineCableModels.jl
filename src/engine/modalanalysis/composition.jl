# A third positional argument composes one line-parameter computation with the modal
# analysis of its result, whatever formulation computed it. The argument is a modal
# formulation or the identifier of one.
function compute(problem::LineParametersProblem,formulation::AbstractFormulation,
        modal::Union{Symbol,ModalAnalysisFormulation};
        options::Union{NamedTuple,ComputationOptions}=ComputationOptions(),
        modal_options::Union{NamedTuple,ComputationOptions}=ComputationOptions())
    phase=compute(problem,formulation;options)
    return compute(ModalAnalysisProblem(phase),ModalAnalysisFormulation(modal);
        options=modal_options)
end
