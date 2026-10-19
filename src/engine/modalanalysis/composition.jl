# A third positional argument composes a line-parameter computation with modal
# analysis of its result, whatever formulation computed it. The argument is a modal
# formulation, a space of them or the identifier of one.
function compute(problem::ScalarOrGrid{LineParametersProblem},
        formulation,modal::Union{Symbol,ScalarOrGrid{ModalAnalysisFormulation}};
        options::Union{NamedTuple,ComputationOptions}=ComputationOptions(),
        modal_options::Union{NamedTuple,ComputationOptions}=ComputationOptions())
    phase=compute(problem,formulation;options)
    return _modal_results(phase,ModalAnalysisFormulation(modal),modal_options)
end

# A parametric study composes with modal analysis of each completed line-parameter result.
function compute(problem::ParametricProblem,
        formulation::ScalarOrGrid{AbstractFormulation},
        modal::ScalarOrGrid{ModalAnalysisFormulation};
        modal_options::Union{NamedTuple,ComputationOptions}=ComputationOptions())
    phase=compute(problem,formulation)
    return _modal_results(phase,modal,modal_options)
end

function _modal_results(phase::LineParameters,modal,options)
    return compute(ModalAnalysisProblem(phase),modal;options)
end
function _modal_results(phase::AbstractResultSpace,modal,options)
    return compute(Gridspace{ModalAnalysisProblem}(phase),modal;options)
end
