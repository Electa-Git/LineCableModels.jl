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

# A parametric problem composes with modal analysis of each completed line-parameter result,
# computed with a core formulation, a space of them or a study whose result space transports
# into modal problems. `options` are the computation options of the parametric run.
function compute(problem::ParametricProblem,
        formulation::Union{ScalarOrGrid{AbstractFormulation}, AbstractTraversal},
        modal::ScalarOrGrid{ModalAnalysisFormulation};
        options::Union{NamedTuple,ComputationOptions}=ComputationOptions(),
        modal_options::Union{NamedTuple,ComputationOptions}=ComputationOptions())
    phase=compute(problem,formulation;options)
    return _modal_results(phase,modal,modal_options)
end

function _modal_results(phase::LineParameters,modal,options)
    return compute(ModalAnalysisProblem(phase),modal;options)
end
function _modal_results(phase::AbstractResultSpace,modal,options)
    return compute(Gridspace{ModalAnalysisProblem}(phase),modal;options)
end
