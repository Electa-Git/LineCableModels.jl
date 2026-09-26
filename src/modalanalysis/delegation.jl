# A third positional argument composes two completed public computations.
# Each upstream wrapper keeps its original no-modal execution path.
const _ModalSelection=Union{ModalAnalysisFormulation,
    Gridspace{<:ModalAnalysisFormulation}}

function compute(problem,formulation,modal::_ModalSelection;
        options::Union{NamedTuple,ComputationOptions}=ComputationOptions(),
        modal_options::Union{NamedTuple,ComputationOptions}=ComputationOptions())
    phase=compute(problem,formulation;options)
    return _modal_results(phase,modal,modal_options)
end

function compute(problem::ParametricProblem,
        formulation::Union{Combinatorial,Gridspace{<:AbstractFormulation}},
        modal::_ModalSelection;
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
function _modal_results(phase::AbstractVector{<:LineParameters},modal,options)
    problems=Gridspace{ModalAnalysisProblem}(ModalAnalysisProblem,(Grid(phase),))
    return compute(problems,modal;options).values
end
