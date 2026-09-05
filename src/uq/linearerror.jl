function compute(problem::ParametricProblem, formulation::LinearError)
    traversed = traverse(problem, formulation)
    return LinearErrorResult(formulation, traversed.values, traversed.details)
end
