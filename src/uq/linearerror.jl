function compute(problem::ParametricProblem, formulation::LinearError)
    ParametricBuilder._traverse(problem, formulation, LinearErrorResult)
end
