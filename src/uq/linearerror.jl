function compute(problem::ParametricProblem, formulation::LinearError)
    traversed = traverse(problem, formulation)
    values=map(traversed.values) do value
        id=get(Grammar.observation_gridpoint(value),:id,nothing)
        Engine.retain_gridpoint(value,id;
            fields=(uncertainty=(estimator=:first_order,representation=:dependency_preserving),
                uncertainty_descriptions=_uncertainty_descriptions(LinearError)))
    end
    return LinearErrorResult(formulation, values, traversed.details)
end
