# Linear propagation traverses a space. `Grid` returns a space as it is and makes a scalar
# problem the space of one point.
function compute(problem::ParametricProblem, formulation::LinearError)
    inner = Combinatorial(formulation.inner;
        options = (retain_details = formulation.options.data.retain_details,))
    traversed = traverse(Grid(problem.space), inner; options = problem.options)
    values = map(traversed.values) do value
        id = get(Commons.observation_gridpoint(value), :id, nothing)
        Engine.retain_gridpoint(value, id;
            fields = (
                uncertainty = (
                    estimator = :first_order, representation = :dependency_preserving),
                uncertainty_descriptions = _uncertainty_descriptions(LinearError)))
    end
    return LinearErrorResult(formulation, values, traversed.details)
end
