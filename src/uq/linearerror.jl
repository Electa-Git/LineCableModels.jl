# Linear propagation traverses a space. `Grid` returns a space as it is and makes a scalar
# problem the space of one point. The `Combinatorial` study records the failed points.
function compute(problem::ParametricProblem, formulation::LinearError)
    inner = Combinatorial(formulation.inner;
        options = (retain_details = formulation.options.data.retain_details,))
    traversed = traverse(Grid(problem.space), inner; options = problem.options)
    values = map(traversed.values) do value
        Engine.retain_gridpoint(value;
            fields = (
                uncertainty = (
                    estimator = :first_order, representation = :dependency_preserving),
                uncertainty_descriptions = _uncertainty_descriptions(LinearError)))
    end
    return LinearErrorResult(values, traversed.failures, traversed.details)
end
