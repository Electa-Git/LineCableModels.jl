# Linear propagation computes each point of the problem's space with its inner formulation,
# from the problem that the point declares with its uncertain values, and each value records
# the estimator and the representation of its uncertainty.
function compute(problem::ParametricProblem, study::LinearError;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    traversed = traverse(problem, study; options)
    values = map(traversed.values) do value
        Engine.retain_gridpoint(value;
            fields = (
                uncertainty = (
                    estimator = :first_order, representation = :dependency_preserving),
                uncertainty_descriptions = _uncertainty_descriptions(LinearError)))
    end
    return LinearErrorResult(values, traversed.failures, traversed.details)
end

# A point of linear propagation fails as a point of the Combinatorial study of its inner
# formulation does.
on_error(study::LinearError, exception) = on_error(Combinatorial(study.inner), exception)

# The formulation of the points of a study of uncertainty is its inner formulation.
points(study::Union{LinearError, MonteCarlo}) = points(study.inner)
