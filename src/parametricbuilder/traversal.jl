"""
$(TYPEDSIGNATURES)

Evaluate one core problem per Gridspace point and optionally retain one details
record per result.

# Arguments

- `problem`: Parametric problem containing a nonempty finite problem space.
- `formulation`: Higher-order formulation with an `inner` formulation and
  `retain_details` computation option.

# Returns

- A named tuple containing one concrete `values` vector and either an empty
  details tuple or `(points=records,)`.
"""
function traverse(problem::ParametricProblem, formulation)
    point_count = length(problem.space)
    point_count > 0 || throw(ArgumentError(
        "higher-order problem space must contain at least one core problem",
    ))

    point_source = points(problem.space)
    first_item = iterate(point_source)
    first_item === nothing && throw(DimensionMismatch(
        "problem-space iteration ended before its declared cardinality",
    ))
    first_point, state = first_item
    first_problem = materialize(first_point)
    first_result = compute(
        first_problem,
        formulation.inner;
        options = problem.options
    )
    check_core_result(typeof(first_result))
    values = Vector{typeof(first_result)}(undef, point_count)
    values[1] = first_result

    details_owner = formulation.options.retain_details ?
                    computation_owner(formulation.inner) : nothing
    retained = if formulation.options.retain_details
        first_record = computation_details(Val(details_owner), first_result)
        records = Vector{typeof(first_record)}(undef, point_count)
        records[1] = first_record
        records
    else
        nothing
    end

    for index in 2:point_count
        item = iterate(point_source, state)
        item === nothing && throw(DimensionMismatch(
            "problem-space iteration ended before its declared cardinality",
        ))
        point, state = item
        core_result = compute(
            materialize(point),
            formulation.inner;
            options = problem.options
        )
        typeof(core_result) === eltype(values) || throw(ArgumentError(
            "higher-order computation produced inconsistent core result types",
        ))
        values[index] = core_result

        if retained !== nothing
            record = computation_details(Val(details_owner), core_result)
            typeof(record) === eltype(retained) || throw(ArgumentError(
                "higher-order computation produced inconsistent details record types",
            ))
            retained[index] = record
        end
    end
    iterate(point_source, state) === nothing || throw(DimensionMismatch(
        "problem-space iteration exceeded its declared cardinality",
    ))

    retained_details = retained === nothing ? (;) : (points = retained,)
    return (; values, details = retained_details)
end

function compute(problem::ParametricProblem, formulation::Combinatorial)
    traversed = traverse(problem, formulation)
    return ParametricResult(formulation, traversed.values, traversed.details)
end
