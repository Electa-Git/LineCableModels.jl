function _evaluate_core(
        point,
        descriptors,
        affine_maps,
        latent,
        formulation::PolynomialChaos,
        core_options
)
    evaluation_count = size(latent, 1)
    physical_values = _physical_values(descriptors, affine_maps, latent, 1)
    # The release API stages supplied leaves into top-level arguments before
    # invoking the point builder; this preserves its established Tuple method.
    realized_arguments = realize_arguments(point, physical_values)
    core_problem = realize(point, realized_arguments)
    first_result = compute(core_problem, formulation.inner; options = core_options)
    results = Vector{typeof(first_result)}(undef, evaluation_count)
    results[1] = first_result

    retained = nothing
    if formulation.options.retain_details
        first_detail = computation_details(typeof(formulation.inner), first_result)
        retained = Vector{typeof(first_detail)}(undef, evaluation_count)
        retained[1] = first_detail
    end

    for evaluation in 2:evaluation_count
        physical_values = _physical_values(
            descriptors, affine_maps, latent, evaluation)
        realized_arguments = realize_arguments(point, physical_values)
        core_problem = realize(point, realized_arguments)
        result = compute(core_problem, formulation.inner; options = core_options)
        typeof(result) === eltype(results) || throw(ArgumentError(
            "polynomial-chaos evaluations produced incompatible core result types",
        ))
        results[evaluation] = result
        if retained !== nothing
            record = computation_details(typeof(formulation.inner), result)
            typeof(record) === eltype(retained) || throw(ArgumentError(
                "polynomial-chaos evaluations produced incompatible details types",
            ))
            retained[evaluation] = record
        end
    end
    return results, retained
end

function _validation_latent(
        formulation::PolynomialChaos,
        point_index::Int,
        dimension::Int
)
    seed = formulation.validation_seed ⊻
           (UInt64(point_index - 1) * 0x9e3779b97f4a7c15)
    rng = Random.Xoshiro(seed)
    if formulation.distribution === :normal
        latent = Base.randn(rng, formulation.validation_points, dimension)
        all(isfinite, latent) || throw(ArgumentError(
            "polynomial-chaos validation coordinates must be finite",
        ))
        return latent
    end
    latent = 2 .* Base.rand(rng, formulation.validation_points, dimension) .- 1
    all(isfinite, latent) || throw(ArgumentError(
        "polynomial-chaos validation coordinates must be finite",
    ))
    return latent
end

function _polynomial_chaos_point(plan, formulation::PolynomialChaos, core_options,
        point_index::Int)
    univariate, multivariate = _coordinate_bases(formulation, plan.dimension)
    PolyChaos.dim(multivariate) == plan.terms || throw(DimensionMismatch(
        "PolyChaos basis term count does not match the cost preflight",
    ))
    latent, product_weights =
        _tensor_rule(multivariate, plan.collocation, plan.dimension)
    affine_maps = _affine_maps(
        plan.descriptors, plan.active, univariate, formulation.distribution)
    collocation_results, collocation_details = _evaluate_core(
        plan.point,
        plan.descriptors,
        affine_maps,
        latent,
        formulation,
        core_options
    )
    validation_latent = _validation_latent(
        formulation, point_index, plan.dimension)
    validation_results, validation_details = _evaluate_core(
        plan.point,
        plan.descriptors,
        affine_maps,
        validation_latent,
        formulation,
        core_options
    )

    polynomial_values = PolyChaos.evaluate(latent, multivariate)
    validation_polynomial_values = PolyChaos.evaluate(validation_latent, multivariate)
    size(polynomial_values) == (plan.terms, plan.collocation) || throw(
        DimensionMismatch("PolyChaos collocation evaluation matrix has an invalid shape"),
    )
    size(validation_polynomial_values) == (plan.terms, plan.validation) || throw(
        DimensionMismatch("PolyChaos validation evaluation matrix has an invalid shape"),
    )
    all(isfinite, polynomial_values) && all(isfinite, validation_polynomial_values) ||
        throw(ArgumentError("PolyChaos basis evaluations must be finite"))
    aggregate = fit_expansion(
        multivariate,
        polynomial_values,
        product_weights,
        collocation_results,
        validation_polynomial_values,
        validation_results,
        formulation,
        point_index,
        plan.dimension
    )
    retained = collocation_details === nothing ? nothing : (
        collocation = collocation_details,
        validation = validation_details,
    )
    return merge(aggregate, (; details = retained))
end

function compute(problem::ParametricProblem, formulation::PolynomialChaos)
    plans = _preflight(problem, formulation)
    first_aggregate = _polynomial_chaos_point(
        first(plans), formulation, problem.options, 1)
    point_count = length(plans)
    values = Vector{typeof(first_aggregate.representation)}(undef, point_count)
    expansion_values = Vector{typeof(first_aggregate.expansion)}(undef, point_count)
    stats = Vector{typeof(first_aggregate.statistics)}(undef, point_count)
    validation_values = Vector{typeof(first_aggregate.validation)}(undef, point_count)
    retained = formulation.options.retain_details ?
               Vector{typeof(first_aggregate.details)}(undef, point_count) : nothing

    values[1] = first_aggregate.representation
    expansion_values[1] = first_aggregate.expansion
    stats[1] = first_aggregate.statistics
    validation_values[1] = first_aggregate.validation
    retained === nothing || (retained[1] = first_aggregate.details)

    for point_index in 2:point_count
        aggregate = _polynomial_chaos_point(
            plans[point_index], formulation, problem.options, point_index)
        typeof(aggregate.representation) === eltype(values) || throw(ArgumentError(
            "polynomial-chaos outer points produced incompatible core result types",
        ))
        typeof(aggregate.expansion) === eltype(expansion_values) || throw(ArgumentError(
            "polynomial-chaos outer points produced incompatible expansion types",
        ))
        typeof(aggregate.statistics) === eltype(stats) || throw(ArgumentError(
            "polynomial-chaos outer points produced incompatible statistics types",
        ))
        typeof(aggregate.validation) === eltype(validation_values) || throw(ArgumentError(
            "polynomial-chaos outer points produced incompatible validation types",
        ))
        values[point_index] = aggregate.representation
        expansion_values[point_index] = aggregate.expansion
        stats[point_index] = aggregate.statistics
        validation_values[point_index] = aggregate.validation
        if retained !== nothing
            typeof(aggregate.details) === eltype(retained) || throw(ArgumentError(
                "polynomial-chaos outer points produced incompatible details types",
            ))
            retained[point_index] = aggregate.details
        end
    end
    retained_details = retained === nothing ? (;) : (points = retained,)
    return PolynomialChaosResult(
        formulation,
        values,
        expansion_values,
        stats,
        validation_values,
        retained_details
    )
end
