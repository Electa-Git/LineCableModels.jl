function _preflight(problem::ParametricProblem, formulation::PolynomialChaos)
    problem.space isa LineCableModels.Gridspace || throw(ArgumentError(
        "PolynomialChaos requires a parameterized Gridspace; use Combinatorial for a deterministic problem",
    ))
    outer_points = collect(points(problem.space))
    declared_count = length(problem.space)
    length(outer_points) == declared_count || throw(DimensionMismatch(
        "problem-space iteration did not match its declared cardinality",
    ))
    isempty(outer_points) && throw(ArgumentError(
        "higher-order problem space must contain at least one core problem",
    ))

    plans = map(enumerate(outer_points)) do (point_index, point)
        descriptors = uncertainties(point)
        for descriptor in descriptors
            mean_value = nominal(descriptor)
            standard_uncertainty = uncertainty(descriptor)
            mean_value isa Real && isfinite(mean_value) || throw(ArgumentError(
                "PolynomialChaos point $point_index has a non-finite or non-real nominal value",
            ))
            standard_uncertainty isa Real && isfinite(standard_uncertainty) &&
                standard_uncertainty >= 0 || throw(ArgumentError(
                "PolynomialChaos point $point_index has an invalid standard uncertainty",
            ))
        end
        active = findall(descriptor -> !iszero(uncertainty(descriptor)), descriptors)
        dimension = length(active)
        dimension > 0 || throw(ArgumentError(
            "PolynomialChaos point $point_index has zero active stochastic dimensions; use Combinatorial for a deterministic space",
        ))
        terms = binomial(big(dimension + formulation.degree), formulation.degree)
        collocation = big(formulation.quadrature_order)^dimension
        return (
            point,
            descriptors,
            active,
            dimension,
            terms,
            collocation,
            validation = big(formulation.validation_points),
        )
    end

    requested = sum(
        plan.collocation + plan.validation for plan in plans;
        init = big(0)
    )
    requested <= formulation.max_evaluations || begin
        point_costs = join(map(enumerate(plans)) do (index, plan)
            "point $index: d=$(plan.dimension), M=$(plan.terms), " *
            "Q=$(plan.collocation), V=$(plan.validation)"
        end, "; ")
        throw(ArgumentError(
            "PolynomialChaos evaluation budget exceeded for $(length(plans)) outer points " *
            "($point_costs); requested total=$requested, budget=$(formulation.max_evaluations)",
        ))
    end

    return map(plans) do plan
        merge(plan, (
            terms = Int(plan.terms),
            collocation = Int(plan.collocation),
            validation = Int(plan.validation),
        ))
    end
end

function _coordinate_bases(formulation::PolynomialChaos, dimension::Int)
    constructor = formulation.distribution === :normal ?
                  PolyChaos.GaussOrthoPoly : PolyChaos.Uniform_11OrthoPoly
    univariate = [
        constructor(
            formulation.degree;
            Nrec = formulation.quadrature_order + 1
        ) for _ in 1:dimension
    ]
    return univariate, PolyChaos.MultiOrthoPoly(univariate, formulation.degree)
end

function _tensor_rule(multivariate_basis, expected_count::Int, dimension::Int)
    node_vectors, weight_vectors = PolyChaos.nw(multivariate_basis)
    length(node_vectors) == length(weight_vectors) == dimension || throw(
        DimensionMismatch("PolyChaos returned inconsistent coordinate quadrature data"),
    )
    all(length(nodes) == length(weights)
    for (nodes, weights) in zip(node_vectors, weight_vectors)) || throw(
        DimensionMismatch("PolyChaos returned misaligned quadrature nodes and weights"),
    )

    latent = Matrix{Float64}(undef, expected_count, dimension)
    product_weights = Vector{Float64}(undef, expected_count)
    coordinate_indices = map(eachindex, node_vectors)
    tensor_indices = Iterators.product(coordinate_indices...)
    row_count = 0
    for (row, indices) in enumerate(tensor_indices)
        row_count = row
        weight = 1.0
        for coordinate in 1:dimension
            index = indices[coordinate]
            latent[row, coordinate] = node_vectors[coordinate][index]
            weight *= weight_vectors[coordinate][index]
        end
        product_weights[row] = weight
    end
    row_count == expected_count || throw(DimensionMismatch(
        "tensor quadrature produced $row_count nodes; expected $expected_count",
    ))
    size(latent) == (expected_count, dimension) || throw(DimensionMismatch(
        "tensor quadrature node matrix has an invalid shape",
    ))
    all(isfinite, latent) || throw(ArgumentError(
        "PolyChaos quadrature nodes must be finite",
    ))
    all(weight -> isfinite(weight) && weight > 0, product_weights) || throw(
        ArgumentError("PolyChaos product quadrature weights must be positive and finite"),
    )
    isapprox(sum(product_weights), 1.0; atol = 256eps(Float64), rtol = 256eps(Float64)) ||
        throw(ArgumentError(
            "PolyChaos product quadrature weights must sum to one",
        ))
    return latent, product_weights
end

function _affine_maps(descriptors, active, univariate, distribution::Symbol)
    return map(zip(active, univariate)) do (descriptor_index, coordinate_basis)
        descriptor = descriptors[descriptor_index]
        coefficients = if distribution === :normal
            PolyChaos.convert2affinePCE(
                nominal(descriptor), uncertainty(descriptor), coordinate_basis)
        else
            PolyChaos.convert2affinePCE(
                nominal(descriptor), uncertainty(descriptor), coordinate_basis;
                kind = "μσ"
            )
        end
        length(coefficients) == 2 || throw(DimensionMismatch(
            "PolyChaos affine maps must contain constant and linear coefficients",
        ))
        return (coefficients[1], coefficients[2])
    end
end

function _physical_values(descriptors, affine_maps, latent, row::Int)
    values = Vector{Any}(undef, length(descriptors))
    active_coordinate = 0
    for descriptor_index in eachindex(descriptors)
        descriptor = descriptors[descriptor_index]
        if iszero(uncertainty(descriptor))
            values[descriptor_index] = nominal(descriptor)
        else
            active_coordinate += 1
            constant, linear = affine_maps[active_coordinate]
            values[descriptor_index] =
                constant + linear * latent[row, active_coordinate]
        end
    end
    active_coordinate == length(affine_maps) || throw(DimensionMismatch(
        "affine coordinate maps did not consume every stochastic dimension",
    ))
    return Tuple(values)
end
