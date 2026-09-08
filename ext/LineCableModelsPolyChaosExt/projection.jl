function _spectral_coefficients(samples, polynomial_values, weights, norms)
    size(samples, 1) == length(weights) == size(polynomial_values, 2) || throw(
        DimensionMismatch("spectral inputs must share the collocation dimension"),
    )
    coefficients = polynomial_values * (samples .* reshape(weights, :, 1))
    coefficients ./= reshape(norms, :, 1)
    return collect(transpose(coefficients))
end

function _spectral_moments(coefficients, norms, constant_index::Int)
    means = copy(@view coefficients[:, constant_index])
    variances = zeros(eltype(coefficients), size(coefficients, 1))
    for term in axes(coefficients, 2)
        term == constant_index && continue
        variances .+= abs2.(@view coefficients[:, term]) .* norms[term]
    end
    tolerance = 256eps(Float64) * max(maximum(abs2, means; init = 0.0), 1.0)
    minimum(variances; init = 0.0) >= -tolerance || throw(ArgumentError(
        "polynomial-chaos variance is negative beyond floating-point roundoff",
    ))
    return means, sqrt.(max.(variances, zero(eltype(variances))))
end

function _relative_errors(predicted, actual)
    residual = predicted .- actual
    relative_rms = norm(residual) / max(norm(actual), eps(Float64))
    relative_max = maximum(abs, residual) /
                   max(maximum(abs, actual), eps(Float64))
    return relative_rms, relative_max, residual
end

function _validation_diagnostics(
        predicted::NamedTuple,
        actual::NamedTuple,
        shape::Tuple,
        frequency_at,
        formulation::PolynomialChaos,
        point_index::Int,
        dimension::Int,
        terms::Int,
        collocation_count::Int
)
    relative_rms_values = map(_relative_errors, predicted, actual)
    relative_rms = map(first, relative_rms_values)
    relative_max = map(value -> value[2], relative_rms_values)
    worst_error, worst_quantity_index = findmax(collect(Base.values(relative_max)))
    quantity = keys(relative_max)[worst_quantity_index]
    residual = Base.values(relative_rms_values)[worst_quantity_index][3]
    location = argmax(abs.(residual))
    native_linear_index = location[2]
    native_index = Tuple(CartesianIndices(shape)[native_linear_index])
    frequency = frequency_at(native_index)
    worst = (
        quantity,
        index = native_index,
        frequency,
        relative_rms_error = Base.values(relative_rms)[worst_quantity_index],
        relative_max_error = worst_error,
    )
    record = (
        stochastic_dimension = dimension,
        polynomial_degree = formulation.degree,
        basis_terms = terms,
        collocation_evaluations = collocation_count,
        validation_evaluations = formulation.validation_points,
        relative_rms_error = relative_rms,
        relative_max_error = relative_max,
        worst,
    )
    accepted = all(error -> error <= formulation.validation_rtol, Base.values(relative_rms)) &&
               all(error -> error <= formulation.validation_rtol, Base.values(relative_max))
    accepted || throw(ArgumentError(
        "PolynomialChaos validation failed at outer point $point_index: " *
        "distribution=:$(formulation.distribution), degree=$(formulation.degree), " *
        "quadrature_order=$(formulation.quadrature_order), d=$dimension, M=$terms, " *
        "Q=$collocation_count, V=$(formulation.validation_points), " *
        "worst_quantity=:$quantity, worst_index=$native_index, frequency=$frequency, " *
        "relative_rms_error=$(worst.relative_rms_error), " *
        "relative_max_error=$(worst.relative_max_error), " *
        "requested_tolerance=$(formulation.validation_rtol)",
    ))
    return record
end

function fit_expansion(
        multivariate_basis,
        polynomial_values,
        product_weights,
        collocation_results::AbstractVector{<:CableConstants},
        validation_polynomial_values,
        validation_results::AbstractVector{<:CableConstants},
        formulation::PolynomialChaos,
        point_index::Int,
        dimension::Int
)
    first_result = first(collocation_results)
    all(result -> typeof(result) === typeof(first_result),
        Iterators.flatten((collocation_results, validation_results))) || throw(
        ArgumentError("polynomial-chaos evaluations produced incompatible result types"),
    )
    all(result -> result.cores == first_result.cores &&
                  result.frequency == first_result.frequency,
        Iterators.flatten((collocation_results, validation_results))) || throw(
        DimensionMismatch(
            "polynomial-chaos cable-constant solves produced incompatible metadata",
        ),
    )

    quadrature_count = length(collocation_results)
    validation_count = length(validation_results)
    assembly_count = length(first_result)
    collocation_samples = map((R, L, C, G)) do selector
        samples = Matrix{Float64}(undef, quadrature_count, assembly_count)
        for sample in eachindex(collocation_results)
            observed = observe(collocation_results[sample], selector)
            all(value -> value isa Real && isfinite(value), observed) || throw(
                ArgumentError("polynomial-chaos cable-constant outputs must be finite reals"),
            )
            samples[sample, :] .= observed
        end
        samples
    end
    validation_samples = map((R, L, C, G)) do selector
        samples = Matrix{Float64}(undef, validation_count, assembly_count)
        for sample in eachindex(validation_results)
            observed = observe(validation_results[sample], selector)
            all(value -> value isa Real && isfinite(value), observed) || throw(
                ArgumentError("polynomial-chaos validation outputs must be finite reals"),
            )
            samples[sample, :] .= observed
        end
        samples
    end
    collocation_product = NamedTuple{(:R, :L, :C, :G)}(collocation_samples)
    actual_product = NamedTuple{(:R, :L, :C, :G)}(validation_samples)

    norms = vec(sum(abs2.(polynomial_values) .* transpose(product_weights); dims = 2))
    all(value -> isfinite(value) && value > 0, norms) || throw(ArgumentError(
        "polynomial-chaos basis norms must be positive and finite",
    ))
    multi_indices = PolyChaos.calculateMultiIndices(dimension, formulation.degree)
    size(multi_indices) == (size(polynomial_values, 1), dimension) || throw(
        DimensionMismatch("PolyChaos basis terms do not match the total-degree index set"),
    )
    constant_terms = findall(row -> all(iszero, row), eachrow(multi_indices))
    length(constant_terms) == 1 || throw(ArgumentError(
        "PolyChaos basis must contain exactly one constant term",
    ))
    constant_index = only(constant_terms)

    coefficient_values = map(Base.values(collocation_product)) do samples
        _spectral_coefficients(samples, polynomial_values, product_weights, norms)
    end
    coefficients = NamedTuple{(:R, :L, :C, :G)}(coefficient_values)
    moment_values = map(Base.values(coefficients)) do coefficient
        mean_values, standard_deviations =
            _spectral_moments(coefficient, norms, constant_index)
        (mean = mean_values, std = standard_deviations)
    end
    stats = NamedTuple{(:R, :L, :C, :G)}(moment_values)
    representation = CableConstants(
        first_result.cores,
        stats.R.mean,
        stats.L.mean,
        stats.C.mean,
        stats.G.mean,
        first_result.frequency
    )
    predicted_values = map(Base.values(coefficients)) do coefficient
        collect(transpose(coefficient * validation_polynomial_values))
    end
    predicted = NamedTuple{(:R, :L, :C, :G)}(predicted_values)
    record = _validation_diagnostics(
        predicted,
        actual_product,
        (assembly_count,),
        _ -> first_result.frequency,
        formulation,
        point_index,
        dimension,
        length(norms),
        quadrature_count
    )
    return (
        representation,
        expansion = (basis = multivariate_basis, coefficients),
        statistics = stats,
        validation = record,
    )
end

function fit_expansion(
        multivariate_basis,
        polynomial_values,
        product_weights,
        collocation_results::AbstractVector{<:LineParameters},
        validation_polynomial_values,
        validation_results::AbstractVector{<:LineParameters},
        formulation::PolynomialChaos,
        point_index::Int,
        dimension::Int
)
    first_result = first(collocation_results)
    all(result -> typeof(result) === typeof(first_result),
        Iterators.flatten((collocation_results, validation_results))) || throw(
        ArgumentError("polynomial-chaos evaluations produced incompatible result types"),
    )
    expected_shape = size(observe(first_result, R))
    expected_frequencies = observe(first_result, LineCableModels.frequencies)
    any(iszero, expected_frequencies) && throw(DomainError(
        expected_frequencies,
        "polynomial-chaos R/L/G/C fitting requires nonzero frequencies",
    ))
    all(result -> basis(result) === basis(first_result) &&
                  result.domain == first_result.domain &&
                  size(observe(result, R)) == expected_shape &&
                  observe(result, LineCableModels.frequencies) == expected_frequencies,
        Iterators.flatten((collocation_results, validation_results))) || throw(
        DimensionMismatch(
            "polynomial-chaos line-parameter solves produced incompatible metadata",
        ),
    )

    quadrature_count = length(collocation_results)
    validation_count = length(validation_results)
    output_count = prod(expected_shape)
    collocation_samples = map((R, L, C, G)) do selector
        samples = Matrix{Float64}(undef, quadrature_count, output_count)
        for sample in eachindex(collocation_results)
            observed = observe(collocation_results[sample], selector)
            all(value -> value isa Real && isfinite(value), observed) || throw(
                ArgumentError("polynomial-chaos line-parameter outputs must be finite reals"),
            )
            samples[sample, :] .= vec(observed)
        end
        samples
    end
    validation_samples = map((R, L, C, G)) do selector
        samples = Matrix{Float64}(undef, validation_count, output_count)
        for sample in eachindex(validation_results)
            observed = observe(validation_results[sample], selector)
            all(value -> value isa Real && isfinite(value), observed) || throw(
                ArgumentError("polynomial-chaos validation outputs must be finite reals"),
            )
            samples[sample, :] .= vec(observed)
        end
        samples
    end
    collocation_product = NamedTuple{(:R, :L, :C, :G)}(collocation_samples)
    actual_product = NamedTuple{(:R, :L, :C, :G)}(validation_samples)

    norms = vec(sum(abs2.(polynomial_values) .* transpose(product_weights); dims = 2))
    all(value -> isfinite(value) && value > 0, norms) || throw(ArgumentError(
        "polynomial-chaos basis norms must be positive and finite",
    ))
    multi_indices = PolyChaos.calculateMultiIndices(dimension, formulation.degree)
    size(multi_indices) == (size(polynomial_values, 1), dimension) || throw(
        DimensionMismatch("PolyChaos basis terms do not match the total-degree index set"),
    )
    constant_terms = findall(row -> all(iszero, row), eachrow(multi_indices))
    length(constant_terms) == 1 || throw(ArgumentError(
        "PolyChaos basis must contain exactly one constant term",
    ))
    constant_index = only(constant_terms)

    coefficient_matrices = map(Base.values(collocation_product)) do samples
        _spectral_coefficients(samples, polynomial_values, product_weights, norms)
    end
    coefficient_values = map(coefficient_matrices) do coefficient
        reshape(coefficient, expected_shape..., length(norms))
    end
    coefficients = NamedTuple{(:R, :L, :C, :G)}(coefficient_values)
    moment_values = map(coefficient_matrices) do coefficient
        mean_values, standard_deviations =
            _spectral_moments(coefficient, norms, constant_index)
        (
            mean = reshape(mean_values, expected_shape),
            std = reshape(standard_deviations, expected_shape),
        )
    end
    stats = NamedTuple{(:R, :L, :C, :G)}(moment_values)
    angular_frequency = reshape(2π .* expected_frequencies, 1, 1, :)
    representation = LineParameters(
        first_result.domain,
        complex.(stats.R.mean, stats.L.mean .* angular_frequency),
        complex.(stats.G.mean, stats.C.mean .* angular_frequency),
        expected_frequencies;
        basis = basis(first_result)
    )
    predicted_values = map(coefficient_matrices) do coefficient
        collect(transpose(coefficient * validation_polynomial_values))
    end
    predicted = NamedTuple{(:R, :L, :C, :G)}(predicted_values)
    record = _validation_diagnostics(
        predicted,
        actual_product,
        expected_shape,
        index -> expected_frequencies[index[3]],
        formulation,
        point_index,
        dimension,
        length(norms),
        quadrature_count
    )
    return (
        representation,
        expansion = (basis = multivariate_basis, coefficients),
        statistics = stats,
        validation = record,
    )
end
