@testitem "UQ / polynomial chaos / dependency-free contracts" tags = [:unit] begin
    using Statistics

    inner = Formulation()
    formulation = PolynomialChaos(inner)
    @test formulation.inner === inner
    @test formulation.degree == 3
    @test formulation.quadrature_order == 4
    @test formulation.distribution === :normal
    @test formulation.validation_points == 16
    @test formulation.validation_rtol == 1.0e-3
    @test formulation.max_evaluations == 10_000
    @test formulation.validation_seed == UInt64(0)
    @test formulation.options == (retain_details = false,)
    @test computation_options(PolynomialChaos, (;)) == (retain_details = false,)

    configured = PolynomialChaos(
        inner;
        degree = 2,
        quadrature_order = 5,
        distribution = :uniform,
        validation_points = 3,
        validation_rtol = 1.0e-7,
        max_evaluations = 99,
        validation_seed = 0x2a,
        options = (retain_details = true,)
    )
    @test configured.degree == 2
    @test configured.quadrature_order == 5
    @test configured.distribution === :uniform
    @test configured.validation_seed == UInt64(42)
    @test configured.options == (retain_details = true,)

    for keyword in (
            (degree = 0,),
            (degree = true,),
            (degree = 2, quadrature_order = 2),
            (distribution = :lognormal,),
            (validation_points = 0,),
            (validation_rtol = 0.0,),
            (validation_rtol = Inf,),
            (max_evaluations = 0,),
            (options = (unknown = true,),))
        @test_throws ArgumentError PolynomialChaos(inner; keyword...)
    end
    @test_throws InexactError PolynomialChaos(inner; validation_seed = -1)

    value = CableConstants(2.0, 3.0, 4.0, 5.0; frequency = 50.0)
    coefficients = (
        R = reshape([2.0, 0.1], 1, 2),
        L = reshape([3.0, 0.2], 1, 2),
        C = reshape([4.0, 0.3], 1, 2),
        G = reshape([5.0, 0.4], 1, 2),
    )
    stats = (
        R = (mean = [2.0], std = [0.1]),
        L = (mean = [3.0], std = [0.2]),
        C = (mean = [4.0], std = [0.3]),
        G = (mean = [5.0], std = [0.4]),
    )
    errors = (R = 1.0e-8, L = 2.0e-8, C = 3.0e-8, G = 4.0e-8)
    record = (
        stochastic_dimension = 1,
        polynomial_degree = 1,
        basis_terms = 2,
        collocation_evaluations = 2,
        validation_evaluations = 2,
        relative_rms_error = errors,
        relative_max_error = errors,
        worst = (
            quantity = :G,
            index = (1,),
            frequency = 50.0,
            relative_rms_error = 4.0e-8,
            relative_max_error = 4.0e-8,
        ),
    )
    result_formulation = PolynomialChaos(
        inner; degree = 1, validation_points = 2, options = (retain_details = true,))
    retained = (
        points = [(
            collocation = [value, value],
            validation = [value, value],
        )],
    )
    result = PolynomialChaosResult(
        result_formulation,
        [value],
        [(basis = :contract_basis, coefficients = coefficients)],
        [stats],
        [record],
        retained
    )

    @test result isa AbstractUncertaintyResult{CableConstants{Float64}}
    @test eltype(result) === CableConstants{Float64}
    @test length(result) == 1
    @test size(result) == (1,)
    @test firstindex(result) == lastindex(result) == 1
    @test only(result) === value
    @test collect(result) == [value]
    @test only(expansions(result)).basis === :contract_basis
    @test only(expansions(result)).coefficients === coefficients
    @test only(statistics(result)) === stats
    @test only(validation(result)) === record
    @test details(result) === retained
    @test basis(result) === :pul
    @test !applicable(samples, result)
    @test !applicable(histograms, result)

    @test observables(typeof(result)) == (
        R, L, C, G,
        (statistics, R), (statistics, L), (statistics, C), (statistics, G)
    )
    @test observe(result, R, 1, 1) == 2.0
    @test observe(result, statistics, R, 1, 1) == (mean = 2.0, std = 0.1)
    @test observe(result, statistics, R, mean, 1, 1) == 2.0
    @test observe(result, statistics, R, std, 1, 1) == 0.1
    publication = observables(
        result,
        (
            (statistics, R, mean, 1, 1),
            (statistics, R, std, 1, 1),
        );
        units = (:base, :base),
        length_unit = :base
    )
    @test keys(publication.columns) == (:mean_R, :std_R)
    factor = scale_factor(R, basis(result), publication[1].unit)
    @test publication.columns.mean_R == [2.0 * factor]
    @test publication.columns.std_R == [0.1 * abs(factor)]
    @test occursin("PolynomialChaos(degree=1", sprint(show, result_formulation))
    @test occursin("PolynomialChaosResult(1 points", sprint(show, result))

    @test_throws DimensionMismatch PolynomialChaosResult(
        result_formulation, [value], [(basis = :basis, coefficients = coefficients)],
        [stats, stats], [record])
    bad_coefficients = merge(coefficients, (R = zeros(1, 3),))
    @test_throws DimensionMismatch PolynomialChaosResult(
        result_formulation, [value], [(basis = :basis, coefficients = bad_coefficients)],
        [stats], [record])
    bad_stats = merge(stats, (R = (mean = [2.0], std = [-0.1]),))
    @test_throws ArgumentError PolynomialChaosResult(
        result_formulation, [value], [(basis = :basis, coefficients = coefficients)],
        [bad_stats], [record])
    bad_record = merge(record, (validation_evaluations = 3,))
    @test_throws DimensionMismatch PolynomialChaosResult(
        result_formulation, [value], [(basis = :basis, coefficients = coefficients)],
        [stats], [bad_record])

    frequencies_axis = [50.0, 100.0]
    resistance = reshape([1.0, 2.0], 1, 1, 2)
    inductance = reshape([3.0, 4.0], 1, 1, 2)
    conductance = reshape([5.0, 6.0], 1, 1, 2)
    capacitance = reshape([7.0, 8.0], 1, 1, 2)
    angular_frequency = reshape(2π .* frequencies_axis, 1, 1, :)
    line_value = LineParameters(
        complex.(resistance, inductance .* angular_frequency),
        complex.(conductance, capacitance .* angular_frequency),
        frequencies_axis
    )
    line_coefficients = map(
        field -> cat(field, zeros(size(field)); dims = 4),
        (R = resistance, L = inductance, C = capacitance, G = conductance)
    )
    line_stats = map(
        field -> (mean = field, std = zeros(size(field))),
        (R = resistance, L = inductance, C = capacitance, G = conductance)
    )
    line_record = merge(record, (
        stochastic_dimension = 0,
        worst = merge(record.worst, (index = (1, 1, 2), frequency = 100.0,)),
    ))
    line_result = PolynomialChaosResult(
        result_formulation,
        [line_value],
        [(basis = :line_basis, coefficients = line_coefficients)],
        [line_stats],
        [line_record]
    )
    @test observe(line_result, frequencies, 1, :) == frequencies_axis
    @test observe(line_result, statistics, L, mean, 1, 1, 1, :) == vec(inductance)
    @test observables(typeof(line_result)) == (
        frequencies, R, L, C, G,
        (statistics, R), (statistics, L), (statistics, C), (statistics, G)
    )
end
