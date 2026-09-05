@testitem "Extensions / PolyChaos boundary / unloaded core" tags = [
    :extension,
    :core_only,
] begin
    import LineCableModels

    @test Base.get_extension(
        LineCableModels, :LineCableModelsPolyChaosExt) === nothing
    @test isdefined(LineCableModels, :PolynomialChaos)
    @test isdefined(LineCableModels, :PolynomialChaosResult)
    @test !isdefined(LineCableModels, :PolyChaos)

    struct CoreOnlyProblem <: LineCableModels.AbstractProblemDefinition end
    space = LineCableModels.Gridspace{CoreOnlyProblem}(CoreOnlyProblem, ())
    problem = LineCableModels.ParametricProblem(space)
    formulation = LineCableModels.PolynomialChaos(LineCableModels.Formulation())
    @test_throws MethodError LineCableModels.compute(problem, formulation)
end

@testitem "Extensions / PolyChaos boundary / activation" tags = [:extension] begin
    import LineCableModels
    @test Base.get_extension(
        LineCableModels, :LineCableModelsPolyChaosExt) === nothing

    using PolyChaos

    extension = Base.get_extension(LineCableModels, :LineCableModelsPolyChaosExt)
    @test extension !== nothing
    @test pkgversion(PolyChaos) == v"2.1.0"
    @test any(
        method -> method.module === extension,
        methods(LineCableModels.compute)
    )
    @test isempty(Test.detect_ambiguities(
        LineCableModels, extension; recursive = true))
end

@testmodule PolynomialChaosAnalyticSupport begin
    using LineCableModels
    import LineCableModels: compute, computation_details, validate

    struct AnalyticProblem <: AbstractProblemDefinition
        x::Float64
        y::Float64
        fixed::Float64
    end
    validate(problem::AnalyticProblem) = problem

    struct CableAnalytic <: AbstractFormulation
        solves::Base.RefValue{Int}
    end
    function compute(problem::AnalyticProblem, formulation::CableAnalytic; options = (;))
        formulation.solves[] += 1
        x, y, fixed = problem.x, problem.y, problem.fixed
        return CableConstants(
            x + 2y + fixed,
            x * y,
            x^2,
            3x - y;
            frequency = 50.0
        )
    end
    function computation_details(::Type{CableAnalytic}, result::CableConstants)
        return (resistance = copy(result.R),)
    end

    struct LineAnalytic <: AbstractFormulation
        solves::Base.RefValue{Int}
    end
    function compute(problem::AnalyticProblem, formulation::LineAnalytic; options = (;))
        formulation.solves[] += 1
        x, y = problem.x, problem.y
        frequencies = [50.0, 100.0]
        factors = reshape(collect(1.0:8.0), 2, 2, 2)
        resistance = factors .* (x + y)
        inductance = factors .* (x * y)
        conductance = factors .* (2x - y)
        capacitance = factors .* x^2
        angular = reshape(2π .* frequencies, 1, 1, :)
        return LineParameters(
            complex.(resistance, inductance .* angular),
            complex.(conductance, capacitance .* angular),
            frequencies;
            basis = :total
        )
    end

    function analytic_space(
            x = 2.0,
            y = 3.0;
            sigma_x = 0.2,
            sigma_y = 0.4,
            fixed = 7.0,
            sigma_fixed = 0.0
    )
        return Gridspace{AnalyticProblem}(
            AnalyticProblem,
            (
                Grid(x, AbsoluteError(sigma_x)),
                Grid(y, AbsoluteError(sigma_y)),
                Grid(fixed, AbsoluteError(sigma_fixed)),
            )
        )
    end

    struct ScalarProblem <: AbstractProblemDefinition
        x::Float64
    end
    validate(problem::ScalarProblem) = problem

    struct PositiveProblem <: AbstractProblemDefinition
        x::Float64
    end
    function validate(problem::PositiveProblem)
        problem.x > 0 || throw(DomainError(problem.x, "x must be positive"))
        return problem
    end

    struct PositiveFormulation <: AbstractFormulation
        solves::Base.RefValue{Int}
    end
    function compute(problem::PositiveProblem, formulation::PositiveFormulation;
            options = (;))
        formulation.solves[] += 1
        return CableConstants(problem.x, 1.0, 1.0, 1.0)
    end

    struct CableVariant <: AbstractFormulation
        solves::Base.RefValue{Int}
        mode::Symbol
    end
    struct UnsupportedResult
        value::Float64
    end
    function compute(problem::ScalarProblem, formulation::CableVariant; options = (;))
        formulation.solves[] += 1
        x = problem.x
        if formulation.mode === :result_type
            return x < 0 ? CableConstants(Float32(1), Float32(1), Float32(1)) :
                   CableConstants(1.0, 1.0, 1.0)
        elseif formulation.mode === :metadata
            core_name = x < 0 ? :negative : :positive
            return CableConstants(1.0, 1.0, 1.0; core = core_name)
        elseif formulation.mode === :nonfinite
            result = CableConstants(1.0, 1.0, 1.0)
            result.R[1] = Inf
            return result
        elseif formulation.mode === :unsupported
            return UnsupportedResult(x)
        elseif formulation.mode === :cubic
            return CableConstants(x^3, x^3 + x, x^3 - x, 2x^3)
        end
        error("unsupported test mode")
    end

    struct VariantDomain <: LineCableModels.Engine.LineParamsDomain
        tag::Int
    end
    LineCableModels.domain(::Type{VariantDomain}) = :variant

    struct LineVariant <: AbstractFormulation
        solves::Base.RefValue{Int}
        mode::Symbol
    end
    function compute(problem::ScalarProblem, formulation::LineVariant; options = (;))
        formulation.solves[] += 1
        x = problem.x
        count = formulation.mode === :shape && x >= 0 ? 2 : 1
        frequency = formulation.mode === :zero_frequency ? [0.0] :
                    formulation.mode === :frequency && x >= 0 ? [60.0] : [50.0]
        values = fill(complex(1.0, 1.0), count, count, 1)
        result_basis = formulation.mode === :basis && x >= 0 ? :total : :pul
        result_domain = formulation.mode === :domain ?
                        VariantDomain(x >= 0 ? 2 : 1) : VariantDomain(1)
        return LineParameters(
            result_domain,
            values,
            values,
            frequency;
            basis = result_basis
        )
    end

    function scalar_space(mean_value = 0.0, sigma = 1.0)
        return Gridspace{ScalarProblem}(
            ScalarProblem,
            (Grid(mean_value, AbsoluteError(sigma)),)
        )
    end

    function positive_space(mean_value = 0.1, sigma = 1.0)
        return Gridspace{PositiveProblem}(
            PositiveProblem,
            (Grid(mean_value, AbsoluteError(sigma)),)
        )
    end
end

@testitem "Extensions / PolyChaos / exact analytic CableConstants" tags = [
    :extension,
    :integration,
] setup = [PolynomialChaosAnalyticSupport] begin
    import PolyChaos
    using Random
    import Statistics

    means = (x = 2.0, y = 3.0, fixed = 7.0)
    sigmas = (x = 0.2, y = 0.4)
    for distribution in (:normal, :uniform)
        solves = Ref(0)
        inner = PolynomialChaosAnalyticSupport.CableAnalytic(solves)
        formulation = PolynomialChaos(
            inner;
            degree = 2,
            distribution,
            validation_points = 7,
            validation_rtol = 1.0e-11,
            validation_seed = 0x55,
            options = (retain_details = true,)
        )
        Random.seed!(104)
        expected_global_draws = rand(5)
        Random.seed!(104)
        result = compute(
            ParametricProblem(PolynomialChaosAnalyticSupport.analytic_space()),
            formulation
        )
        @test rand(5) == expected_global_draws
        @test solves[] == 3^2 + 7
        @test result isa PolynomialChaosResult{<:CableConstants}
        @test only(result).cores == [:core]
        @test only(result).frequency == 50.0

        expected_mean = (
            R = means.x + 2 * means.y + means.fixed,
            L = means.x * means.y,
            C = means.x^2 + sigmas.x^2,
            G = 3 * means.x - means.y,
        )
        expected_std = (
            R = sqrt(sigmas.x^2 + 4sigmas.y^2),
            L = sqrt(
                (means.x^2 + sigmas.x^2) * (means.y^2 + sigmas.y^2) -
                means.x^2 * means.y^2
            ),
            C = distribution === :normal ?
                sqrt(2 * sigmas.x^4 + 4 * means.x^2 * sigmas.x^2) :
                sqrt(0.8 * sigmas.x^4 + 4 * means.x^2 * sigmas.x^2),
            G = sqrt(9 * sigmas.x^2 + sigmas.y^2),
        )
        stats = only(statistics(result))
        for quantity in (:R, :L, :C, :G)
            @test only(getproperty(stats, quantity).mean) ≈
                  getproperty(expected_mean, quantity) rtol = 2.0e-14
            @test only(getproperty(stats, quantity).std) ≈
                  getproperty(expected_std, quantity) rtol = 2.0e-13
        end

        expansion = only(expansions(result))
        @test size(expansion.coefficients.R) == (1, 6)
        for quantity in (:R, :L, :C, :G)
            coefficients = vec(getproperty(expansion.coefficients, quantity))
            @test PolyChaos.mean(coefficients, expansion.basis) ≈
                  only(getproperty(stats, quantity).mean) rtol = 2.0e-14
            @test PolyChaos.std(coefficients, expansion.basis) ≈
                  only(getproperty(stats, quantity).std) rtol = 2.0e-13
        end
        latent_check = [0.2 0.3; -0.1 0.5]
        coefficient_check = vec(expansion.coefficients.R)
        matrix_prediction = collect(transpose(
            reshape(coefficient_check, 1, :) *
            PolyChaos.evaluate(latent_check, expansion.basis)
        ))
        @test vec(matrix_prediction) ≈
              PolyChaos.evaluatePCE(coefficient_check, latent_check, expansion.basis)
        diagnostic = only(validation(result))
        @test diagnostic.stochastic_dimension == 2
        @test diagnostic.basis_terms == 6
        @test diagnostic.collocation_evaluations == 9
        @test diagnostic.validation_evaluations == 7
        @test all(<=(formulation.validation_rtol), diagnostic.relative_rms_error)
        @test all(<=(formulation.validation_rtol), diagnostic.relative_max_error)
        @test keys(details(result)) == (:points,)
        @test length(only(details(result).points).collocation) == 9
        @test length(only(details(result).points).validation) == 7
        @test keys(first(only(details(result).points).collocation)) == (:resistance,)
        @test !applicable(samples, result)
        @test !applicable(histograms, result)
        @test @inferred(statistics(result)) === result.stats
        @test @inferred(expansions(result)) === result.expansion_values
        @test @inferred(validation(result)) === result.validation_values
        @test @inferred(observe(
            result, statistics, R, Statistics.mean, 1, 1)) ≈ expected_mean.R
        @test @inferred(observe(
            result, statistics, R, Statistics.std, 1, 1)) ≈ expected_std.R
    end
end

@testitem "Extensions / PolyChaos / exact analytic LineParameters and observations" tags = [
    :extension,
    :integration,
] setup = [PolynomialChaosAnalyticSupport] begin
    import PolyChaos
    import Statistics

    solves = Ref(0)
    formulation = PolynomialChaos(
        PolynomialChaosAnalyticSupport.LineAnalytic(solves);
        degree = 2,
        validation_points = 7,
        validation_rtol = 1.0e-11,
        validation_seed = 91
    )
    result = compute(
        ParametricProblem(PolynomialChaosAnalyticSupport.analytic_space()),
        formulation
    )
    @test solves[] == 16
    @test result isa PolynomialChaosResult{<:LineParameters}
    @test basis(only(result)) === :total
    @test frequencies(only(result)) == [50.0, 100.0]
    @test size(only(result).Z) == (2, 2, 2)

    factors = reshape(collect(1.0:8.0), 2, 2, 2)
    expected_mean = (R = 5.0, L = 6.0, C = 4.04, G = 1.0)
    expected_std = (
        R = sqrt(0.2^2 + 0.4^2),
        L = sqrt((2.0^2 + 0.2^2) * (3.0^2 + 0.4^2) - 6.0^2),
        C = sqrt(2 * 0.2^4 + 4 * 2.0^2 * 0.2^2),
        G = sqrt(4 * 0.2^2 + 0.4^2),
    )
    stats = only(statistics(result))
    expansion = only(expansions(result))
    for quantity in (:R, :L, :C, :G)
        @test getproperty(stats, quantity).mean ≈
              factors .* getproperty(expected_mean, quantity) rtol = 2.0e-14
        @test getproperty(stats, quantity).std ≈
              factors .* getproperty(expected_std, quantity) rtol = 2.0e-13
        @test size(getproperty(expansion.coefficients, quantity)) == (2, 2, 2, 6)
        @test all(isfinite, getproperty(expansion.coefficients, quantity))
    end
    @test all(isfinite, only(validation(result)).relative_rms_error)
    @test all(isfinite, only(validation(result)).relative_max_error)

    stored_statistics = deepcopy(stats)
    requests = (
        (statistics, R, Statistics.mean, 1, Colon(), Colon(), Colon()),
        (statistics, R, Statistics.std, 1, Colon(), Colon(), Colon()),
        (statistics, L, Statistics.mean, 1, Colon(), Colon(), Colon()),
        (statistics, L, Statistics.std, 1, Colon(), Colon(), Colon()),
        (statistics, C, Statistics.mean, 1, Colon(), Colon(), Colon()),
        (statistics, C, Statistics.std, 1, Colon(), Colon(), Colon()),
        (statistics, G, Statistics.mean, 1, Colon(), Colon(), Colon()),
        (statistics, G, Statistics.std, 1, Colon(), Colon(), Colon()),
    )
    publication = observables(result, requests; units = ntuple(_ -> :milli, 8))
    @test keys(publication.columns) == (
        :mean_R, :std_R, :mean_L, :std_L,
        :mean_C, :std_C, :mean_G, :std_G,
    )
    @test length(first(values(publication.columns))) == 8
    @test stats == stored_statistics
    @test publication[1].values != stats.R.mean
    @test observe(result, frequencies, 1, :) == [50.0, 100.0]
end

@testitem "Extensions / PolyChaos / cost and failure boundaries" tags = [
    :extension,
] setup = [PolynomialChaosAnalyticSupport] begin
    using PolyChaos

    solves = Ref(0)
    analytic = PolynomialChaosAnalyticSupport.CableAnalytic(solves)
    budgeted = PolynomialChaos(
        analytic;
        degree = 4,
        validation_points = 5,
        max_evaluations = 10
    )
    budget_error = try
        compute(
            ParametricProblem(PolynomialChaosAnalyticSupport.analytic_space()),
            budgeted
        )
        nothing
    catch exception
        exception
    end
    @test budget_error isa ArgumentError
    @test solves[] == 0
    budget_message = sprint(showerror, budget_error)
    for token in ("outer points", "d=2", "M=15", "Q=25", "V=5",
            "requested total=30", "budget=10")
        @test occursin(token, budget_message)
    end

    deterministic_solves = Ref(0)
    deterministic = PolynomialChaosAnalyticSupport.CableAnalytic(deterministic_solves)
    @test_throws ArgumentError compute(
        ParametricProblem(PolynomialChaosAnalyticSupport.analytic_space(
            sigma_x = 0.0, sigma_y = 0.0)),
        PolynomialChaos(deterministic)
    )
    @test deterministic_solves[] == 0
    @test_throws MethodError PolynomialChaos(analytic; covariance = [1.0 0.0; 0.0 1.0])

    positive_solves = Ref(0)
    positive_space = PolynomialChaosAnalyticSupport.positive_space(0.1, 1.0)
    @test_throws DomainError compute(
        ParametricProblem(positive_space),
        PolynomialChaos(
            PolynomialChaosAnalyticSupport.PositiveFormulation(positive_solves);
            degree = 1,
            validation_points = 2
        )
    )
    @test positive_solves[] == 0

    cable_failures = (
        result_type = ArgumentError,
        metadata = DimensionMismatch,
        nonfinite = ArgumentError,
    )
    for mode in keys(cable_failures)
        variant_solves = Ref(0)
        variant = PolynomialChaosAnalyticSupport.CableVariant(variant_solves, mode)
        @test_throws getproperty(cable_failures, mode) compute(
            ParametricProblem(PolynomialChaosAnalyticSupport.scalar_space()),
            PolynomialChaos(variant; degree = 1, validation_points = 2)
        )
        @test variant_solves[] > 0
    end
    unsupported = PolynomialChaosAnalyticSupport.CableVariant(Ref(0), :unsupported)
    @test_throws MethodError compute(
        ParametricProblem(PolynomialChaosAnalyticSupport.scalar_space()),
        PolynomialChaos(unsupported; degree = 1, validation_points = 2)
    )

    line_failures = (
        frequency = DimensionMismatch,
        basis = ArgumentError,
        domain = DimensionMismatch,
        shape = DimensionMismatch,
        zero_frequency = DomainError,
    )
    for mode in keys(line_failures)
        variant = PolynomialChaosAnalyticSupport.LineVariant(Ref(0), mode)
        @test_throws getproperty(line_failures, mode) compute(
            ParametricProblem(PolynomialChaosAnalyticSupport.scalar_space()),
            PolynomialChaos(variant; degree = 1, validation_points = 2)
        )
    end

    cubic = PolynomialChaosAnalyticSupport.CableVariant(Ref(0), :cubic)
    validation_error = try
        compute(
            ParametricProblem(PolynomialChaosAnalyticSupport.scalar_space(2.0, 0.5)),
            PolynomialChaos(
                cubic;
                degree = 1,
                validation_points = 4,
                validation_rtol = 1.0e-14,
                validation_seed = 7
            )
        )
        nothing
    catch exception
        exception
    end
    @test validation_error isa ArgumentError
    validation_message = sprint(showerror, validation_error)
    for token in (
            "outer point 1", "distribution=:normal", "degree=1",
            "quadrature_order=2", "d=1", "M=2", "Q=2", "V=4",
            "worst_quantity", "worst_index", "frequency",
            "relative_rms_error", "relative_max_error", "requested_tolerance")
        @test occursin(token, validation_message)
    end
end

@testitem "Extensions / PolyChaos / real CableConstants smoke" tags = [
    :extension,
    :integration,
] begin
    using PolyChaos

    conductor = Material(
        kind = :conductor,
        rho = 1.7241e-8,
        mu_r = Grid(0.999994, AbsoluteError(1.0e-4)),
        alpha = 0.00393
    )
    dielectric = Material(
        kind = :insulator,
        rho = 1.0e14,
        eps_r = Grid(2.3, AbsoluteError(0.01)),
        tan_delta = Grid(0.01, AbsoluteError(1.0e-4))
    )
    design = build(
        CableDesign,
        "pce-smoke",
        terminal(:core, core(conductor; r = 0.01)),
        insulation(dielectric; t = 0.005)
    )
    problem_space = CableConstantsProblem(design; frequency = 50.0)
    point = only(LineCableModels.points(problem_space))
    @test length(LineCableModels.uncertainties(point)) == 3
    nominal_values = map(nominal, LineCableModels.uncertainties(point))
    nominal_problem = LineCableModels.realize(
        point,
        LineCableModels.realize_arguments(point, nominal_values)
    )
    deterministic_before = compute(nominal_problem, CableConstantsFormulation())

    result = compute(
        ParametricProblem(problem_space),
        PolynomialChaos(
            CableConstantsFormulation();
            degree = 1,
            validation_points = 4,
            validation_rtol = 1.0e-6,
            validation_seed = 19
        )
    )
    monte_carlo = compute(
        ParametricProblem(problem_space),
        MonteCarlo(CableConstantsFormulation(); trials = 16, seed = 19)
    )
    deterministic_after = compute(nominal_problem, CableConstantsFormulation())

    @test deterministic_after == deterministic_before
    @test only(result) isa CableConstants
    @test only(result).cores == deterministic_before.cores
    @test only(result).frequency == deterministic_before.frequency
    @test size(only(expansions(result)).coefficients.R) == (1, 4)
    @test size(only(statistics(result)).R.mean) == (1,)
    diagnostic = only(validation(result))
    @test diagnostic.stochastic_dimension == 3
    @test diagnostic.basis_terms == 4
    @test diagnostic.collocation_evaluations == 8
    @test diagnostic.validation_evaluations == 4
    @test all(isfinite, diagnostic.relative_rms_error)
    @test all(isfinite, diagnostic.relative_max_error)
    @test diagnostic.collocation_evaluations + diagnostic.validation_evaluations <
          trial_count(monte_carlo, 1)
end
