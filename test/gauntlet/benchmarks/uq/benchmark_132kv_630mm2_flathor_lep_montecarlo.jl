@testitem "UQ benchmark / 132 kV 630 mm² / LEP versus Monte Carlo" tags=[:gauntlet, :uq] setup=[GauntletSupport] begin
    using Test
    using Measurements
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport

    model = load_case(
        :cable_132kv_630mm2_flathor;
        variation = RelativeStandardUncertainty(
            10.0; tags = (:geometry, :cable_layer)
        )
    )
    inner = uq_inner_formulation()
    problem = ParametricProblem(model.problem)
    reference = benchmark_calculation(
        :linear_error,
        :uq,
        problem,
        LinearError(inner)
    )
    # A nested fixed-seed pilot at 128, 256, and 512 trials reduced the worst
    # meaningful std discrepancy from 9.83% to 7.90% to 3.58%. At 512 trials
    # the worst meaningful mean discrepancy was 1.23%, so 512 is the locked
    # count beneath the shared 5% mean and 10% std engineering gates.
    monte_carlo_trials = UQ_MONTE_CARLO_TRIALS
    candidate = benchmark_calculation(
        :monte_carlo,
        :uq,
        problem,
        MonteCarlo(
            inner;
            trials = monte_carlo_trials,
            seed = 0x132630,
            distribution = :normal,
            return_samples = false,
            return_histograms = false
        )
    )
    tolerances = uq_moment_tolerances()
    benchmark = benchmark_definition(
        :benchmark_132kv_630mm2_flathor_lep_montecarlo,
        :cable_132kv_630mm2_flathor,
        :uq,
        @__FILE__,
        model,
        reference,
        candidate,
        UQMomentPolicy(),
        tolerances
    )
    outcome = run_benchmark(benchmark)
    @info "LEP versus Monte Carlo engineering comparison" trials=monte_carlo_trials seed=UInt64(0x132630) errors=moment_error_summary(outcome.comparison, tolerances.reference) monte_carlo_over_lep_speedup=outcome.performance.speedup

    @test outcome.reference_result isa LinearErrorResult
    @test outcome.candidate_result isa MonteCarloResult
    @test reference.problem === candidate.problem
    @test reference.formulation.inner === candidate.formulation.inner
    @test reference.options == candidate.options == (;)
    @test only(outcome.candidate_result.trial_counts) == monte_carlo_trials
    @test outcome.candidate_result.root_seed == UInt64(0x132630)
    @test size(outcome.reference.values.R.mean) == model.expected_size
    @test size(outcome.candidate.values.R.mean) == model.expected_size
    @test outcome.passes
    @test outcome.performance !== nothing
    @test outcome.performance.passes !== false
end
