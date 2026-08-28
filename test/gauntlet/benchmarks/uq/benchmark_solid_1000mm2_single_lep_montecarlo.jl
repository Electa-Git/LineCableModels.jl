@testitem "UQ benchmark / solid 1000 mm² single phase / LEP versus Monte Carlo" tags=[
    :gauntlet, :uq] setup=[GauntletSupport] begin
    using Test
    using Measurements
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport

    model=load_case(
        :solid_1000mm2_single;
        variation = RelativeStandardUncertainty(
            10.0; tags = (:geometry, :cable_layer)
        )
    )
    inner=uq_inner_formulation()
    problem=ParametricProblem(model.problem)
    reference=benchmark_calculation(
        :linear_error, :uq, problem, LinearError(inner)
    )
    monte_carlo_trials=UQ_MONTE_CARLO_TRIALS
    monte_carlo_seed=UInt64(0x501000)
    candidate=benchmark_calculation(
        :monte_carlo,
        :uq,
        problem,
        MonteCarlo(
            inner;
            trials = monte_carlo_trials,
            seed = monte_carlo_seed,
            distribution = :normal,
            return_samples = false,
            return_histograms = false
        )
    )
    tolerances=uq_moment_tolerances()
    benchmark=benchmark_definition(
        :benchmark_solid_1000mm2_single_lep_montecarlo,
        :solid_1000mm2_single,
        :uq,
        @__FILE__,
        model,
        reference,
        candidate,
        UQMomentPolicy(),
        tolerances
    )
    outcome=run_benchmark(benchmark)
    @info "LEP versus Monte Carlo engineering comparison" trials=monte_carlo_trials seed=monte_carlo_seed errors=moment_error_summary(
        outcome.comparison, tolerances.reference) monte_carlo_over_lep_speedup=outcome.performance.speedup

    @test outcome.reference_result isa LinearErrorResult
    @test outcome.candidate_result isa MonteCarloResult
    @test reference.problem === candidate.problem
    @test reference.formulation.inner === candidate.formulation.inner
    @test only(outcome.candidate_result.trial_counts) == monte_carlo_trials
    @test outcome.candidate_result.root_seed == monte_carlo_seed
    @test size(outcome.reference.values.R.mean) == model.expected_size
    @test size(outcome.candidate.values.R.mean) == model.expected_size
    @test outcome.passes
    @test outcome.performance !== nothing
    @test outcome.performance.passes !== false
end
