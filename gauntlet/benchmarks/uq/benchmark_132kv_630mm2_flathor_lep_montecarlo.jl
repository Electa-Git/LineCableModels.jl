using Measurements
# Ordinary benchmark declaration. Construction performs no solve.
(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) -> begin
    frequencies === nothing ||
        (variation=compose_variations(variation, ExactOverrides(; frequencies)))
    model=load_case(
        :cable_132kv_630mm2_flathor;
        variation = compose_variations(
            variation, RelativeStandardUncertainty(
                10.0; tags = (:geometry, :cable_layer)))
    )
    inner=uq_inner_formulation()
    problem=ParametricProblem(model.problem)
    reference=BenchmarkCalculation(
        :linear_error,
        problem,
        LinearError(inner)
    )
    # A nested fixed-seed pilot at 128, 256, and 512 trials reduced the worst
    # meaningful std discrepancy from 9.83% to 7.90% to 3.58%. At 512 trials
    # the worst meaningful mean discrepancy was 1.23%, so 512 is the locked
    # count beneath the shared 5% mean and 10% std engineering gates.
    monte_carlo_trials=UQ_MONTE_CARLO_TRIALS
    candidate=BenchmarkCalculation(
        :monte_carlo,
        problem,
        MonteCarlo(
            inner;
            trials = monte_carlo_trials,
            seed = 0x132630,
            distribution = :normal,
            return_samples = false,
            return_histograms = false
        ); options = candidate_options
    )
    tolerances=uq_moment_tolerances()
    return benchmark_definition(
        :benchmark_132kv_630mm2_flathor_lep_montecarlo,
        :cable_132kv_630mm2_flathor,
        :uq,
        @__FILE__,
        model,
        reference,
        candidate,
        (quantities=(:R, :L, :C, :G), statistics=(:mean, :std)),
        tolerances
    )
end
