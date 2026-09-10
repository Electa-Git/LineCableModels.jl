using Measurements
# Ordinary benchmark declaration. Construction performs no solve.
(; frequencies = nothing,
    reference_options = (;),
    candidate_options = (;),
    variation = NoVariation()) -> begin
    frequencies === nothing ||
        (variation=compose_variations(variation, ExactOverrides(; frequencies)))
    model=load_case(
        :cable_640kv_2000mm2_bipole;
        variation = compose_variations(
            variation, RelativeStandardUncertainty(
                10.0; tags = (:geometry, :cable_layer)))
    )
    inner=uq_inner_formulation()
    problem=ParametricProblem(model.problem)
    reference=BenchmarkCalculation(
        :linear_error, problem, LinearError(inner); options = reference_options
    )
    monte_carlo_trials=UQ_MONTE_CARLO_TRIALS
    monte_carlo_seed=UInt64(0x640200)
    candidate=BenchmarkCalculation(
        :monte_carlo,
        problem,
        MonteCarlo(
            inner;
            trials = monte_carlo_trials,
            seed = monte_carlo_seed,
            distribution = :normal,
            return_samples = false,
            return_histograms = false
        ); options = candidate_options
    )
    tolerances=uq_moment_tolerances()
    return benchmark_definition(
        :benchmark_640kv_2000mm2_bipole_lep_montecarlo,
        :cable_640kv_2000mm2_bipole,
        :uq,
        @__FILE__,
        model,
        reference,
        candidate,
        (quantities=(:R, :L, :C, :G), statistics=(:mean, :std)),
        tolerances
    )
end
