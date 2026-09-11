using Measurements
(;
    frequencies = 10.0 .^ range(-1, 7; length = 101)) -> begin
    model=load_case(:two_insulated_wires;
        variation = compose_variations(ExactOverrides(; frequencies),
            RelativeStandardUncertainty(2.0; tags = (:geometry, :cable_layer))))
    problem=ParametricProblem(model.problem)
    inner=uq_inner_formulation()
    reference=BenchmarkCalculation(:linear, problem, LinearError(inner))
    candidate=BenchmarkCalculation(:sampled,
        problem,
        MonteCarlo(inner; trials = 8, seed = UInt64(0x1234),
            return_samples = false, return_histograms = false))
    benchmark_definition(:compare_uncertainty, model.id, :manual, @__FILE__,
        model, reference, candidate, (quantities=(:R, :L, :C, :G), statistics=(:mean, :std)), (;))
end
