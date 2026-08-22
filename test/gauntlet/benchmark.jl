function benchmark_local(case; samples::Int = 10, seconds::Real = 10)
    compute!(case.problem, case.formulation)
    compute!(case.problem, case.formulation)
    problem = case.problem
    formulation = case.formulation
    benchmark = BenchmarkTools.@benchmarkable compute!($problem, $formulation) samples=samples seconds=seconds evals=1
    trial = BenchmarkTools.run(benchmark)
    times = Float64.(trial.times) .* 1.0e-9
    return (
        minimum_seconds = minimum(times),
        median_seconds = median(times),
        bytes = trial.memory,
        allocations = trial.allocs,
        samples = length(times)
    )
end
