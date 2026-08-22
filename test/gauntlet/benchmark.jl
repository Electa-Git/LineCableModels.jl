function _performance_identity()
    return (
        julia_version = string(VERSION),
        kernel = string(Sys.KERNEL),
        architecture = string(Sys.ARCH),
        threads = Threads.nthreads(),
        blas = sprint(show, BLAS.get_config())
    )
end

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
        samples = length(times),
        environment = _performance_identity()
    )
end

function performance_comparison(accepted, current, tolerance)
    same_environment = accepted.environment == current.environment
    ratios = (
        median_time = current.median_seconds / accepted.median_seconds,
        bytes = accepted.bytes == 0 ? (current.bytes == 0 ? 1.0 : Inf) :
                current.bytes / accepted.bytes,
        allocations = accepted.allocations == 0 ?
                      (current.allocations == 0 ? 1.0 : Inf) :
                      current.allocations / accepted.allocations
    )
    passes = same_environment ?
             ratios.median_time <= tolerance.median_time_ratio &&
             ratios.bytes <= tolerance.bytes_ratio &&
             ratios.allocations <= tolerance.allocations_ratio : nothing
    return (; comparable = same_environment, ratios, passes)
end
