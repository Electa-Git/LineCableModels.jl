function _blas_identity()
    configuration = LinearAlgebra.BLAS.get_config()
    return sprint(show, configuration)
end

function _ratio(value, baseline, key)
    baseline === nothing ? nothing :
    value / Float64(baseline[key])
end

function _performance_baseline(id::AbstractString)
    path = joinpath(_GAUNTLET_ROOT, "config", "performance-baseline.toml")
    isfile(path) || return nothing
    baselines = TOML.parsefile(path)
    return get(get(baselines, "cases", Dict{String, Any}()), String(id), nothing)
end

"""Return a pinned performance record, or `nothing` when no baseline exists."""
function performance_baseline(id::AbstractString)
    baseline = _performance_baseline(id)
    baseline === nothing && return nothing
    median_seconds = Float64(baseline["median_seconds"])
    frequencies_value = Int(baseline["frequencies"])
    return Perf(
        Float64(baseline["minimum_seconds"]),
        median_seconds,
        Int(baseline["bytes"]),
        Int(baseline["allocations"]),
        frequencies_value,
        Int(baseline["conductors"]),
        frequencies_value / median_seconds,
        String(baseline["julia"]),
        String(baseline["os"]),
        String(baseline["arch"]),
        Int(baseline["threads"]),
        String(baseline["blas"]),
        nothing,
        nothing,
        nothing,
        nothing
    )
end

"""Benchmark one materialized case after two explicit warm-up evaluations."""
function benchmark(
        case::Case;
        samples::Int = 10,
        seconds::Real = 10,
        baseline = _performance_baseline(case.id)
)
    case.problem === nothing && throw(ArgumentError("rejected cases cannot be benchmarked"))
    compute!(case.problem, case.formulation)
    compute!(case.problem, case.formulation)
    problem = case.problem
    formulation = case.formulation
    bench = BenchmarkTools.@benchmarkable compute!($problem, $formulation) samples=samples seconds=seconds evals=1
    trial = BenchmarkTools.run(bench)
    times = Float64.(trial.times) .* 1e-9
    result = compute!(case.problem, case.formulation)
    conductors = size(Z(result), 1)
    count = length(frequencies(result))
    median_seconds = median(times)
    return Perf(
        minimum(times),
        median_seconds,
        trial.memory,
        trial.allocs,
        count,
        conductors,
        count / median_seconds,
        string(VERSION),
        string(Sys.KERNEL),
        string(Sys.ARCH),
        Threads.nthreads(),
        _blas_identity(),
        _ratio(minimum(times), baseline, "minimum_seconds"),
        _ratio(median_seconds, baseline, "median_seconds"),
        _ratio(trial.memory, baseline, "bytes"),
        _ratio(trial.allocs, baseline, "allocations")
    )
end

function _performance_comparison(performance::Perf)
    ratios = filter(!isnothing, (performance.bytes_ratio, performance.allocations_ratio))
    verdict = isempty(ratios) ? Diagnostic() :
              maximum(ratios) <= 1.05 ? Pass() : Fail()
    detail = isempty(ratios) ?
             "no pinned allocation baseline is available" :
             @sprintf("bytes ratio %.4f; allocations ratio %.4f; wall-clock ratio is informational",
        performance.bytes_ratio,
        performance.allocations_ratio)
    return Comparison(
        PerformanceCheck(),
        verdict,
        Tolerance(0.0, 0.0),
        nothing,
        detail
    )
end
