"""
$(TYPEDSIGNATURES)

Repeat a caller-supplied computation and collect its retained scan measurements.
The callable must request `timing=true` through ordinary computation options.
No additional stopwatch surrounds the callable.

# Arguments

- `f`: Zero-argument callable returning line parameters, an ordinary batch,
  a parametric or linear-error result, or a Monte Carlo result.

# Keywords

- `samples=1`: Positive number of measured repetitions; Boolean values are rejected.
- `warmup=0`: Nonnegative number of discarded repetitions; Boolean values are rejected.

# Returns

- `(; result, timings)`: The last scientific result and detached timing payloads
  for each measured repetition. Batches and parametric/linear-error results keep
  their value order; Monte Carlo measurements keep population/trial order.
  Reused results retain empty timing records.

# Notes

Owned scans retain wall, GC, compilation, and recompilation durations in seconds,
and Julia allocation volume in bytes. Native scans retain caller wall time and
their backend measurements in seconds. These describe complete frequency scans,
excluding batch-shared preparation, measurement attachment, and callbacks.
The callable controls verbosity, seeds, callbacks, and reuse. No options are changed.

# Errors

- `ArgumentError`: Invalid counts, unsupported results, or absent requested timing.
- Exceptions from `f` propagate immediately without retry.

# Examples

```julia
measurement = LineCableModels.benchmark(; samples=3, warmup=1) do
    compute(problem, formulation; options=(timing=true, verbosity=(default=0,)))
end
```
"""
function benchmark(f; samples = 1, warmup = 0)
    samples isa Integer && !(samples isa Bool) && 0 < samples <= typemax(Int) ||
        throw(ArgumentError("samples must be a positive integer representable as Int"))
    warmup isa Integer && !(warmup isa Bool) && 0 <= warmup <= typemax(Int) ||
        throw(ArgumentError("warmup must be a nonnegative integer representable as Int"))
    for _ in 1:warmup
        f()
    end
    # Projection consumes completed outputs only; it never visits input spaces.
    function measurement(value)
        if value isa Union{LineParameters, MonteCarloResult}
            haskey(details(value).data, :timing) || throw(ArgumentError(
                "benchmark requires timing=true in the computation options"))
            recorded = details(value).data.timing
            return value isa MonteCarloResult ?
                   [NamedTuple[deepcopy(record) for record in population]
                    for population in recorded] : deepcopy(recorded)
        elseif value isa
               Union{AbstractVector{<:LineParameters}, ParametricResult, LinearErrorResult}
            isempty(value) &&
                throw(ArgumentError("benchmark requires nonempty timing=true results"))
            return NamedTuple[measurement(core) for core in value]
        end
        throw(ArgumentError("benchmark requires supported line-parameter results with timing=true; got $(typeof(value))"))
    end
    result = f()
    first_timing = measurement(result)
    timings = first_timing isa NamedTuple ? Vector{NamedTuple}(undef, samples) :
              Vector{typeof(first_timing)}(undef, samples)
    timings[1] = first_timing
    for index in 2:samples
        # Drop the previous large result before starting another repetition.
        result = nothing
        result = f()
        timings[index] = measurement(result)
    end
    return (; result, timings)
end
