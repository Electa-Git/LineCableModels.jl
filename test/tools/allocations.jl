# Allocation counts and bytes of the preservation corpus, for the allocation ceilings in
# `test/quality/preservation.toml`. Run `julia --project=test test/tools/allocations.jl`
# in a fresh process from the repository root.
#
# It prints the `[allocations]` rows of `test/quality/preservation.toml`. For each
# scenario and frequency count, each `compute` call is warmed up twice and then measured
# three times inside a function. `allocations` is the number of allocated objects (pool,
# big and allocated by `malloc`), which must be identical across the measured calls.
# `bytes` is the minimum over them, because the runtime's byte accounting of buffers from
# `malloc` adds a few bytes on some calls. The counts depend on what the process computed before, so the
# corpus always runs in the same order in a fresh process. The quality item
# `Quality / preservation / allocation ceilings` runs this file.
#
# It also prints `[allowances]` and `[derivative]`, which `preservation.toml` does not
# record. Measurements stores a partial derivative only when it is nonzero. Whether a
# round-off derivative is exactly zero depends on the machine's last bits. A row's
# allowance is the number of uncertain real scalars that the scenario's result publishes.
# Each of them can store one derivative entry more or less. `[derivative]` gives the bytes
# of one entry.
using LineCableModels, Measurements
using LineCableModels.Commons: AbstractParametricResult, AbstractUncertaintyResult

Base.include(@__MODULE__, joinpath(@__DIR__, "..", "support", "scenarios.jl"))
using .CurrentScenarios: preservation_corpus

const FREQUENCIES = (2, 4)
const WARMUP = 2
const CALLS = 3

# The objects and bytes of `CALLS` warmed-up calls of `compute` on `arguments`, and the
# result of the last warm-up call. The measured calls only compute.
function measure(arguments)
    result = nothing
    for _ in 1:WARMUP
        result = compute(arguments...)
    end
    objects, bytes = Int[], Int[]
    for _ in 1:CALLS
        before = Base.gc_num()
        compute(arguments...)
        difference = Base.GC_Diff(Base.gc_num(), before)
        push!(objects, difference.poolalloc + difference.bigalloc + difference.malloc)
        push!(bytes, difference.allocd)
    end
    return (; objects, bytes, result)
end

# The uncertain real scalars that a result publishes. Real and imaginary parts count
# separately.
allowance(::Measurement) = 1
allowance(::Real) = 0
allowance(value::Complex) = allowance(real(value)) + allowance(imag(value))
allowance(values::AbstractArray) = sum(allowance, values; init = 0)
allowance(result::LineParameters) =
    allowance(result.Z) + allowance(result.Y) + allowance(result.f)
allowance(result::CableConstants) = allowance(result.R) + allowance(result.L) +
    allowance(result.C) + allowance(result.G) + allowance(result.frequency)
allowance(result::Union{AbstractParametricResult, AbstractUncertaintyResult}) =
    allowance(result.values)

# One derivative entry takes the difference in bytes between an uncertain value and an
# exact one.
function derivative_bytes()
    measurement(1.0, 0.1)
    measurement(1.0, 0.0)
    return @allocated(measurement(1.0, 0.1)) - @allocated(measurement(1.0, 0.0))
end

# The rows of one measured scenario, and a description of its varied object counts or
# nothing.
function scenario_rows(name, n, measured)
    prefix = "$name | $n frequencies"
    varied = allequal(measured.objects) ? nothing :
        "$prefix: allocations $(measured.objects) across calls"
    return (["$prefix | allocations" => first(measured.objects),
        "$prefix | bytes" => minimum(measured.bytes)], varied)
end

# The table rows as key => value pairs, the scenarios whose object counts varied, and the
# allowance of each scenario. The allowances come after the whole corpus, so that the
# measured process computes only the corpus.
function rows(counts = FREQUENCIES)
    found = Pair{String, Int}[]
    varied = String[]
    results = Pair{String, Any}[]
    for n in counts
        corpus = preservation_corpus(n)
        for name in keys(corpus)
            measured = measure(corpus[name])
            pairs, variation = scenario_rows(name, n, measured)
            append!(found, pairs)
            variation === nothing || push!(varied, variation)
            push!(results, "$name | $n frequencies" => measured.result)
        end
    end
    allowances = [scenario => allowance(result) for (scenario, result) in results]
    return (; rows = found, varied, allowances)
end

function main()
    result = rows()
    println("[allocations]")
    foreach(row -> println(repr(first(row)), " = ", last(row)), result.rows)
    println("\n[allowances]")
    foreach(row -> println(repr(first(row)), " = ", last(row)), result.allowances)
    println("\n[derivative]\nbytes = ", derivative_bytes())
    isempty(result.varied) && return 0
    println(stderr, "Nondeterministic allocation counts:")
    foreach(line -> println(stderr, "  ", line), result.varied)
    return 1
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main())
