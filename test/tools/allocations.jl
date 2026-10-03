# Allocation counts and bytes of the preservation corpus (L2). Run in a fresh process,
# from the repository root:
#
#     julia --project=test test/tools/allocations.jl
#
# It prints the `[allocations]` rows of `test/quality/preservation.toml`. For each
# scenario and frequency count, each `compute` call is warmed up twice and then measured
# three times inside a function. `allocations` is the number of allocated objects (pool,
# big and malloc'd), which must be identical across the measured calls; `bytes` is the
# minimum over them, because the runtime's byte accounting of malloc'd buffers adds a few
# bytes on some calls. The counts depend on what the process computed before, so the
# corpus always runs in the same order in a fresh process. The quality item
# `Quality / preservation / L2 allocation ceilings` runs this file.
using LineCableModels, Measurements

Base.include(@__MODULE__, joinpath(@__DIR__, "..", "support", "scenarios.jl"))
using .CurrentScenarios: preservation_corpus

const FREQUENCIES = (2, 4)
const WARMUP = 2
const CALLS = 3

# The objects and bytes of `CALLS` warmed-up calls of `compute` on `arguments`.
function measure(arguments)
    for _ in 1:WARMUP
        compute(arguments...)
    end
    objects, bytes = Int[], Int[]
    for _ in 1:CALLS
        before = Base.gc_num()
        compute(arguments...)
        difference = Base.GC_Diff(Base.gc_num(), before)
        push!(objects, difference.poolalloc + difference.bigalloc + difference.malloc)
        push!(bytes, difference.allocd)
    end
    return (; objects, bytes)
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

# The table rows, as key => value pairs, and the scenarios whose object counts varied.
function rows(frequencies = FREQUENCIES)
    found = Pair{String, Int}[]
    varied = String[]
    for n in frequencies
        corpus = preservation_corpus(n)
        for name in keys(corpus)
            pairs, variation = scenario_rows(name, n, measure(corpus[name]))
            append!(found, pairs)
            variation === nothing || push!(varied, variation)
        end
    end
    return (; rows = found, varied)
end

function main()
    result = rows()
    println("[allocations]")
    foreach(row -> println(repr(first(row)), " = ", last(row)), result.rows)
    isempty(result.varied) && return 0
    println(stderr, "Nondeterministic allocation counts:")
    foreach(line -> println(stderr, "  ", line), result.varied)
    return 1
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main())
