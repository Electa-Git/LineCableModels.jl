# L3. Timing of the preservation corpus at a git revision and on the working tree, on
# this machine. Run from the repository root:
#
#     julia --project=test test/tools/performance.jl REF [--frequencies N] [--seconds S]
#
# REF is checked out in a temporary git worktree. Its test environment starts from this
# repository's `Manifest.toml`, resolved again for REF's projects, so both sides use the
# same dependency versions wherever REF allows them; the report lists any that differ.
# Each side runs its own revision's corpus (`preservation_corpus` in
# `test/support/scenarios.jl`, `N` frequencies, default 4) in a fresh process; a
# revision without one uses the working tree's, and the report says which copy each side
# used and whether the two corpus files differ. Each scenario gets one warm-up call, then
# BenchmarkTools samples for `S` seconds (default 3).
# The report gives the median time per scenario and its change. A slowdown above 5 % in
# any scenario fails, and so does a scenario that fails on the working tree; a scenario
# that cannot run at REF is reported as not comparable. Not run in CI: timings belong to
# one machine.
const REPOSITORY = dirname(dirname(@__DIR__))
const CORPUS = joinpath("test", "support", "scenarios.jl")
const THRESHOLD = 0.05

# Each side prints one line per scenario: name, median, minimum [ns] and sample count,
# or name, "failed" and the error.
const PROGRAM = raw"""
using LineCableModels, Measurements, BenchmarkTools
corpus, n, seconds = ARGS[1], parse(Int, ARGS[2]), parse(Float64, ARGS[3])
include(corpus)
using .CurrentScenarios: preservation_corpus

function timing(arguments, seconds)
    compute(arguments...)
    trial = run(@benchmarkable(compute($arguments...), evals = 1); seconds, samples = 100_000)
    return (median(trial).time, minimum(trial).time, length(trial.times))
end

function main(n, seconds)
    corpus = try
        preservation_corpus(n)
    catch error
        println("corpus\tfailed\t", replace(sprint(showerror, error), r"\s+" => " "))
        return
    end
    for name in keys(corpus)
        result = try
            timing(corpus[name], seconds)
        catch error
            println(name, "\tfailed\t", first(replace(sprint(showerror, error), r"\s+" => " "), 200))
            continue
        end
        println(name, "\t", join(result, "\t"))
        flush(stdout)
    end
end
main(n, seconds)
"""

Base.include(@__MODULE__, joinpath(@__DIR__, "worktree.jl"))

# Scenario => (median, minimum, samples) in ns, or the error text of a failed scenario,
# in corpus order.
function timings(project, corpus, n, seconds)
    found = Pair{String, Any}[]
    for line in eachline(IOBuffer(output(julia(project, PROGRAM, corpus, n, seconds))))
        fields = split(line, '\t')
        push!(found, String(fields[1]) => (fields[2] == "failed" ? String(fields[3]) :
            (parse(Float64, fields[2]), parse(Float64, fields[3]), parse(Int, fields[4]))))
    end
    return found
end

milliseconds(t) = string(round(t / 1e6; sigdigits = 4), " ms")

function compare(reference; n = 4, seconds = 3.0)
    working = joinpath(REPOSITORY, "test")
    return at_revision(reference) do project
        corpus = own_copy(dirname(project), CORPUS; marker = "function preservation_corpus(")
        println(reference, ": ", corpus.own ? "its own corpus" :
            "the working tree's corpus ($reference has no preservation corpus)")
        println("working tree: its own corpus")
        corpus.own && read(corpus.path) != read(joinpath(REPOSITORY, CORPUS)) &&
            println("Notice: the corpus file differs between the two sides.")
        println("Timing $n frequencies, $(seconds) s per scenario: $reference, then the working tree.")
        old = Dict(timings(project, corpus.path, n, seconds))
        new = timings(working, joinpath(REPOSITORY, CORPUS), n, seconds)
        slower = String[]
        println(rpad("scenario", 17), lpad(reference, 14), lpad("working", 14), lpad("change", 10))
        for (name, measured) in new
            if measured isa String
                println(rpad(name, 17), "fails on the working tree: ", measured)
                push!(slower, name)
            elseif !(get(old, name, "absent") isa Tuple)
                println(rpad(name, 17), lpad("—", 14), lpad(milliseconds(measured[1]), 14),
                    "  not comparable at $reference: ", get(old, name, "absent"))
            else
                change = measured[1] / old[name][1] - 1
                change > THRESHOLD && push!(slower, name)
                println(rpad(name, 17), lpad(milliseconds(old[name][1]), 14),
                    lpad(milliseconds(measured[1]), 14),
                    lpad(string(round(100change; digits = 1), " %"), 10),
                    change > THRESHOLD ? "  slower than $(round(Int, 100THRESHOLD)) %" : "")
            end
        end
        return slower
    end
end

function main(arguments)
    usage = "Usage: julia --project=test test/tools/performance.jl REF [--frequencies N] [--seconds S]"
    isempty(arguments) && error(usage)
    reference, options = arguments[1], arguments[2:end]
    n, seconds = 4, 3.0
    while !isempty(options)
        length(options) >= 2 || error(usage)
        option, value = popfirst!(options), popfirst!(options)
        option == "--frequencies" ? (n = parse(Int, value)) :
            option == "--seconds" ? (seconds = parse(Float64, value)) : error(usage)
    end
    slower = compare(reference; n, seconds)
    isempty(slower) && return 0
    println(stderr, "Slower than ", round(Int, 100THRESHOLD), " % or failing: ", join(slower, ", "))
    return 1
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
