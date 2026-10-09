# Timing comparison of the preservation corpus at a git revision and on the working
# tree, on this machine. Run
# `julia --project=test test/tools/performance.jl REF [--frequencies N] [--rounds R] [--seconds S]`
# from the repository root.
#
# REF is checked out in a temporary git worktree. Its test environment starts from this
# repository's `Manifest.toml`, resolved again for REF's projects, so both sides use the
# same dependency versions wherever REF allows them. The report lists any that differ.
# Each side runs its own revision's corpus, `preservation_corpus` in
# `test/support/scenarios.jl` with `N` frequencies (default 4). A revision without one
# uses the working tree's corpus. The report says which copy each side used and whether
# the two corpus files differ.
#
# One worker process per side holds its corpus for the whole run. `taskset` pins both
# workers to the same logical CPU, a performance core on a hybrid machine
# (`/sys/devices/cpu_core/cpus`), never CPU 0 or its sibling. The sides alternate and
# never run together, and they share identical hardware. Without `taskset` they run
# unpinned, with a notice. Each scenario is warmed up once on each side and then sampled
# with BenchmarkTools in `R` rounds (default 6) of `S` seconds per side (default 0.5),
# with a garbage collection before each round. The side that goes first alternates. Load
# on the machine only ever adds time. Each side keeps its minimum over all samples, and
# the medians are printed for information only. The allocation ceilings in
# `test/quality/preservation.toml` guard the allocations behind garbage collection
# exactly.
#
# A scenario whose round minima spread by more than 2 % on either side is unstable and
# has no verdict. Separate builds of identical code differ by a few percent, because each
# worktree builds REF's package image again. A stable scenario more than 10 % slower
# than REF is a slowdown, and one 5 to 10 % slower is a possible slowdown below the
# resolution of the tool. A slowdown, or a scenario that fails on the working tree, exits
# with 1. A possible slowdown or an unstable scenario exits with 2: no verdict, and not a
# regression. Otherwise the tool exits with 0. A scenario that cannot run at REF is
# reported as not comparable. A failing scenario is never retried. The tool prints the
# 1-minute load average at the start and the end. Both revisions run on the machine that
# runs the tool, a developer machine or a CI runner. Other load on that machine makes
# scenarios unstable. The CI workflow `.github/workflows/timing.yml` runs this tool.
const REPOSITORY = dirname(dirname(@__DIR__))
const CORPUS = joinpath("test", "support", "scenarios.jl")
const POSSIBLE = 0.05
const SLOWER = 0.10
const STABILITY = 0.02

# The worker builds its corpus once and prints the scenario names on one line, or
# `failed` and the error. It then answers each request `scenario seconds` with one line:
# the sample times [ns] of that many seconds of samples, or `failed` and the error.
const WORKER = raw"""
using LineCableModels, Measurements, BenchmarkTools
include(ARGS[1])
using .CurrentScenarios: preservation_corpus

flat(error) = first(replace(sprint(showerror, error), r"\s+" => " "), 200)

# A scenario's benchmark, compiled and warmed up on first use.
function benchmark(arguments)
    compute(arguments...)
    # BenchmarkTools collects garbage before each round (`gctrial`), so the scenarios that
    # allocate heavily start every round with samples free of collections.
    trial = @benchmarkable(compute($arguments...), evals = 1, samples = 100_000)
    BenchmarkTools.warmup(trial; verbose = false)
    return trial
end

function serve(corpus)
    benchmarks = Dict{Symbol, Any}()
    for line in eachline(stdin)
        name, seconds = split(line)
        reply = try
            scenario = Symbol(name)
            trial = get!(() -> benchmark(corpus[scenario]), benchmarks, scenario)
            join(run(trial; seconds = parse(Float64, seconds)).times, ' ')
        catch error
            "failed " * flat(error)
        end
        println(reply)
        flush(stdout)
    end
end

corpus = try
    preservation_corpus(parse(Int, ARGS[2]))
catch error
    println("failed ", flat(error))
    exit()
end
println(join(keys(corpus), ' '))
flush(stdout)
serve(corpus)
"""

Base.include(@__MODULE__, joinpath(@__DIR__, "worktree.jl"))

# Logical CPU ids from a list such as "0-15,18".
function cpu_list(text)
    ids = Int[]
    for part in split(strip(text), ','; keepempty = false)
        bounds = parse.(Int, split(part, '-'))
        append!(ids, first(bounds):last(bounds))
    end
    return ids
end

# The logical CPU both workers share, or nothing when `taskset` is missing. It is a
# performance core of a hybrid machine, or any CPU otherwise, and never CPU 0 or its
# sibling.
function shared_cpu(; taskset = Sys.which("taskset"),
        performance = "/sys/devices/cpu_core/cpus",
        siblings = "/sys/devices/system/cpu/cpu0/topology/thread_siblings_list",
        cpus = Sys.CPU_THREADS)
    taskset === nothing && return nothing
    candidates = isfile(performance) ? cpu_list(read(performance, String)) : collect(0:cpus-1)
    avoided = isfile(siblings) ? union(cpu_list(read(siblings, String)), [0]) : [0]
    eligible = setdiff(candidates, avoided)
    return isempty(eligible) ? nothing : first(eligible)
end

struct Worker
    process::Base.Process
    errors::String
    scenarios::Vector{String}
    failure::Union{Nothing, String}
end

# Starts a worker on `corpus` in `project`, pinned to `cpu` unless it is nothing.
# `finish` collects its first line.
function launch(project, corpus, n, cpu)
    errors = tempname()
    command = julia(project, WORKER, corpus, n)
    cpu === nothing || (command = `taskset -c $cpu $command`)
    process = open(pipeline(command; stderr = errors), "r+")
    return (process, errors)
end

function finish((process, errors))
    line = readline(process)
    isempty(line) &&
        error("The timing worker stopped:\n" * read(errors, String))
    startswith(line, "failed ") && return Worker(process, errors, String[], line[8:end])
    return Worker(process, errors, split(line), nothing)
end

function request(worker, scenario, seconds)
    println(worker.process, scenario, " ", seconds)
    flush(worker.process)
    line = readline(worker.process)
    isempty(line) &&
        error("The timing worker stopped:\n" * read(worker.errors, String))
    return startswith(line, "failed ") ? line[8:end] : parse.(Float64, split(line))
end

milliseconds(t) = string(round(t / 1e6; sigdigits = 4), " ms")
median(samples) = (s = sort(samples); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] :
    (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)

# Scenario => side (1 = REF, 2 = working tree) => all samples and the minimum of each
# round, or the error of the first failure. Rounds alternate which side goes first. A
# failed scenario is not retried.
function interleave(workers, rounds, seconds)
    measured = Dict{String, Vector{Any}}()
    for scenario in last(workers).scenarios
        found = Any[(samples = Float64[], minima = Float64[]) for _ in 1:2]
        if !(scenario in first(workers).scenarios)
            found[1] = something(first(workers).failure, "absent")
        end
        for round in 1:rounds, side in (isodd(round) ? (1, 2) : (2, 1))
            found[side] isa String && continue
            reply = request(workers[side], scenario, seconds)
            if reply isa String
                found[side] = reply
            else
                append!(found[side].samples, reply)
                push!(found[side].minima, minimum(reply))
            end
        end
        measured[scenario] = found
    end
    return measured
end

# The verdict on a scenario from its change against REF and the round spreads of both
# sides: `:unstable` when a spread exceeds 2 %, else `:slower` above 10 %, `:possible`
# between 5 and 10 %, otherwise `:same`.
verdict(change, spreads) = maximum(spreads) > STABILITY ? :unstable :
    change > SLOWER ? :slower : change > POSSIBLE ? :possible : :same

const VERDICTS = Dict(:unstable => "  unstable: no verdict",
    :slower => "  slower than $(round(Int, 100SLOWER)) %",
    :possible => "  possible slowdown, below the tool's resolution", :same => "")

# 1 for a slowdown or a failure on the working tree, otherwise 2 for a possible slowdown
# or an unstable scenario, otherwise 0.
status(verdicts) = any(in((:slower, :failed)), verdicts) ? 1 :
    any(in((:possible, :unstable)), verdicts) ? 2 : 0

# The relative spread of the round minima of one side.
spread(minima) = (maximum(minima) - minimum(minima)) / minimum(minima)
percent(x) = string(round(100x; digits = 1), " %")
load() = round(first(Sys.loadavg()); digits = 2)

function compare(reference; n = 4, rounds = 6, seconds = 0.5)
    working = joinpath(REPOSITORY, "test")
    println("Load average (1 min) at the start: ", load())
    return at_revision(reference) do project
        corpus = own_copy(dirname(project), CORPUS; marker = "function preservation_corpus(")
        println(reference, ": ", corpus.own ? "its own corpus" :
            "the working tree's corpus ($reference has no preservation corpus)")
        println("working tree: its own corpus")
        corpus.own && read(corpus.path) != read(joinpath(REPOSITORY, CORPUS)) &&
            println("Notice: the corpus file differs between the two sides.")
        cpu = shared_cpu()
        println(cpu === nothing ?
            "Notice: taskset or an eligible CPU is missing; the workers run unpinned." :
            "Both workers are pinned to logical CPU $cpu.")
        println("Timing $n frequencies in $rounds interleaved rounds of $(seconds) s per side.")
        started = (launch(project, corpus.path, n, cpu),
            launch(working, joinpath(REPOSITORY, CORPUS), n, cpu))
        workers = map(finish, started)
        try
            last(workers).failure === nothing ||
                error("The working tree's corpus fails: " * last(workers).failure)
            measured = interleave(workers, rounds, seconds)
            verdicts = Pair{String, Symbol}[]
            println(rpad("scenario", 17), lpad("$reference min", 14), lpad("working min", 14),
                lpad("change", 10), lpad("spread", 16), "   medians, for information")
            for scenario in last(workers).scenarios
                old, new = measured[scenario]
                if new isa String
                    println(rpad(scenario, 17), "fails on the working tree: ", new)
                    push!(verdicts, scenario => :failed)
                elseif old isa String
                    println(rpad(scenario, 17), lpad("—", 14),
                        lpad(milliseconds(minimum(new.samples)), 14),
                        "  not comparable at $reference: ", old)
                else
                    change = minimum(new.samples) / minimum(old.samples) - 1
                    spreads = (spread(old.minima), spread(new.minima))
                    found = verdict(change, spreads)
                    push!(verdicts, scenario => found)
                    println(rpad(scenario, 17), lpad(milliseconds(minimum(old.samples)), 14),
                        lpad(milliseconds(minimum(new.samples)), 14), lpad(percent(change), 10),
                        lpad(percent(first(spreads)) * " / " * percent(last(spreads)), 16), "   ",
                        milliseconds(median(old.samples)), " / ", milliseconds(median(new.samples)),
                        VERDICTS[found])
                end
            end
            return verdicts
        finally
            for worker in workers
                close(worker.process)
                wait(worker.process)
            end
            println("Load average (1 min) at the end: ", load())
        end
    end
end

function main(arguments)
    usage = "Usage: julia --project=test test/tools/performance.jl REF " *
        "[--frequencies N] [--rounds R] [--seconds S]"
    isempty(arguments) && error(usage)
    reference, options = arguments[1], arguments[2:end]
    n, rounds, seconds = 4, 6, 0.5
    while !isempty(options)
        length(options) >= 2 || error(usage)
        option, value = popfirst!(options), popfirst!(options)
        option == "--frequencies" ? (n = parse(Int, value)) :
            option == "--rounds" ? (rounds = parse(Int, value)) :
            option == "--seconds" ? (seconds = parse(Float64, value)) : error(usage)
    end
    verdicts = compare(reference; n, rounds, seconds)
    named(kinds) = join([s for (s, v) in verdicts if v in kinds], ", ")
    code = status(last.(verdicts))
    code == 1 && println(stderr, "Slower than ", round(Int, 100SLOWER), " % or failing: ",
        named((:slower, :failed)))
    code == 2 && println(stderr, "Possible slowdown or unstable, no verdict: ",
        named((:possible, :unstable)))
    return code
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
