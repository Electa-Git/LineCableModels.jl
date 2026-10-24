# Per-item durations from the logs of complete test runs. Run
# `julia --project=test test/tools/durations.jl LOG...` from the repository root.
#
# The runner prints `Starting [t s] NAME` as each item starts and the elapsed time
# at the end. An item lasts until the next one starts. Its duration includes
# compilation and the effects of the items run before it. The report lists the
# items above 10 s without the `slow` tag, the `slow` items below 5 s and the
# minutes per owner tag. Findings are reported, never applied. Run it on the full
# ordinary log and on the environment logs at every track end.
using TestItemRunner
include(joinpath(@__DIR__, "..", "support", "runner.jl"))

const REPOSITORY = dirname(dirname(@__DIR__))
const SLOW = 10.0
const FAST = 5.0
const STARTING = r"^Starting \[([0-9.]+)s\] (.+)$"
const ELAPSED = r"maintained test items in \d+ files; run \w+; ([0-9.]+) s elapsed$"

# Item name => seconds. The last item of a log without a final line has no duration.
function durations(log)
    starts = Tuple{String, Float64}[]
    finish = nothing
    for line in eachline(log)
        if (m = match(STARTING, line)) !== nothing
            push!(starts, (m[2], parse(Float64, m[1])))
        elseif (m = match(ELAPSED, line)) !== nothing
            finish = parse(Float64, m[1])
        end
    end
    isempty(starts) && error("No item starts in $log")
    finish === nothing && println("No final line in $log; its last item has no duration.")
    found = Dict{String, Float64}()
    for (i, (name, start)) in enumerate(starts)
        stop = i < length(starts) ? starts[i+1][2] : finish
        stop === nothing || (found[name] = stop - start)
    end
    return found
end

minutes(seconds) = round(seconds / 60; digits = 1)

function report(logs)
    items = Dict(item.name => item for item in ValidationTestRunner.inventory(REPOSITORY).items)
    measured = Dict{String, Float64}()
    for log in logs
        found = durations(log)
        println(log, ": ", length(found), " items, ", minutes(sum(values(found))), " min")
        for (name, seconds) in found
            haskey(measured, name) && error("$name appears in more than one log")
            measured[name] = seconds
        end
    end
    unknown = sort!([name for name in keys(measured) if !haskey(items, name)])
    isempty(unknown) || println("\nNot in the test files (renamed or removed): ",
        join(unknown, "; "))
    row(name) = string(lpad(round(measured[name]; digits = 1), 7), " s  ",
        items[name].file, " | ", name)
    known = sort!([name for name in keys(measured) if haskey(items, name)];
        by = name -> -measured[name])
    untagged = [name for name in known if measured[name] > SLOW && :slow ∉ items[name].tags]
    println("\nAbove $(Int(SLOW)) s without `slow` (", length(untagged), ", ",
        minutes(sum((measured[name] for name in untagged); init = 0.0)), " min):")
    foreach(name -> println(row(name)), untagged)
    tagged = [name for name in known if measured[name] < FAST && :slow in items[name].tags]
    println("\n`slow` below $(Int(FAST)) s (", length(tagged), "):")
    foreach(name -> println(row(name)), tagged)
    println("\nMinutes per owner tag, `slow` items and the others:")
    owner(name) = something(ValidationTestRunner.owner_tag(items[name].tags), :none)
    for tag in (ValidationTestRunner.OWNERS..., :none)
        names = [name for name in known if owner(name) === tag]
        isempty(names) && continue
        slow = [name for name in names if :slow in items[name].tags]
        println(rpad(tag, 14), lpad(length(names), 4), " items  ", lpad(length(slow), 3),
            " slow  ", lpad(minutes(sum((measured[name] for name in slow); init = 0.0)), 5),
            " min  ", lpad(minutes(sum((measured[name] for name in setdiff(names, slow));
                init = 0.0)), 5), " min")
    end
    return 0
end

function main(arguments)
    isempty(arguments) && error("Usage: julia --project=test test/tools/durations.jl LOG...")
    return report(arguments)
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
