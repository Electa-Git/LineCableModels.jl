# Derive the owner of each test item from the `src/` and `ext/` code it executes.
#
# Record a selection under coverage, from the repository root, once per environment:
#
#     julia --project=test --code-coverage=@$PWD --code-coverage=/tmp/trace.info \
#         test/tools/owners.jl record DIR [SELECTORS...]
#
# Each item start dumps and clears the coverage counters, so `DIR/items.tsv` holds one
# line per item: name, test file and the executed files with their executed line counts.
# The package is loaded before the first item, so loading is not attributed to it. A
# `@testmodule` is evaluated once, by the first item that uses it, and its executions
# count for that item only. The `--code-coverage=FILE.info` flag keeps Julia from
# writing `.cov` files next to the sources at exit.
using TestItemRunner
include(joinpath(@__DIR__, "..", "support", "runner.jl"))

const REPOSITORY = dirname(dirname(@__DIR__))

# The files under `src/` and `ext/` with at least one executed line, and how many.
function executed(tracefile)
    found = Dict{String, Int}()
    file = nothing
    for line in eachline(tracefile)
        if startswith(line, "SF:")
            path = relpath(normpath(line[4:end]), REPOSITORY)
            file = startswith(path, "src/") || startswith(path, "ext/") ? path : nothing
        elseif file !== nothing && startswith(line, "DA:")
            hits = parse(Int, split(line[4:end], ',')[2])
            hits > 0 && (found[file] = get(found, file, 0) + 1)
        end
    end
    return found
end

function record(directory, selectors)
    Base.JLOptions().code_coverage == 0 &&
        error("Run `record` with --code-coverage=@<repository> --code-coverage=<trace>.info")
    mkpath(directory)
    items = joinpath(directory, "items.tsv")
    isfile(items) && error("$items exists; record into an empty directory")
    tracefile = joinpath(directory, "item.info")
    @eval Main using LineCableModels
    files = Dict{String, String}()
    current = Ref{Union{Nothing, String}}(nothing)
    function dump()
        ccall(:jl_write_coverage_data, Cvoid, (Cstring,), tracefile)
        if current[] !== nothing
            found = executed(tracefile)
            line = join(sort!(["$file=$count" for (file, count) in found]), ';')
            open(io -> println(io, current[], '\t', files[current[]], '\t', line), items, "a")
        end
        rm(tracefile; force = true)
        ccall(:jl_clear_coverage_data, Cvoid, ())
    end
    select = ValidationTestRunner.selection(selectors, joinpath(REPOSITORY, "test"))
    filter = item -> begin
        accepted = select(item)
        accepted && (files[String(item.name)] = relpath(item.filename, REPOSITORY))
        accepted
    end
    dump()
    try
        ValidationTestRunner.run_tests(REPOSITORY; filter, verbose = true,
            on_start = name -> (dump(); current[] = name))
    finally
        dump()
    end
end

function main(arguments)
    isempty(arguments) && error("Usage: owners.jl record DIR [SELECTORS...]")
    command, rest = arguments[1], arguments[2:end]
    if command == "record"
        isempty(rest) && error("Usage: owners.jl record DIR [SELECTORS...]")
        record(rest[1], rest[2:end])
    else
        error("Unknown command $command; supported: record")
    end
    return 0
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
