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
# writing `.cov` files next to the sources at exit. Then report:
#
#     julia --project=test test/tools/owners.jl report DIR...
#
# Each executed file counts at its load position (`path_owner` in
# `test/support/taxonomy.jl`); the item's executed owner is the latest of them. The
# report lists every recorded item with its executed owner, its owner tag and the
# executed files at that owner, then the items whose owner tag is earlier than the
# code they execute. Items that run code only in child processes record nothing.
# The report never changes tags; it is optional at a track end.
using TestItemRunner
include(joinpath(@__DIR__, "..", "support", "runner.jl"))

const REPOSITORY = dirname(dirname(@__DIR__))

# The files under `src/` and `ext/` with at least one executed line, and how many.
function executed(tracefile)
    found = Dict{String, Int}()
    file = nothing
    for line in eachline(tracefile)
        if startswith(line, "SF:")
            path = ValidationTestRunner.relative(normpath(line[4:end]), REPOSITORY)
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
        accepted && (files[String(item.name)] =
            ValidationTestRunner.relative(item.filename, REPOSITORY))
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

# Item name => test file, executed owner, and the executed files at that owner with
# their executed line counts, most first.
function executed_owners(directories, owners)
    found = Dict{String, @NamedTuple{file::String, owner::Union{Nothing, Symbol},
        evidence::Vector{Pair{String, Int}}}}()
    for directory in directories, line in eachline(joinpath(directory, "items.tsv"))
        name, file, entries = split(line, '\t')
        files = map(split(entries, ';'; keepempty = false)) do entry
            path, count = rsplit(entry, '='; limit = 2)
            String(path) => parse(Int, count)
        end
        owner = ValidationTestRunner.latest(
            [ValidationTestRunner.path_owner(owners, first(f)) for f in files])
        evidence = sort!([f for f in files
            if ValidationTestRunner.path_owner(owners, first(f)) === owner]; by = last, rev = true)
        haskey(found, name) && error("$name is recorded more than once")
        found[name] = (; file, owner, evidence)
    end
    return found
end

function report(directories)
    owners = ValidationTestRunner.loaded_owners(REPOSITORY)
    recorded = executed_owners(directories, owners)
    items = ValidationTestRunner.inventory(REPOSITORY).items
    tagged = Dict(item.name => ValidationTestRunner.owner_tag(item.tags) for item in items)
    label(owner) = owner === nothing ? "-" : String(owner)
    evidence(record) = join(["$path=$count" for (path, count) in first(record.evidence, 3)], ", ")
    row(name, record) = string(rpad(label(record.owner), 14),
        rpad(label(get(tagged, name, nothing)), 14), record.file, " | ", name,
        isempty(record.evidence) ? "" : "  [" * evidence(record) * "]")
    names = sort!(collect(keys(recorded)); by = name -> (recorded[name].file, name))
    println("Recorded items (", length(names), "): executed owner, owner tag, ",
        "file | name  [executed files at that owner]")
    foreach(name -> println(row(name, recorded[name])), names)
    unknown = [name for name in names if !haskey(tagged, name)]
    isempty(unknown) ||
        println("\nNot in the test files (renamed or removed): ", join(unknown, "; "))
    earlier = [name for name in names if haskey(tagged, name) && recorded[name].owner !== nothing &&
        ValidationTestRunner.rank(something(tagged[name], :none)) <
        ValidationTestRunner.rank(recorded[name].owner)]
    println("\nOwner tag earlier than the executed code (", length(earlier), "):")
    foreach(name -> println(row(name, recorded[name])), earlier)
    return 0
end

const USAGE = "Usage: owners.jl record DIR [SELECTORS...] | owners.jl report DIR..."

function main(arguments)
    length(arguments) >= 2 || error(USAGE)
    command, rest = arguments[1], arguments[2:end]
    command == "record" && return (record(rest[1], rest[2:end]); 0)
    command == "report" && return report(rest)
    error(USAGE)
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
