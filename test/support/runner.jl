module ValidationTestRunner

using TestItemRunner, TOML, Test

include(joinpath(@__DIR__, "taxonomy.jl"))

# Tags describe purpose or environment, never the observed outcome. Explicit
# selectors can reach every item, including those excluded from ordinary runs.
const ORDINARY_EXCLUDED_TAGS = Set((:quality, :aqua, :visual, :core_only,
    :fem_numerical, :pscad_native))
# The tags of the additional environments. `changed:REF` excludes them too.
const ENVIRONMENT_TAGS = setdiff(ORDINARY_EXCLUDED_TAGS, (:quality,))

# The paths changed since the git revision `reference`: modified, added, deleted and
# untracked files, relative to `root`. Both sides of a rename count.
function changed_paths(reference, root)
    git(arguments...) = Cmd(["git", "-C", root, arguments...])
    success(pipeline(git("rev-parse", "--verify", "--quiet", reference * "^{commit}");
        stdout=devnull, stderr=devnull)) || error("Unknown git revision for changed: $reference")
    tracked = readlines(git("diff", "--name-only", "--no-renames", "--relative", reference))
    untracked = readlines(git("ls-files", "--others", "--exclude-standard"))
    return sort!(unique!([tracked; untracked]))
end

# The earliest owner of the changed `src/` and `ext/` paths, or nothing.
source_owner(paths, owners) =
    earliest(Symbol[owner for owner in (path_owner(owners, p) for p in paths) if owner !== nothing])

# The items that changes to `paths` can affect, as (file, name) pairs:
# - every `quality` item.
# - every non-slow item whose owner is the earliest owner of a changed `src/` or `ext/`
#   path, or a later owner.
# - every item of a changed test file, slow or not, and of the files that use a setup
#   defined in a changed test file.
# - after a change under `test/support/`, every non-slow item and every item under
#   `test/unit/core/`, which checks the runner.
# Items with an environment tag remain excluded, as in ordinary runs.
function changed_items(paths, items, setups, owners)
    floor = source_owner(paths, owners)
    support = any(startswith("test/support/"), paths)
    files = Set(paths)
    shared = Set(name for (name, file) in setups
        if file in files && !startswith(file, "test/support/"))
    selected = Set{Tuple{String, String}}()
    for item in items
        any(in(ENVIRONMENT_TAGS), item.tags) && continue
        owner = owner_tag(item.tags)
        affected = :quality in item.tags || item.file in files ||
            any(in(shared), item.setups) ||
            (support && startswith(item.file, "test/unit/core/")) ||
            (:slow ∉ item.tags && (support || (floor !== nothing && owner !== nothing &&
                rank(owner) >= rank(floor))))
        affected && push!(selected, (item.file, item.name))
    end
    return selected
end

# The non-quality items `changed:REF` selects. Its quality items are all of them.
function changed_selection(reference, root)
    paths = changed_paths(reference, root)
    (; items, setups) = inventory(root)
    owners = loaded_owners(root)
    selected = changed_items(paths, items, setups, owners)
    println("changed:", reference, ": ", length(paths), " changed paths; earliest source owner ",
        something(source_owner(paths, owners), "none"), "; ",
        count(p -> any(i -> i.file == p, items), paths), " changed test files; test/support ",
        any(startswith("test/support/"), paths) ? "changed" : "unchanged")
    return setdiff(selected, Set((i.file, i.name) for i in items if :quality in i.tags))
end

# The revision of a `changed:REF` selector, or nothing. It takes no other selector.
function changed_reference(arguments)
    queries = filter(!=("--list"), arguments)
    any(startswith("changed:"), queries) || return nothing
    length(queries) == 1 || error("changed:REF cannot be combined with other selectors")
    reference = chop(only(queries); head=8, tail=0)
    isempty(reference) && error("Empty changed: revision")
    return reference
end

# Runs the `ordinary` items, as (file, name) pairs, in this process, then the quality
# selectors in a fresh process, exactly as `tag:quality` runs alone: the reflection
# guards read the live method table, which earlier items may extend. Fails if either
# part fails.
function run_changed(root, ordinary; list=false, quality=["tag:quality"])
    failure = nothing
    if isempty(ordinary)
        println("No items outside quality selected")
    else
        try
            run_tests(root; list, verbose=true, filter=item ->
                (relative(item.filename, root), String(item.name)) in ordinary)
        catch caught
            caught isa InterruptException && rethrow()
            failure = caught
        end
    end
    arguments = [joinpath(root, "test", "runtests.jl"); quality; "--workers=1";
        list ? ["--list"] : String[]]
    command = `$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) $arguments`
    println("Quality items in a fresh process: ", join(quality, " "))
    flush(stdout)
    passed = success(run(ignorestatus(command)))
    failure === nothing || throw(failure)
    passed || error("The quality items failed in their fresh process")
    return nothing
end

# Packs the test files of `items` onto at most `count` workers, longest first by the
# summed `weight` of their items. `items` are inventory items in inventory order. The
# items of a file stay on one worker. The files with a `:fem` item share one worker,
# because they use the FEM runs root, which has no setting. Each worker's items keep the
# order of `items`. One worker runs today's order.
function partition(items, weight, count)
    files = unique(item.file for item in items)
    fem = Set(item.file for item in items if :fem in item.tags)
    units = [[file] for file in files if file ∉ fem]
    isempty(fem) || push!(units, [file for file in files if file in fem])
    load = Dict(file => 0.0 for file in files)
    for item in items
        load[item.file] += weight(item)
    end
    workers = [Set{String}() for _ in 1:min(count, length(units))]
    totals = zeros(length(workers))
    # The sort is stable: units of equal load keep the order of their files.
    for unit in sort(units; by = unit -> sum(file -> load[file], unit), rev = true)
        k = argmin(totals)
        union!(workers[k], unit)
        totals[k] += sum(file -> load[file], unit)
    end
    return [[item for item in items if item.file in worker] for worker in workers]
end

# Runs the `selected` inventory items in fresh processes that start together. `count`
# workers take the ordinary items, which `partition` packs by the item durations that
# earlier runs measured, or by item count without them. One more process takes the
# quality items, as `tag:quality` runs alone, or two when one would take longer than the
# longest worker.
#
# Each process has its own TMPDIR, made with `mktempdir()` and removed afterwards. Its log
# under `.linecablemodels/test/<run>/` keeps the `Starting` lines that
# `test/tools/durations.jl` reads. The console shows each process's item starts and
# failures, then a summary. The measured durations go to
# `.linecablemodels/test/durations.tsv`. Fails if any process fails.
function run_workers(root, selected, count)
    isempty(selected) &&
        error("No test items selected in $root; check selectors and excluded tags")
    started = time()
    output = joinpath(root, ".linecablemodels", "test")
    record = joinpath(output, "durations.tsv")
    measured = Dict{Tuple{String, String}, Float64}()
    if isfile(record)
        for line in eachline(record)
            file, name, seconds = split(line, '\t')
            measured[(file, name)] = parse(Float64, seconds)
        end
    end
    # An item without a measurement weighs as much as the median measured item.
    typical = isempty(measured) ? 1.0 : sort!(collect(values(measured)))[cld(length(measured), 2)]
    weight = item -> get(measured, (item.file, item.name), typical)
    estimate = items -> sum(weight, items; init = 0.0)
    workers = partition([item for item in selected if :quality ∉ item.tags], weight, count)
    groups = [("worker $k", items) for (k, items) in enumerate(workers)]
    # The quality items take a second process only when one would outlast every worker.
    quality = [item for item in selected if :quality in item.tags]
    divided = estimate(quality) > maximum(estimate, workers; init = 0.0)
    for (k, items) in enumerate(partition(quality, weight, divided ? 2 : 1))
        push!(groups, (divided ? "quality $k" : "quality", items))
    end
    directory = mkpath(joinpath(output, Libc.strftime("%Y%m%d-%H%M%S", time()) * "-$(getpid())"))
    println("Fresh processes: ", length(groups), "; logs in ", directory)
    for (label, items) in groups
        println(label, ": ", length(items), " items in ", length(unique(i.file for i in items)),
            " files", isempty(measured) ? "" : "; estimated $(round(Int, estimate(items))) s")
    end
    flush(stdout)
    project = dirname(Base.active_project())
    runner = joinpath(root, "test", "support", "runner.jl")
    temporary = mktempdir()
    processes = Base.Process[]
    results = try
        tasks = map(groups) do (label, items)
            name = replace(label, ' ' => '-')
            log = joinpath(directory, name * ".log")
            program = """
                using TestItemRunner
                include($(repr(runner)))
                const SELECTED = $(repr(Set((item.file, item.name) for item in items)))
                ValidationTestRunner.run_tests($(repr(root)); filter = item ->
                    (ValidationTestRunner.relative(item.filename, $(repr(root))),
                        String(item.name)) in SELECTED)
                """
            command = addenv(`$(Base.julia_cmd()) --project=$project -e $program`,
                "TMPDIR" => mkdir(joinpath(temporary, name)))
            pipe = Pipe()
            process = run(pipeline(command; stdout = pipe, stderr = pipe); wait = false)
            push!(processes, process)
            close(pipe.in)
            @async begin
                named = Dict(item.name => item for item in items)
                current, since = nothing, time()
                open(log, "w") do io
                    for line in eachline(pipe)
                        println(io, line)
                        flush(io)
                        m = match(r"^Starting \[[0-9.]+s\] (.+)$", line)
                        if m !== nothing
                            now = time()
                            current === nothing ||
                                (measured[(current.file, current.name)] = now - since)
                            current, since = get(named, m[1], nothing), now
                        end
                        if m !== nothing || startswith(line, "Test Failed at") ||
                           startswith(line, "Error During Test at")
                            println("[", label, "] ", line)
                            flush(stdout)
                        end
                    end
                end
                wait(process)
                current === nothing || (measured[(current.file, current.name)] = time() - since)
                (; label, items, log, passed = success(process), seconds = time() - started)
            end
        end
        fetch.(tasks)
    finally
        foreach(process -> process_running(process) && kill(process), processes)
        rm(temporary; recursive = true, force = true)
    end
    partial = tempname(output)
    open(partial, "w") do io
        for ((file, name), seconds) in sort!(collect(measured))
            println(io, file, '\t', name, '\t', round(seconds; digits = 2))
        end
    end
    mv(partial, record; force = true)
    for result in results
        println(result.label, ": ", length(result.items), " items in ",
            length(unique(i.file for i in result.items)), " files; ",
            result.passed ? "run completed" : "run failed", "; ",
            round(result.seconds; digits = 2), " s; ", result.log)
    end
    if !isempty(workers)
        seconds = [result.seconds for result in results[eachindex(workers)]]
        println("Longest worker ", round(maximum(seconds); digits = 2), " s; mean ",
            round(sum(seconds) / length(seconds); digits = 2), " s")
    end
    passed = all(result -> result.passed, results)
    println("Selected ", length(selected), " maintained test items in ",
        length(unique(item.file for item in selected)), " files; ",
        passed ? "run completed" : "run failed", "; ", round(time() - started; digits = 2),
        " s elapsed")
    flush(stdout)
    passed || error("Failed in ", join((r.label for r in results if !r.passed), ", "),
        "; see their logs in ", directory)
    return nothing
end

# Runs the selection that `arguments` names, from the test directory `directory`.
# `--workers=N` runs it in N fresh worker processes (`run_workers`). The default is one
# under CI, which runs it in this process as before, and otherwise a quarter of the CPU
# threads, at most 8. `--list` always lists in this process.
function execute(arguments, directory)
    root = dirname(rstrip(abspath(directory), ['/', '\\']))
    options = filter(startswith("--workers="), arguments)
    length(options) <= 1 || error("Give --workers=N once")
    count = isempty(options) ?
            (haskey(ENV, "CI") ? 1 : max(1, min(8, Sys.CPU_THREADS ÷ 4))) :
            something(tryparse(Int, chop(only(options); head = 10, tail = 0)), 0)
    count >= 1 || error("--workers=N takes a positive integer")
    arguments = filter(∉(options), arguments)
    list = "--list" in arguments
    reference = changed_reference(arguments)
    if count == 1 || list
        reference === nothing &&
            return run_tests(root; verbose=true, list, filter=selection(arguments, directory))
        return run_changed(root, changed_selection(reference, root); list)
    end
    items = inventory(root).items
    selected = if reference === nothing
        select = selection(arguments, directory)
        [item for item in items if select((; filename = joinpath(root, item.file),
            item.name, item.tags))]
    else
        ordinary = changed_selection(reference, root)
        [item for item in items if :quality in item.tags || (item.file, item.name) in ordinary]
    end
    return run_workers(root, selected, count)
end

function selection(arguments, directory; excluded=ORDINARY_EXCLUDED_TAGS)
    queries = filter(!=("--list"), arguments)
    any(startswith("--"), queries) &&
        error("Unknown test option; supported options: --list and --workers=N")
    changed_reference(arguments) === nothing || error("changed:REF runs through `execute`, ",
        "which runs its quality items in a fresh process")
    tags = [Symbol(chop(q; head=4, tail=0)) for q in queries if startswith(q, "tag:")]
    names = lowercase.(filter(q -> !startswith(q, "tag:"), queries))
    any(isempty, names) && error("Empty test selector")
    Symbol("") in tags && error("Empty test tag")
    prefix = abspath(directory) * Base.Filesystem.path_separator
    return item -> begin
        startswith(abspath(item.filename), prefix) || return false
        # A search by filename or item name must not accidentally launch a live station.
        :pscad_native in item.tags && :pscad_native ∉ tags && return false
        isempty(queries) && return isempty(intersect(excluded, item.tags))
        tag_match = isempty(tags) || any(in(item.tags), tags)
        name_match = isempty(names) || any(names) do query
            occursin(query, lowercase(relpath(item.filename, directory))) ||
                occursin(query, lowercase(String(item.name)))
        end
        tag_match && name_match
    end
end

# The items and setups under `root/test`, read without running them. Each item
# has its file relative to `root`, name, tags and setups. `setups` maps each setup
# name to the file that defines it.
function inventory(root)
    root = abspath(root)
    prefix = joinpath(root, "test", "")
    items = @NamedTuple{file::String, name::String, tags::Vector{Symbol},
        setups::Vector{Symbol}}[]
    setups = Dict{Symbol, String}()
    for file in TestItemRunner.find_test_files(root)
        startswith(abspath(file), prefix) || continue
        tree = TestItemRunner.JuliaSyntax.parseall(TestItemRunner.JuliaSyntax.SyntaxNode,
            read(file, String); filename=file)
        found, defined, errors = [], [], []
        TestItemRunner.TestItemDetection.find_test_detail!(tree, found, defined, errors)
        isempty(errors) || error("Invalid test item or setup definition in $file")
        path = relative(file, root)
        for item in found
            push!(items, (; file=path, name=String(item.name), tags=item.option_tags,
                setups=item.option_setup))
        end
        for setup in defined
            setups[Symbol(setup.name)] = path
        end
    end
    return (; items, setups)
end

function validate_config(path)
    isfile(path) || return
    document = TOML.parsefile(path)
    unknown = setdiff(keys(document), ("config-version", "include", "exclude"))
    isempty(unknown) || error("Invalid test discovery keys in $path: $(join(unknown, ", "))")
    version = get(document, "config-version", 1)
    version isa Integer && !(version isa Bool) && version == 1 ||
        error("Unsupported test discovery config-version in $path: $version")
    for key in ("include", "exclude")
        patterns = get(document, key, String[])
        patterns isa AbstractVector && all(p -> p isa String, patterns) ||
            error("Test discovery $key must be an array of strings in $path")
    end
    return
end

# `on_start` receives each item name as the item starts.
# `test/tools/owners.jl` uses it to record coverage per item.
function run_tests(root; filter=(_ -> true), verbose=true, list=false, on_start=(_ -> nothing))
    started = time_ns()
    root = abspath(root)
    validate_config(joinpath(root, "JuliaTestItems.toml"))
    # File selection belongs to the installed runner. Excluded
    # source cannot introduce a setup even when an included item requests it.
    files = TestItemRunner.find_test_files(root)
    checked = Set([root])
    for file in files
        directory = dirname(file)
        while directory != root
            if directory ∉ checked
                validate_config(joinpath(directory, "JuliaTestItems.toml"))
                push!(checked, directory)
            end
            directory = dirname(directory)
        end
        # Item extraction alone can overlook invalid source after an item.
        TestItemRunner.JuliaSyntax.parseall(TestItemRunner.JuliaSyntax.SyntaxNode,
            read(file, String); filename=file)
    end
    selected = NamedTuple[]
    finished = false
    # Keep native Test reporting and failure propagation. This factory only
    # announces item starts, which identifies its last active scope.
    function testset(description; verbose)
        if any(item -> item.name == description, selected)
            println("Starting [", round((time_ns() - started) / 1e9; digits=2),
                "s] ", description)
            flush(stdout)
            on_start(description)
        end
        Test.DefaultTestSet(description; verbose)
    end
    try
        result = TestItemRunner.run_tests(root; verbose, testset, filter=item -> begin
            accepted = filter(item)
            if accepted
                push!(selected, (; name=String(item.name), file=relpath(item.filename, root)))
                list && println(relpath(item.filename, root), " | ", item.name,
                    " | ", join(string.(item.tags), ","))
            end
            accepted && !list
        end)
        isempty(selected) && error("No test items selected in $root; check selectors and excluded tags")
        finished = true
        return result
    finally
        println(list ? "Listed " : "Selected ", length(selected), " maintained test items in ",
            length(unique(item.file for item in selected)), " files; ",
            list ? "no test bodies executed" : finished ? "run completed" : "run failed", "; ",
            round((time_ns() - started) / 1e9; digits=2), " s elapsed")
        flush(stdout)
    end
end

end
