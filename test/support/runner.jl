module ValidationTestRunner

using TestItemRunner, TOML, Test

include(joinpath(@__DIR__, "taxonomy.jl"))

# Tags describe purpose or environment, never the observed outcome. Explicit
# selectors can reach every item, including those excluded from ordinary runs.
const ORDINARY_EXCLUDED_TAGS = Set((:quality, :aqua, :visual, :core_only,
    :fem_numerical, :pscad_native))
# The tags of the additional environments; `changed:REF` excludes them too.
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
# - every `quality` item;
# - every non-slow item whose owner is the earliest owner of a changed `src/` or `ext/`
#   path, or a later owner;
# - every item of a changed test file, slow or not, and of the files that use a setup
#   defined in a changed test file;
# - after a change under `test/support/`, every non-slow item and every item under
#   `test/unit/core/`, which checks the runner.
# Items with an environment tag stay excluded, as in ordinary runs.
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
    arguments = [joinpath(root, "test", "runtests.jl"); quality; list ? ["--list"] : String[]]
    command = `$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) $arguments`
    println("Quality items in a fresh process: ", join(quality, " "))
    flush(stdout)
    passed = success(run(ignorestatus(command)))
    failure === nothing || throw(failure)
    passed || error("The quality items failed in their fresh process")
    return nothing
end

# Runs the selection that `arguments` names, from the test directory `directory`.
function execute(arguments, directory)
    root = dirname(rstrip(abspath(directory), ['/', '\\']))
    reference = changed_reference(arguments)
    reference === nothing && return run_tests(root; verbose=true, list="--list" in arguments,
        filter=selection(arguments, directory))
    return run_changed(root, changed_selection(reference, root); list="--list" in arguments)
end

function selection(arguments, directory; excluded=ORDINARY_EXCLUDED_TAGS)
    queries = filter(!=("--list"), arguments)
    any(startswith("--"), queries) && error("Unknown test option; supported option: --list")
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
