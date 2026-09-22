module ValidationTestRunner

using TestItemRunner, TOML, Test

# Tags describe purpose or environment, never the observed outcome. Explicit
# selectors can reach every item, including those excluded from ordinary runs.
const ORDINARY_EXCLUDED_TAGS = Set((:quality, :aqua, :visual, :core_only,
    :fem_numerical, :pscad_native))

function selection(arguments, directory; excluded=ORDINARY_EXCLUDED_TAGS)
    queries = filter(!=("--list"), arguments)
    any(startswith("--"), queries) && error("Unknown test option; supported option: --list")
    tags = [Symbol(chop(q; head=4, tail=0)) for q in queries if startswith(q, "tag:")]
    names = lowercase.(filter(q -> !startswith(q, "tag:"), queries))
    any(isempty, names) && error("Empty test selector")
    Symbol("") in tags && error("Empty test tag")
    prefix = abspath(directory) * Base.Filesystem.path_separator
    return item -> begin
        startswith(abspath(item.filename), prefix) || return false
        # A filename/name search must not accidentally launch a live station.
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

function run_tests(root; filter=(_ -> true), verbose=true, list=false)
    started = time_ns()
    root = abspath(root)
    validate_config(joinpath(root, "JuliaTestItems.toml"))
    # File selection belongs to the installed runner. In particular, excluded
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
    # announces item starts, so a terminated run reveals its last active scope.
    function testset(description; verbose)
        if any(item -> item.name == description, selected)
            println("Starting [", round((time_ns() - started) / 1e9; digits=2),
                "s] ", description)
            flush(stdout)
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
