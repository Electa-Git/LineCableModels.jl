module ValidationTestRunner

using TestItemRunner, TOML

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

function run_tests(root; filter=(_ -> true), verbose=true)
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
    count = Ref(0)
    result = TestItemRunner.run_tests(root; verbose, filter=item -> begin
        accepted = filter(item)
        accepted && (count[] += 1)
        accepted
    end)
    count[] > 0 || error("No test items selected in $root; check selectors and excluded tags")
    println("Selected ", count[], " maintained test items")
    return result
end

end
