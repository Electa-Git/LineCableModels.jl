using TestItemRunner, LineCableModels
TestItemRunner.run_tests(pkgdir(LineCableModels);
    filter = item -> (:gauntlet in item.tags || :gauntlet_toolkit in item.tags) &&
        (isempty(ARGS) || any(ARGS) do query
            occursin(lowercase(query), lowercase(relpath(item.filename, @__DIR__))) ||
                occursin(lowercase(query), lowercase(String(item.name)))
        end),
    verbose = true)
