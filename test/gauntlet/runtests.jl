using TestItemRunner, LineCableModels
include(joinpath(@__DIR__, "..", "support", "runner.jl"))
ValidationTestRunner.run_tests(pkgdir(LineCableModels);
    filter = item -> startswith(abspath(item.filename), abspath(@__DIR__) * Base.Filesystem.path_separator) &&
        (:gauntlet in item.tags || :gauntlet_toolkit in item.tags) &&
        (isempty(ARGS) || any(ARGS) do query
            occursin(lowercase(query), lowercase(relpath(item.filename, @__DIR__))) ||
                occursin(lowercase(query), lowercase(String(item.name)))
        end),
    verbose = true)
