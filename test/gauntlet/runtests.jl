using TestItemRunner, LineCableModels
include(joinpath(@__DIR__, "..", "support", "runner.jl"))
ValidationTestRunner.run_tests(pkgdir(LineCableModels);
    filter=ValidationTestRunner.selection(ARGS, @__DIR__; excluded=Set{Symbol}()),
    list="--list" in ARGS, verbose=true)
