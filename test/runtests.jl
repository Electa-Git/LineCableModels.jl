using TestItemRunner
include(joinpath(@__DIR__, "support", "runner.jl"))

ValidationTestRunner.run_tests(dirname(@__DIR__);
    filter=ValidationTestRunner.selection(ARGS, @__DIR__),
    list="--list" in ARGS, verbose=true)
