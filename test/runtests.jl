using TestItemRunner
include(joinpath(@__DIR__, "support", "runner.jl"))

ValidationTestRunner.execute(ARGS, @__DIR__)
