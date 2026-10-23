using TestItemRunner
include(joinpath(@__DIR__, "support", "runner.jl"))

# The Python environment of the PSCAD extension stays out of the tracked folders. Workers
# inherit it.
haskey(ENV, "JULIA_CONDAPKG_ENV") ||
    (ENV["JULIA_CONDAPKG_ENV"] = joinpath(dirname(@__DIR__), ".linecablemodels", "condapkg"))

ValidationTestRunner.execute(ARGS, @__DIR__)
