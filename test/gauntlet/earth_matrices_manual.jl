# The superseded model-comparison generator is retired. Current closure and
# independent equal-medium controls live with the maintained engine tests.
# Existing research output under .linecablemodels/qa remains user data.
append!(ARGS,["unit/engine/unified_earth_return"])
include(joinpath(@__DIR__,"../runtests.jl"))
