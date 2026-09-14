# The frozen-table comparison/generator is retired. Run the maintained current
# controls through the strict runner; their output is scoped test evidence.
# Performance measurements have no current numerical-acceptance authority.
push!(ARGS,"unit/engine/unified_earth_return")
include(joinpath(@__DIR__,"..","runtests.jl"))
