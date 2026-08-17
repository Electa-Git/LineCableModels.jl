using Gauntlet

length(ARGS) in (1, 2) || error(
    "usage: stage_smoke.jl <normalized-corpus> [destination]",
)
destination = length(ARGS) == 2 ? ARGS[2] :
              joinpath(@__DIR__, "..", "fixtures", "smoke")
display(stage_smoke(ARGS[1], destination))
