using Gauntlet

length(ARGS) == 2 || error("usage: ingest.jl <PSCAD-reference-v1> <destination>")
result = ingest(PSCAD(), ARGS[1], ARGS[2])
display(result)
