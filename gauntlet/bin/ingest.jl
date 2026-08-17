using Gauntlet

length(ARGS) in (2, 4) || error(
    "usage: ingest.jl <PSCAD-reference-v1> <destination> [--amendments <path>]",
)
amendments = length(ARGS) == 4 && ARGS[3] == "--amendments" ? ARGS[4] : nothing
length(ARGS) == 2 || amendments !== nothing || error("invalid ingestion arguments")
result = ingest(:pscad, ARGS[1], ARGS[2]; amendments)
display(result)
