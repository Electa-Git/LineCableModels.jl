using Gauntlet

isempty(ARGS) && error(
    "usage: stage_artifacts.jl <staging-dir> [--repeat-normalized <path>]",
)
repeat_normalized = if length(ARGS) == 3 && ARGS[2] == "--repeat-normalized"
    ARGS[3]
elseif length(ARGS) == 1
    nothing
else
    error("invalid artifact-staging arguments")
end
display(stage_artifacts(:pscad, ARGS[1]; repeat_normalized))
