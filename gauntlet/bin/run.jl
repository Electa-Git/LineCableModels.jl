using Gauntlet
using TOML

function usage()
    error("usage: run.jl smoke | run.jl full --artifact-dir <path>")
end

isempty(ARGS) && usage()
mode = Symbol(first(ARGS))
corpus = if mode === :smoke
    length(ARGS) == 1 || usage()
    Corpus(:smoke)
elseif mode === :full
    length(ARGS) == 3 && ARGS[2] == "--artifact-dir" || usage()
    Corpus(ARGS[3])
else
    usage()
end

suite = Suite(mode; corpus)
report = gauntlet(suite)
output = joinpath(@__DIR__, "..", "reports")
mkpath(output)
for extension in ("json", "csv", "md")
    write_report(joinpath(output, "$(mode).$(extension)"), report)
end
display(report)
any(
    comparison -> comparison.verdict isa Gauntlet.Fail,
    (comparison for trial in report.trials for comparison in trial.comparisons)
) && error("Gauntlet suite contains failing comparisons; inspect the emitted reports")
