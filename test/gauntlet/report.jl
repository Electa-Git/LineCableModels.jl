using DataFrames

include("artifacts.jl")
using .GauntletArtifacts

length(ARGS) == 1 || throw(ArgumentError("usage: report.jl BACKEND"))
frame = report(Symbol(only(ARGS)))
show(stdout, MIME("text/plain"), frame; allrows = true, allcols = true)
println()
