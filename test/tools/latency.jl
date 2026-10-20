# First-call latency of the preservation corpus at 2 frequencies. Run
# `julia --project=test --startup-file=no test/tools/latency.jl` from the repository root.
#
# One fresh process times the first call of the scalar coaxial computation and then the
# first call of the chained modal computation. A second fresh process times the first
# call of the product space. This file starts both processes. Each first call includes
# its compilation, and the times leave out the construction of the corpus. The times are
# samples of the machine's load, and checks do not read them.
const SEQUENCES = ("scalar", "space")

if isempty(ARGS)
    project = dirname(Base.active_project())
    for sequence in SEQUENCES
        run(`$(Base.julia_cmd()) --project=$project --startup-file=no $(@__FILE__) $sequence`)
    end
else
    only(ARGS) in SEQUENCES || error("Usage: julia --project=test test/tools/latency.jl")
    using LineCableModels
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "support", "scenarios.jl"))
    using .CurrentScenarios: preservation_corpus
    corpus = preservation_corpus(2)
    if only(ARGS) == "scalar"
        scalar = @elapsed compute(corpus.coaxial...)
        modal = @elapsed compute(corpus.modal...)
        println("scalar first call: ", round(scalar; digits = 1), " s")
        println("chained modal first call after it: ", round(modal; digits = 1), " s")
    else
        space = @elapsed compute(corpus.product...)
        println("product space first call: ", round(space; digits = 1), " s")
    end
end
