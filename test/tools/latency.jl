# First-call latency of the preservation corpus at 2 frequencies, and of a space and the
# modal analysis of its results. Run `julia --project=test --startup-file=no
# test/tools/latency.jl` from the repository root, or add one sequence name to run that
# sequence alone.
#
# One fresh process times the first call of the scalar coaxial computation and then the
# first call of the chained modal computation. A second fresh process times the first
# call of the product space. A third fresh process, `chain`, builds a formulation that
# returns fixed line parameters and a space of two problems, then times the first call of
# that space and then the first call of the explicit two-step,
# `compute(Gridspace{ModalAnalysisProblem}(result), modal)` on that space's result. This
# file starts the three processes. Each first call includes its compilation, and the
# times leave out the construction of the corpus. The times are samples of the machine's
# load. Only the `chain` sequence has a check: `test/quality/chain_latency.jl` fails when
# its two times together exceed 30 s.
const SEQUENCES = ("scalar", "space", "chain")

if isempty(ARGS)
    project = dirname(Base.active_project())
    for sequence in SEQUENCES
        run(`$(Base.julia_cmd()) --project=$project --startup-file=no $(@__FILE__) $sequence`)
    end
else
    only(ARGS) in SEQUENCES ||
        error("Usage: julia --project=test test/tools/latency.jl [$(join(SEQUENCES, " | "))]")
    using LineCableModels
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "support", "scenarios.jl"))
    using .CurrentScenarios: preservation_corpus, three_bare_wires_problem
    if only(ARGS) == "scalar"
        corpus = preservation_corpus(2)
        scalar = @elapsed compute(corpus.coaxial...)
        modal = @elapsed compute(corpus.modal...)
        println("scalar first call: ", round(scalar; digits = 1), " s")
        println("chained modal first call after it: ", round(modal; digits = 1), " s")
    elseif only(ARGS) == "space"
        corpus = preservation_corpus(2)
        space = @elapsed compute(corpus.product...)
        println("product space first call: ", round(space; digits = 1), " s")
    else
        # A formulation that returns its stored line parameters. Both calls then compile
        # the runner and the modal analysis, not the coaxial backend.
        struct FixedPhase{P <: LineParameters} <: AbstractFormulation
            parameters::P
        end
        LineCableModels.computation_options(::Type{<:FixedPhase}, options::ComputationOptions) =
            options
        LineCableModels.compute(::LineCableModels.Resolved{<:LineParametersProblem},
            formulation::FixedPhase;
            options::Union{NamedTuple, ComputationOptions} = ComputationOptions()) =
            formulation.parameters
        frequencies = [50.0, 80.0]
        fixed = FixedPhase(compute(three_bare_wires_problem(; frequencies), Formulation()))
        space = Gridspace{LineParametersProblem}(
            rho -> three_bare_wires_problem(; rho, frequencies), (Grid((10.0, 100.0)),))
        modal = ModalAnalysisFormulation(:default)
        computed = @timed compute(space, fixed)
        twostep = @elapsed compute(Gridspace{ModalAnalysisProblem}(computed.value), modal)
        println("space first call: ", round(computed.time; digits = 2), " s")
        println("two-step modal analysis of its result, first call after it: ",
            round(twostep; digits = 2), " s")
    end
end
