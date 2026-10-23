@testitem "ModalAnalysis / composition / any formulation that delivers line parameters" tags=[:unit, :parametric] setup=[TestFixtures] begin
    # A test-owned formulation returns the same phase-domain result for every problem.
    struct FixedPhase{P <: LineParameters} <: AbstractFormulation
        parameters::P
    end
    LineCableModels.computation_options(::Type{<:FixedPhase}, options::ComputationOptions) = options
    function LineCableModels.compute(::LineCableModels.Resolved{<:LineParametersProblem},
            formulation::FixedPhase;
            options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
        return formulation.parameters
    end
    frequencies = [50.0, 80.0]
    problem = TestFixtures.three_bare_wires_problem(; frequencies)
    phase = compute(problem, Formulation())
    fixed = FixedPhase(phase)
    modal = ModalAnalysisFormulation(:default)
    direct = compute(ModalAnalysisProblem(phase), modal)
    same(value) = Z(value) == Z(direct) && Y(value) == Y(direct) && gamma(value) == gamma(direct)

    @test same(compute(problem, fixed, modal))
    space = Gridspace{LineParametersProblem}(
        rho -> TestFixtures.three_bare_wires_problem(; rho, frequencies),
        (Grid((10.0, 100.0)),))
    @test same(compute(LineCableModels.materialize(first(LineCableModels.points(space))),
        fixed, modal))
    # The modal analysis of the results of a space is a computation of its own.
    twostep = compute(Gridspace{ModalAnalysisProblem}(compute(space, fixed)), modal)
    @test twostep isa ParametricResult
    @test length(twostep) == 2
    @test all(same, twostep)
end
