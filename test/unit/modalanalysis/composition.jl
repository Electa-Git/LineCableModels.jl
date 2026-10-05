@testitem "ModalAnalysis / composition / any formulation that delivers line parameters" tags=[:unit, :parametric] setup=[TestFixtures] begin
    # A test-owned formulation returns the same phase-domain result for every problem.
    struct FixedPhase{P <: LineParameters} <: AbstractFormulation
        parameters::P
    end
    function LineCableModels.compute(::LineParametersProblem, formulation::FixedPhase;
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
    @test same(compute(first(LineCableModels.points(space)), fixed, modal))
    for composed in (compute(space, fixed; modal = :default), compute(space, fixed, modal),
            compute(ParametricProblem(space), Combinatorial(fixed), modal),
            compute(ParametricProblem(space), Combinatorial(fixed); modal = :default))
        @test composed isa ParametricResult
        @test length(composed) == 2
        @test all(same, composed)
    end
    study = compute(ParametricProblem(space), Combinatorial(fixed))
    twostep = compute(Gridspace{ModalAnalysisProblem}(study), modal)
    composed = compute(ParametricProblem(space), Combinatorial(fixed), modal)
    @test [Z(value) for value in composed] == [Z(value) for value in twostep]
    @test [gamma(value) for value in composed] == [gamma(value) for value in twostep]
end
