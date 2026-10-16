@testitem "LineModels / propagation problem validation" tags=[:unit, :engine] setup=[TestFixtures] begin
    import LineCableModels: validate

    phase = TestFixtures.two_conductor_results()
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation())
    segment = PropagationParameters(modal; line_length=1000.0)
    problem = @inferred LineModelProblem(segment)
    @test problem.segment === segment
    @test validate(problem) === problem
    @test eltype(problem) === eltype(modal)
    @test eltype(typeof(problem)) === eltype(problem)
    @test_throws DomainError LineModelProblem(PropagationParameters(segment; line_length=0.0))
    for f in ([0.0, 10.0, 100.0], [-1.0, 10.0, 100.0],
            [10.0, 10.0, 100.0], [100.0, 10.0, 1000.0],
            [10.0, 100.0, Inf], [10.0, 100.0, NaN])
        changed = deepcopy(segment)
        changed.parameters.f .= f
        @test_throws Union{ArgumentError, DomainError} LineModelProblem(changed)
    end
    empty = segment[Int[]]
    @test_throws ArgumentError LineModelProblem(empty)
    for component in (:Z, :Y)
        changed = deepcopy(segment)
        getproperty(changed.parameters, component).values[1] = Inf
        @test_throws DomainError LineModelProblem(changed)
    end
    for component in (:Tv, :Ti)
        changed = deepcopy(segment)
        getproperty(changed.parameters.domain.operators, component)[1] = NaN
        @test_throws DomainError LineModelProblem(changed)
    end
    changed = deepcopy(segment)
    changed.parameters.domain.gamma[1] = NaN
    @test_throws DomainError LineModelProblem(changed)
    shifted = PropagationParameters(segment.parameters, segment.line_length,
        ComputationDetails(merge(details(segment).data, (inputs=(sigma=1.0,),))))
    @test_throws ArgumentError("line models are fitted on the real frequency axis") LineModelProblem(shifted)
end

@testitem "LineModels / completed segment transport and fitting fallbacks" tags=[:unit, :parametric] setup=[TestFixtures] begin
    import LineCableModels.LineModels as LM

    source = Gridspace{LineParametersProblem}(rho ->
        TestFixtures.three_bare_wires_problem(; rho, frequencies=[10.0, 100.0]),
        (Grid((10.0, 100.0)),))
    phase = compute(source, Formulation())
    modal = compute(Gridspace{ModalAnalysisProblem}(phase), ModalAnalysisFormulation())
    segments = ParametricResult(Combinatorial(ModalAnalysisFormulation()),
        collect(Gridspace{PropagationParameters}(modal)))
    problems = Gridspace{LineModelProblem}(segments)
    @test length(problems) == length(segments) == 2
    @test [point.segment for point in problems] == collect(segments)
    for identifier in (:ulm, :fdcm)
        selected = LineModelFormulation(identifier)
        for point in problems
            @test_throws ArgumentError("load UniversalLineModel to fit :$identifier") invoke(
                LM.fit, Tuple{LM.Formula{identifier}, Vararg{Any}},
                selected.formula, point, selected.formula.parameters, selected.formula.options)
        end
    end
    @test_throws ArgumentError("not implemented yet") compute(first(problems), LineModelFormulation(:marti))
    @test_throws ArgumentError("not implemented yet") LM.fit(LM.Formula(:marti), nothing, (;), FormulationOptions())
    @test computation_options(LineModelFormulation, ComputationOptions()).data ==
        (verbosity=(default=0,), timing=false)
    @test_throws ArgumentError computation_options(LineModelFormulation, ComputationOptions(; unknown=true))
    @test_throws TypeError computation_options(LineModelFormulation, ComputationOptions(; timing=1))
    @test_throws ArgumentError computation_options(LineModelFormulation, ComputationOptions(; verbosity=(default=3,)))
    buffers = (values=zeros(2),)
    @test LM.initialize_buffers(LM.Formula(:ulm), Float64, nothing, nothing, buffers) === buffers
end
