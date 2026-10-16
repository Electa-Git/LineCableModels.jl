@testitem "LineModels / observations, reports and text display" tags=[:unit, :report] setup=[TestFixtures] begin
    using DataFrames
    import LineCableModels.LineModels as LM

    payload = WidebandModel((poles=ComplexF64[-1], residues=ones(ComplexF64, 2, 2, 1), direct=zeros(2, 2)),
        [(delay=1e-6, poles=ComplexF64[-2, -3], residues=ones(ComplexF64, 2, 2, 2))])
    errors = (rms=[0.001 0.002; 0.003 0.004], max=fill(0.005, 2, 2), f_max=fill(100.0, 2, 2))
    result = LineModel(payload, ["west", "east"], 1000.0, [10.0, 100.0],
        (admittance=errors, propagation=deepcopy(errors)),
        (checked=false, passive=nothing, bands=Tuple{Float64, Float64}[], worst=nothing,
            f_worst=nothing, perturbation=nothing),
        (targets=nothing, achieved=nothing), ComputationDetails())
    @test observe(result, poles) == poles(result)
    @test observe(result, delays) == delays(result)
    observed = observables(result, (fitting_error,))
    @test length(observed.quantities) == 2
    for product in observed.quantities
        @test product.values == errors.rms
        @test product.values !== errors.rms
        @test size(product.values) == (2, 2)
    end
    @test !isempty(sprint(show, MIME"text/plain"(), report(observed)))
    @test !isempty(sprint(show, MIME"text/html"(), report(observed)))
    @test length(observables(result).quantities) == 4
    @test !isempty(sprint(show, MIME"text/plain"(), report(observables(result, (poles, delays)))))
    phase = TestFixtures.two_conductor_results()
    segment = PropagationParameters(compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation()); line_length=1000.0)
    problem = LineModelProblem(segment)
    formulation = LineModelFormulation()
    workspace = LM.LineModelWorkspace((coordinates=result.coordinates, length=result.line_length), (;), (;))
    for (value, name) in ((problem, "LineModelProblem"), (formulation, "LineModelFormulation"),
            (workspace, "LineModelWorkspace"), (result, "LineModel"))
        @test occursin(name, sprint(show, value))
        @test !isempty(sprint(summary, value))
        @test !isempty(sprint(show, MIME"text/plain"(), value))
        @test sprint(show, MIME"text/plain"(), value; context=:compact=>true) == sprint(show, value)
    end
end
