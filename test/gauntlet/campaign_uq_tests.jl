@testitem "Gauntlet / UQ calculations retain deterministic seeds and recover" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels, Measurements
    model=load_case(:two_insulated_wires;variation=compose_variations(
        ExactOverrides(frequencies=[1.,37.]),RelativeStandardUncertainty(2.;tags=(:geometry,:cable_layer))))
    problem=ParametricProblem(model.problem)
    inner=uq_inner_formulation()
    reference=BenchmarkCalculation(:linear,problem,LinearError(inner))
    candidate=BenchmarkCalculation(:sampled,problem,MonteCarlo(inner;trials=8,seed=UInt64(0x1234),
        return_samples=false,return_histograms=false))
    definition=benchmark_definition(:uncertainty,model.id,:fixture,@__FILE__,model,reference,candidate,(quantities=(:R, :L, :C, :G), statistics=(:mean, :std)),(;))
    mktempdir() do directory
        result=only(run_campaign(directory,[definition];on_error=:fail)).result
        @test frequencies(only(result.reference)) == [1.,37.]
        @test result.candidate_result.root_seed == 0x1234
        @test result.candidate_result.trial_counts == [8]
        recovered=only(resume_campaign(directory))
        @test recovered.skipped && recovered.result===nothing
        retained=read_campaign(directory)[1].candidate
        raw=Gauntlet._read_execution(retained.metadata.path,"result")
        @test statistics(raw) == statistics(result.candidate_result)
        @test retained.metadata.sampling.root_seed == 0x1234
    end
end

@testitem "Gauntlet / joint law calculations retain identity and recover" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using Measurements
    const G = GauntletSupport.Gauntlet
    declaration = G._catalogue_uq_benchmark(:two_insulated_wires,@__FILE__;
        frequencies=[0.1,50.,1e7])
    model = declaration.model
    problem = declaration.reference.problem
    @test declaration.reference.formulation.options.data.distribution === :uniform
    @test declaration.candidate.problem === problem
    reference = BenchmarkCalculation(:monte_carlo,problem,MonteCarlo(uq_inner_formulation();
        trials=8,seed=0x1234,distribution=:uniform))
    candidate = BenchmarkCalculation(:linear_error,problem,LinearError(uq_inner_formulation()))
    definition = benchmark_definition(:joint,model.id,:fixture,@__FILE__,model,
        reference,candidate,(quantities=(:R,:L,:C,:G),statistics=(:mean,:std)),(;))
    old = load_case(:two_insulated_wires;variation=compose_variations(
        ExactOverrides(frequencies=[0.1,50.,1e7]),
        RelativeStandardUncertainty(10.;tags=(:geometry,:cable_layer))))
    old_candidate = BenchmarkCalculation(:linear_error,ParametricProblem(old.problem),candidate.formulation)
    @test calculation_record(old_candidate).input_sha256 != calculation_record(candidate).input_sha256
    mktempdir() do directory
        result = run_benchmark(definition;directory)
        @test result.reference_result.trial_counts == [8]
        for value in (result.reference_result,result.candidate_result),
                quantity in (R,L,C,LineCableModels.G)
            @test all(isfinite,observe(only(value.values),quantity))
        end
        restored = run_benchmark(definition;directory)
        @test restored.timings.execution.reference.reused
        @test restored.timings.execution.candidate.reused
        @test_throws r"inputs changed" G._execute(old_candidate;
            directory=joinpath(directory,"candidate"),model=old)
    end
end
