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
        @test result.reference.frequencies == [1.,37.]
        @test result.candidate_result.root_seed == 0x1234
        @test result.candidate_result.trial_counts == [8]
        recovered=only(resume_campaign(directory)).result
        @test statistics(recovered.candidate_result) == statistics(result.candidate_result)
        @test recovered.timings.execution.candidate.reused
        @test read_campaign(directory)[1].candidate.metadata.sampling.root_seed == 0x1234
    end
end
