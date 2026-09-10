@testitem "UQ benchmark / 380 kV 2000 mm² flat vertical / LEP versus Monte Carlo" tags=[:gauntlet, :uq] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels
    # Catalogue construction is independent from native availability. Live campaigns are explicit CLI actions.
    benchmark=benchmark_definition(:benchmark_380kv_2000mm2_flatver_lep_montecarlo)
    @test benchmark isa BenchmarkDefinition
    @test benchmark.reference.id != benchmark.candidate.id
    @test benchmark.case_id === :cable_380kv_2000mm2_flatver
    @test benchmark.reference.problem === benchmark.candidate.problem
end
