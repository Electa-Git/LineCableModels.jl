@testitem "UQ benchmark / 525 kV 1600 mm² bipole / LEP versus Monte Carlo" tags=[:gauntlet, :uq] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels
    # Catalogue construction is independent from native availability. Live campaigns are explicit CLI actions.
    benchmark=benchmark_definition(:benchmark_525kv_1600mm2_bipole_lep_montecarlo)
    @test benchmark isa BenchmarkDefinition
    @test benchmark.reference.id != benchmark.candidate.id
    @test benchmark.case_id === :cable_525kv_1600mm2_bipole
    @test benchmark.reference.problem === benchmark.candidate.problem
end
