@testitem "UQ benchmark / 132 kV 630 mm² / LEP versus Monte Carlo" tags=[:gauntlet, :uq] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels
    # Catalogue construction is independent from native availability. Live campaigns are explicit CLI actions.
    benchmark=benchmark_definition(:benchmark_132kv_630mm2_flathor_lep_montecarlo)
    @test benchmark isa BenchmarkDefinition
    @test benchmark.reference.id != benchmark.candidate.id
    @test benchmark.case_id === :cable_132kv_630mm2_flathor
    @test benchmark.reference.problem === benchmark.candidate.problem
end
