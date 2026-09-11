@testitem "PSCAD benchmark / 18 kV 1000 mm² trefoil" tags=[:gauntlet, :pscad] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels
    # Catalogue construction is independent from native availability. Live campaigns are explicit CLI actions.
    benchmark=benchmark_definition(:benchmark_18kv_1000mm2_trefoil_pscad)
    @test benchmark isa BenchmarkDefinition
    @test benchmark.reference.id != benchmark.candidate.id
    @test benchmark.case_id === :cable_18kv_1000mm2_trefoil
    @test benchmark.reference.problem === benchmark.candidate.problem
end
