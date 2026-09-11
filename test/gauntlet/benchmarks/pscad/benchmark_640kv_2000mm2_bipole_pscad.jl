@testitem "PSCAD benchmark / 640 kV 2000 mm² compact-stranded bipole" tags=[:gauntlet, :pscad] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels
    # Catalogue construction is independent from native availability. Live campaigns are explicit CLI actions.
    benchmark=benchmark_definition(:benchmark_640kv_2000mm2_bipole_pscad)
    @test benchmark isa BenchmarkDefinition
    @test benchmark.reference.id != benchmark.candidate.id
    @test benchmark.case_id === :cable_640kv_2000mm2_bipole
    @test benchmark.reference.problem === benchmark.candidate.problem
end
