@testitem "PSCAD benchmark / solid 1000 mm² single phase" tags=[:gauntlet, :pscad] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels
    # Catalogue construction is independent from native availability. Live campaigns are explicit CLI actions.
    benchmark=benchmark_definition(:benchmark_solid_1000mm2_single_pscad)
    @test benchmark isa BenchmarkDefinition
    @test benchmark.reference.id != benchmark.candidate.id
    @test benchmark.case_id === :solid_1000mm2_single
    @test benchmark.reference.problem === benchmark.candidate.problem
end
