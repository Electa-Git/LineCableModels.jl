@testitem "UQ benchmark / two insulated wires / LEP versus Monte Carlo" tags=[:gauntlet, :uq] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels
    # Catalogue construction is independent from native availability. Live campaigns are explicit CLI actions.
    benchmark=benchmark_definition(:benchmark_two_insulated_wires_lep_montecarlo)
    @test benchmark isa BenchmarkDefinition
    @test benchmark.reference.id != benchmark.candidate.id
    @test benchmark.case_id === :two_insulated_wires
    @test benchmark.reference.problem === benchmark.candidate.problem
end
