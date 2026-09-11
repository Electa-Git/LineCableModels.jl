@testitem "Gauntlet / clearance / NA2XS2Y touching trefoil analytical and UQ regression" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    using LineCableModels
    using Measurements
    model = Gauntlet.load_case(:cable_30kv_na2xs2y_630mm2_trefoil)
    problem = model.problem
    @test length(problem.frequencies) == 101
    @test first(problem.frequencies) == 0.1
    @test last(problem.frequencies) == 1e7
    inner = Gauntlet.uq_inner_formulation()
    result = compute(problem, inner)
    @test size(observe(result, Z)) == (6, 6, 101)
    @test all(isfinite, observe(result, Z))
    @test all(isfinite, observe(result, Y))

    benchmark = Gauntlet.benchmark_definition(
        :benchmark_30kv_na2xs2y_630mm2_trefoil_lep_montecarlo;
        frequencies = [0.1, 50.0, 1e7])
    uncertain_problem = benchmark.reference.problem
    mc = compute(uncertain_problem, MonteCarlo(inner; trials = 3, seed = 100,
        options = (retain_details = true,)))
    @test mc.trial_counts == [3]
    @test isempty(only(mc.details.failures))
    @test only(mc.details.clearance).adjustments == 3
    lep = compute(uncertain_problem, LinearError(inner))
    for quantity in (Z, Y)
        values = observe(only(lep.values), quantity)
        @test all(isfinite, nominal.(values))
        @test all(isfinite, uncertainty.(real.(values)))
        @test all(isfinite, uncertainty.(imag.(values)))
    end
end
