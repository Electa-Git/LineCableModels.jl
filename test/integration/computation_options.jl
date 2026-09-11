@testitem "Execution options / inner owner validation through traversal" tags=[:integration] setup=[TestFixtures] begin
    inner = CableConstantsFormulation()
    problem = CableConstantsProblem(TestFixtures.mv_cable_design())
    space = Gridspace{CableConstantsProblem}(identity, (Grid((problem,)),))
    point = first(LineCableModels.points(space))
    invalid = (unsupported=true,)
    pending = ParametricProblem(space, invalid)
    @test pending.options === invalid
    @test_throws ArgumentError compute(point, inner; options=invalid)
    @test_throws ArgumentError compute(problem, [inner, inner]; options=invalid)
    formulas = Gridspace{CableConstantsFormulation}(identity, (Grid((inner,)),))
    @test_throws ArgumentError compute(space, formulas; options=invalid)
    for outer in (Combinatorial(inner), LinearError(inner),
            MonteCarlo(inner; options=(trials=1, seed=1)))
        @test_throws ArgumentError compute(pending, outer)
        @test isempty(details(compute(ParametricProblem(space), outer)))
    end
end

@testitem "Execution options / callbacks and output basis survive all traversal paths" tags=[:integration] setup=[TestFixtures] begin
    problem = TestFixtures.line_parameters_problem()
    inner = Formulation()
    space = Gridspace{LineParametersProblem}(identity, (Grid((problem,)),))
    calls = Int[]
    callback = (problem, index, result) -> push!(calls, index)
    options = (output_basis=:total, on_result=callback)
    # The selected basis is represented by Val, so a runtime Symbol determines
    # that field's type; the normalized record must still retain concrete fields.
    execution = LineCableModels.computation_options(LineCableModelsCoaxial, options)
    @test isconcretetype(typeof(execution))
    @test execution.on_result === callback
    @test fieldtype(typeof(execution), :on_result) === typeof(callback)
    batch = compute(problem, [inner, inner]; options)
    @test calls == [1, 2]
    @test all(value -> basis(value) === :total, batch)
    empty!(calls)
    value = compute(first(LineCableModels.points(space)), inner; options)
    @test calls == [1]
    @test basis(value) === :total
    for (outer, count) in ((Combinatorial(inner), 1), (LinearError(inner), 1),
            (MonteCarlo(inner; options=(trials=2, seed=1)), 2))
        empty!(calls)
        result = compute(ParametricProblem(space, options), outer)
        @test calls == fill(1, count)
        @test basis(first(result)) === :total
    end
end

@testitem "UQ / tuple and shorthand controls produce identical seeded calculations" tags=[:integration] setup=[TestFixtures] begin
    using Measurements
    using Distributions
    using Random
    design = TestFixtures.mv_cable_design()
    inner = CableConstantsFormulation()
    space = Gridspace{CableConstantsProblem}(
        temperature -> CableConstantsProblem(design; temperature),
        (Grid((20.0, 40.0), AbsoluteError(0.5)),))
    problem = ParametricProblem(space)
    options = (trials=4, seed=71, distribution=Normal(10, 3),
        return_samples=true, return_histograms=true, bins=2, retain_details=true)
    formulation = @inferred MonteCarlo(inner; options)
    @test formulation.options.distribution === options.distribution
    sampled = compute(problem, formulation)
    replay = compute(problem, MonteCarlo(inner; options...))
    @test samples(sampled) == samples(replay)
    @test statistics(sampled) == statistics(replay)
    @test sampled.point_seeds == replay.point_seeds
    @test sampled.trial_counts == replay.trial_counts == [4, 4]
    @test length.(sampled.details.trials) == [4, 4]
    @test LineCableModels.sampling_distribution(sampled) === options.distribution
    for quantity in (R, L, C, G), point in 1:2
        a = only(observe(sampled, histograms, quantity, point))
        b = only(observe(replay, histograms, quantity, point))
        @test a.edges == b.edges
        @test a.density == b.density
        values = observe(sampled, samples, quantity, point)
        @test length(a.density) == (all(==(first(values)), values) ? 1 : 2)
    end

    sampler = (rng, center, deviation) -> center + deviation * randn(rng)
    custom = compute(problem, MonteCarlo(inner;
        options=merge(options, (distribution=sampler,))))
    normal = compute(problem, MonteCarlo(inner;
        options=merge(options, (distribution=:normal,))))
    @test samples(custom) == samples(normal)
    @test LineCableModels.sampling_distribution(custom) === sampler

    automatic = compute(problem, MonteCarlo(inner;
        options=(seed=71, confidence=0.8, cdf_tol=0.5)))
    # One concentric assembly contributes four scalar R/L/C/G marginals.
    expected = ceil(Int, log(2 * 4 / (1 - 0.8)) / (2 * 0.5^2))
    @test automatic.trial_counts == [expected, expected]
    @test LineCableModels.confidence(automatic) == 0.8
    @test LineCableModels.cdf_tolerance(automatic) == 0.5
    @test samples(automatic) === nothing
    @test histograms(automatic) === nothing
end
