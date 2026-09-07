@testitem "UQ / histogram observations / explicit binning and retention" tags=[:unit] setup=[TestFixtures] begin
    cable = TestFixtures.cable_monte_carlo_result()
    frequency = [50.0, 1000.0]
    scale = reshape(collect(1.0:8.0), 2, 2, 2)
    sample_arrays = (R=cat((scale .* trial .* 1e-4 for trial in 1:4)...; dims=4),
        L=cat((scale .* trial .* 1e-7 for trial in 1:4)...; dims=4),
        C=cat((scale .* trial .* 1e-10 for trial in 1:4)...; dims=4),
        G=cat((scale .* trial .* 1e-9 for trial in 1:4)...; dims=4))
    summaries = map(values -> map(CartesianIndices(scale)) do index
        SampleSummary(values[Tuple(index)..., :])
    end, sample_arrays)
    models = map(values -> map(CartesianIndices(scale)) do index
        HistogramDensity(values[Tuple(index)..., :]; bins=2)
    end, sample_arrays)
    omega = reshape(2pi .* frequency, 1, 1, :)
    line = LineParameters(2.5e-4 .* scale .+ im .* omega .* 2.5e-7 .* scale,
        2.5e-9 .* scale .+ im .* omega .* 2.5e-10 .* scale, frequency)
    matrix = MonteCarloResult(cable.formulation, [line], [summaries], [sample_arrays],
        [models], cable.root_seed, cable.point_seeds, cable.trial_counts)

    for (source, indices) in ((cable, (1,)), (matrix, (2, 1, 2)))
        samples_only = MonteCarloResult(source.formulation, source.values, source.stats,
            source.sample_values, nothing, source.root_seed, source.point_seeds,
            source.trial_counts)
        model_only = MonteCarloResult(source.formulation, source.values, source.stats,
            nothing, source.histogram_values, source.root_seed, source.point_seeds,
            source.trial_counts)
        summaries_only = MonteCarloResult(source.formulation, source.values, source.stats,
            nothing, nothing, source.root_seed, source.point_seeds, source.trial_counts)
        for selector in (R, L, C, G)
            retained = observe(source, histograms, selector, 1, indices...)
            sample = observe(source, samples, selector, 1, indices..., :)
            original_edges, original_density = copy(retained.edges), copy(retained.density)
            @test @inferred(observe(source, histograms, selector, 1, indices..., nothing)) === retained
            @test observe(model_only, histograms, selector, 1, indices..., nothing) === retained
            @test observe(model_only, histograms, selector, 1, indices..., 2) === retained
            for input in (source, samples_only), bins in (1, 3, 4)
                model = @inferred observe(input, histograms, selector, 1, indices..., bins)
                expected = HistogramDensity(sample; bins)
                @test length(model.density) == bins
                @test model.edges == expected.edges
                @test model.density == expected.density
                @test sum(model.density .* diff(model.edges)) ≈ 1.0
            end
            automatic = observe(samples_only, histograms, selector, 1, indices..., nothing)
            @test automatic.edges == HistogramDensity(sample).edges
            @test automatic.density == HistogramDensity(sample).density
            @test_throws r"return_samples=true" observe(model_only, histograms, selector, 1, indices..., 3)
            @test_throws ArgumentError observe(summaries_only, histograms, selector, 1, indices..., nothing)
            for input in (source, model_only, samples_only), bins in (0, -1)
                @test_throws ArgumentError observe(input, histograms, selector, 1, indices..., bins)
            end
            @test retained.edges == original_edges
            @test retained.density == original_density
            @test observe(source, samples, selector, 1, indices..., :) == sample
        end
    end

    constant_samples = map(_ -> fill(2.0, 1, 4), only(cable.sample_values))
    constant_stats = map(_ -> [SampleSummary(fill(2.0, 4))], only(cable.stats))
    constant = MonteCarloResult(cable.formulation, [CableConstants(2.0, 2.0, 2.0, 2.0)],
        [constant_stats], [constant_samples], nothing, cable.root_seed,
        cable.point_seeds, cable.trial_counts)
    for bins in (nothing, 1, 4)
        model = observe(constant, histograms, R, 1, 1, bins)
        @test length(model.density) == 1
        @test only(diff(model.edges)) > 0
        @test first(model.edges) < 2 < last(model.edges)
        @test sum(model.density .* diff(model.edges)) ≈ 1.0
    end
end
