@testitem "Makie addons / native Monte Carlo verbs" tags=[:visual] setup=[
    NativePlotTestSupport, UseNativePlotSupport, TestFixtures
] begin
    get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false")=="true"||
    error("set LINECABLEMODELS_TEST_PLOTTING=true to run the visual contract")
    using CairoMakie

    result=TestFixtures.cable_monte_carlo_result()
    plots=(
        Makie.hist(
            result,
            R;
            backend = :cairo,
            display_plot = false,
            figure_title = "Retained samples",
            panel_titles = ("Sample histogram",),
            legend_title = "Distribution",
            legend_labels = ("observations",),
            legend_position = :inside,
            legend_anchor = :lt
        ),
        Makie.stairs(result, R; backend = :cairo, display_plot = false),
        Makie.ecdfplot(result, R; backend = :cairo, display_plot = false),
        Makie.lines(result, R; backend = :cairo, display_plot = false),
        Makie.qqplot(result, R; backend = :cairo, display_plot = false)
    )
    @test all(plot -> plot isa UIPlot, plots)
    @test all(plot -> length(plot.axes) == 1, plots)
    @test all(
        plot -> Set(keys(plot.controls)) ==
                Set((:reset, :export_svg, :legend)), plots)
    @test all(plot -> plot.legend !== nothing, plots)
    @test first(plots).title.text[] == "Retained samples"
    @test only(first(plots).axes).title[] == "Sample histogram"
    @test first(last(first(first(plots).legend.entrygroups[]))).label[] ==
          "observations"
    @test length(last(plots).legend.entrygroups[][1][2]) == 2
    @test_throws ArgumentError Makie.qqplot(
        result, R; qqline = :invalid, backend = :cairo, display_plot = false
    )
end

@testitem "Makie addons / histogram bin requests reach the retained marginal" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    result = TestFixtures.cable_monte_carlo_result()
    retained = observe(result, histograms, R, 1, 1)
    expected = HistogramDensity(observe(result, samples, R, 1, 1, :); bins=3)
    plot = Makie.stairs(result, R; bins=3, backend=:cairo, display_plot=false,
        controls=false, length_unit=:base, quantity_units=:base, clip=false)
    native = only(filter(item -> item isa Makie.Stairs, only(plot.axes).scene.plots))
    @test first.(native[1][]) ≈ expected.edges
    @test last.(native[1][]) ≈ [expected.density; last(expected.density)]
    @test length(retained.density) == 2
    @test retained.edges == [1.0, 3.0, 5.0]
    @test_throws ArgumentError Makie.stairs(result, R; bins=0, backend=:cairo,
        display_plot=false)
    without_samples = MonteCarloResult(result.formulation, result.values, result.stats,
        nothing, result.histogram_values, result.root_seed, result.point_seeds,
        result.trial_counts)
    @test_throws r"return_samples=true" Makie.stairs(without_samples, R; bins=3,
        backend=:cairo, display_plot=false)
end
