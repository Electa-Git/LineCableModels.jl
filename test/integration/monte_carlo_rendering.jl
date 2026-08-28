@testitem "PlotBuilder / Cairo / native Monte Carlo verbs" tags=[:visual] setup=[
    PlotBuilderTestSupport, UsePlotBuilderSupport, TestFixtures] begin
    if get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false")!="true"
        error("set LINECABLEMODELS_TEST_PLOTTING=true to run the visual contract")
    else
        using CairoMakie

        result=TestFixtures.cable_monte_carlo_result()
        calls=(
            Makie.hist(result, R;
                bins = 2, normalization = :pdf,
                backend = :cairo, display_plot = false),
            Makie.stairs(result, R;
                backend = :cairo, display_plot = false),
            Makie.ecdfplot(result, R;
                backend = :cairo, display_plot = false),
            Makie.lines(result, R;
                backend = :cairo, display_plot = false),
            Makie.qqplot(result, R;
                backend = :cairo, display_plot = false)
        )

        @test all(plot -> plot isa UIPlot, calls)
        @test all(plot -> length(plot.panels) == 1, calls)
        @test all(plot -> plot.context.backend === :cairo, calls)
        @test all(plot -> haskey(plot.controls, :reset), calls)
        @test all(plot -> haskey(plot.controls, :export_svg), calls)
        @test all(plot -> !haskey(plot.controls, :xlog), calls)
        @test all(plot -> !haskey(plot.controls, :ylog), calls)
        @test first(calls).panels[1].metadata.series[1].xdata ==
              [1_000.0, 2_000.0, 3_000.0, 4_000.0]
        @test length(calls[4].panels[1].metadata.series[1].xdata) == 500
        @test length(calls[5].panels[1].metadata.series) == 2

        sample_only=MonteCarloResult(
            result.formulation,
            result.values,
            result.stats,
            result.sample_values,
            nothing,
            result.root_seed,
            result.point_seeds,
            result.trial_counts
        )
        derived=Makie.stairs(sample_only, R;
            bins = 2, backend = :cairo, display_plot = false)
        @test derived isa UIPlot

        histogram_only=MonteCarloResult(
            result.formulation,
            result.values,
            result.stats,
            nothing,
            result.histogram_values,
            result.root_seed,
            result.point_seeds,
            result.trial_counts
        )
        @test_throws ArgumentError Makie.hist(
            histogram_only, R; backend = :cairo, display_plot = false)
        @test_throws ArgumentError Makie.qqplot(
            result, R; qqline = :fit, backend = :cairo, display_plot = false)

        multi_point=MonteCarloResult(
            result.formulation,
            repeat(result.values, 2),
            repeat(result.stats, 2),
            repeat(result.sample_values, 2),
            repeat(result.histogram_values, 2),
            result.root_seed,
            repeat(result.point_seeds, 2),
            repeat(result.trial_counts, 2)
        )
        @test_throws ArgumentError Makie.hist(
            multi_point, R; backend = :cairo, display_plot = false)
    end
end
