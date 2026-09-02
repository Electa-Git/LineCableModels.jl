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
