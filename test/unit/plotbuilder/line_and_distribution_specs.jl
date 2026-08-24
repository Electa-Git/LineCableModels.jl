@testitem "Engine / plot specification / line selection and physical units" tags=[:unit] setup=[
    PlotBuilderTestSupport,
    UsePlotBuilderSupport,
    TestFixtures
] begin
    const E=LineCableModels.Engine
    const PB=LineCableModels.PlotBuilder

    parameters=TestFixtures.two_conductor_results()
    series=parameters.Z
    shunt=parameters.Y

    @test E._indices(nothing, 3) == [1, 2, 3]
    @test E._indices(:, 3) == [1, 2, 3]
    @test E._indices(2, 3) == [2]
    @test E._indices(1:2, 3) == [1, 2]
    @test E._indices([3, 1], 3) == [3, 1]
    @test_throws ArgumentError E._indices((1, 2), 3)
    @test E._conductor_pairs(series, (1, [1, 2])) == [(1, 1), (1, 2)]
    @test_throws ArgumentError E._conductor_pairs(series, 1)
    @test_throws BoundsError E._conductor_pairs(series, ([0], :))

    @test E._line_components(series, R) == (:R,)
    @test E._line_components(series, X) == (:X,)
    @test E._line_components(series, L) == (:L,)
    @test E._line_components(series, real) == (:Z_re,)
    @test E._line_components(series, imag) == (:Z_im,)
    @test E._line_components(series, abs) == (:Z_abs,)
    @test E._line_components(series, angle) == (:Z_angle,)
    @test E._line_components(series, Z) == (:Z_re, :Z_im)
    @test_throws ArgumentError E._line_components(series, C)

    @test E._line_components(shunt, G) == (:G,)
    @test E._line_components(shunt, B) == (:B,)
    @test E._line_components(shunt, C) == (:C,)
    @test E._line_components(shunt, real) == (:Y_re,)
    @test E._line_components(shunt, imag) == (:Y_im,)
    @test E._line_components(shunt, abs) == (:Y_abs,)
    @test E._line_components(shunt, angle) == (:Y_angle,)
    @test E._line_components(shunt, Y) == (:Y_re, :Y_im)
    @test_throws ArgumentError E._line_components(shunt, R)

    @test E._line_components(parameters, real) == (:Z_re, :Y_re)
    @test E._line_components(parameters, imag) == (:Z_im, :Y_im)
    @test E._line_components(parameters, abs) == (:Z_abs, :Y_abs)
    @test E._line_components(parameters, angle) == (:Z_angle, :Y_angle)
    @test_throws ArgumentError E._line_components(parameters, identity)

    render=PB.make_render(
        E.LineParameterPlotDefinition,
        parameters;
        quantities = (R, L, G, C),
        con = (1:2, [1, 2])
    )
    @test length(render.figures) == 2
    @test sort([page.title for page in render.figures]) ==
          ["Series impedance", "Shunt admittance"]
    @test sum(length(view.series) for page in render.figures for view in page.views) == 16
    @test all(
        item -> item.xaxis.allowed_scales == (:linear, :log10),
        (view for page in render.figures for view in page.views)
    )
    @test Set(
        series_spec.label
    for page in render.figures
    for view in page.views
    for series_spec in view.series
    ) == Set(["Z[1,1]", "Z[1,2]", "Z[2,1]", "Z[2,2]",
        "Y[1,1]", "Y[1,2]", "Y[2,1]", "Y[2,2]"])

    @test_throws ArgumentError PB.make_render(
        E.LineParameterPlotDefinition,
        parameters;
        quantities = [R]
    )
    @test_throws ArgumentError PB.make_render(
        E.LineParameterPlotDefinition,
        parameters;
        quantities = (Z, real)
    )
    @test_throws DimensionMismatch PB.make_render(
        E.LineParameterPlotDefinition,
        series;
        frequencies = [50.0]
    )
    @test_throws DomainError PB.make_render(
        E.LineParameterPlotDefinition,
        series;
        frequencies = [0.0, 50.0],
        quantities = (L,)
    )

    residual_conductance=fill(1.0e-17, size(shunt))
    lossless=LineParameters(
        copy(Z(parameters)),
        complex.(residual_conductance, imag.(Y(parameters))),
        frequencies(parameters)
    )
    lossless_render=PB.make_render(
        E.LineParameterPlotDefinition,
        lossless;
        quantities = (G,)
    )
    lossless_series=only(only(lossless_render.figures).views).series
    @test all(series_spec -> all(iszero, series_spec.ydata), lossless_series)
    @test all(==(1.0e-17), G(lossless))

    small_conductance=fill(1.0e-12, size(shunt))
    lossy=LineParameters(
        copy(Z(parameters)),
        complex.(small_conductance, imag.(Y(parameters))),
        frequencies(parameters)
    )
    lossy_render=PB.make_render(E.LineParameterPlotDefinition, lossy; quantities = (G,))
    lossy_series=only(only(lossy_render.figures).views).series
    @test all(
        series_spec -> all(==(1.0e-9), series_spec.ydata),
        lossy_series
    )
    @test all(==(1.0e-12), G(lossy))
end

@testitem "Engine / plot specification / line-parameter comparison grid" tags=[:unit] setup=[
    PlotBuilderTestSupport,
    UsePlotBuilderSupport
] begin
    const E=LineCableModels.Engine
    const PB=LineCableModels.PlotBuilder

    frequency=[10.0, 100.0, 1_000.0]
    resistance=[1.0 0.2 0.1; 0.2 1.2 0.15; 0.1 0.15 1.4] .* 1.0e-4
    inductance=[2.0 0.4 0.3; 0.4 2.2 0.35; 0.3 0.35 2.4] .* 1.0e-7
    conductance=[3.0 0.5 0.4; 0.5 3.2 0.45; 0.4 0.45 3.4] .* 1.0e-9
    capacitance=[4.0 0.8 0.6; 0.8 4.2 0.7; 0.6 0.7 4.4] .* 1.0e-10
    function result(scale; result_frequency = frequency, result_domain = PhaseDomain,
            result_basis = :per_length)
        impedance=repeat(scale .* resistance, 1, 1, length(result_frequency)) .+
                  im .* repeat(scale .* inductance, 1, 1, length(result_frequency)) .*
                  reshape(2π .* result_frequency, 1, 1, :)
        admittance=repeat(scale .* conductance, 1, 1, length(result_frequency)) .+
                   im .* repeat(scale .* capacitance, 1, 1, length(result_frequency)) .*
                   reshape(2π .* result_frequency, 1, 1, :)
        return LineParameters(
            result_domain,
            impedance,
            admittance,
            result_frequency;
            basis = result_basis
        )
    end

    parameters=Tuple(result(1+0.01index) for index in 0:4)
    labels=Tuple("Formulation $index" for index in eachindex(parameters))
    render=PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        parameters;
        legend = labels,
        quantities = (R, C),
        xscale = :log10
    )

    @test length(render.figures) == 2
    @test all(page -> page.size == (1200, 800), render.figures)
    @test all(page -> length(page.views) == 9, render.figures)
    @test first(first(render.figures).views).title ==
          "Z[1,1] · Series resistance"
    @test all(
        view -> length(view.series) == length(parameters),
        (view for page in render.figures for view in page.views)
    )
    @test Set(
        (first(view.placement.area.rows), first(view.placement.area.columns))
    for view in first(render.figures).views
    ) == Set((row, column) for row in 1:3 for column in 1:3)
    @test Set(
        series.label
    for view in first(render.figures).views
    for series in view.series
    ) == Set(labels)
    @test Set(
        series.group
    for view in first(render.figures).views
    for series in view.series
    ) == Set(Symbol("line_parameters_$index") for index in eachindex(parameters))
    @test all(
        series -> series.attributes.linestyle === :solid,
        (series for page in render.figures for view in page.views
        for series in view.series)
    )
    first_panel=first(first(render.figures).views)
    @test length(unique(series.attributes.color for series in first_panel.series)) ==
          length(parameters)
    @test all(
        page -> [series.attributes.color for series in first(page.views).series] ==
                [series.attributes.color for series in first_panel.series],
        render.figures
    )

    page=first(render.figures)
    root=only(page.layout.grids)
    @test page.layout.name === :line_parameters_comparison
    @test root.rows isa Vector{PB.AbstractTrackSize}
    @test root.rows[2] isa PB.RelativeTrack
    @test root.columns[1] isa PB.RelativeTrack
    @test root.columns[2] isa PB.ContentTrack
    @test Set(slot.name for slot in page.layout.slots) ==
          Set((:toolbar, :canvas, :legend, :status))
    @test page.controls.reset
    @test page.controls.export_svg
    @test page.legend.interactive
    @test page.legend.overflow === :ellipsis
    @test page.status.initial == "Ready."

    default_render=PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        parameters[1:2];
        legend = labels[1:2]
    )
    @test getfield.(getproperty.(default_render.figures, :key), :component) ==
          [:Z_re, :Z_im, :Y_re, :Y_im]

    @test_throws ArgumentError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        (first(parameters),);
        legend = ("one",)
    )
    @test_throws ArgumentError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        parameters
    )
    @test_throws ArgumentError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        parameters;
        legend = collect(labels)
    )
    @test_throws DimensionMismatch PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        parameters;
        legend = labels[1:2]
    )
    smaller=LineParameters(
        PhaseDomain,
        Z(first(parameters))[1:2, 1:2, :],
        Y(first(parameters))[1:2, 1:2, :],
        frequency
    )
    @test_throws DimensionMismatch PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        (first(parameters), smaller);
        legend = ("three conductors", "two conductors")
    )
    @test_throws ArgumentError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        (parameters[1], result(1.0; result_frequency = [10.0, 100.0, 2_000.0]));
        legend = ("first", "frequency mismatch")
    )
    @test_throws ArgumentError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        (parameters[1], result(1.0; result_domain = ModalDomain));
        legend = ("phase", "modal")
    )
    @test_throws ArgumentError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        (parameters[1], result(1.0; result_basis = :total));
        legend = ("per length", "total")
    )
    @test_throws DomainError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        (
            result(1.0; result_frequency = [0.0, 10.0, 100.0]),
            result(1.1; result_frequency = [0.0, 10.0, 100.0])
        );
        legend = ("first", "second"),
        quantities = (L,)
    )
end

@testitem "UQ / plot specification / empirical and model distributions" tags=[:unit] setup=[
    PlotBuilderTestSupport,
    UsePlotBuilderSupport,
    TestFixtures
] begin
    const Cmp=LineCableModels.UQ
    const PB=LineCableModels.PlotBuilder
    const Spec=Cmp.MCDistributionPlotDefinition

    result=TestFixtures.cable_monte_carlo_result()
    expected_kinds=Dict(
        (:hist, :samples)=>[:histogram],
        (:hist, :pdf)=>[:stairs],
        (:hist, :both)=>[:histogram, :stairs],
        (:pdf, :samples)=>[:stairs],
        (:ecdf, :samples)=>[:line],
        (:ecdf, :pdf)=>[:line],
        (:ecdf, :both)=>[:line, :line],
        (:qq, :samples)=>[:scatter, :line]
    )
    for ((mode, data), kinds) in expected_kinds
        render=PB.make_render(Spec, result; mode, data)
        view=only(only(render.figures).views)
        @test getfield.(view.series, :kind) == kinds
        for series_spec in view.series
            series_spec.xdata===nothing||@test all(isfinite, series_spec.xdata)
            series_spec.ydata===nothing||@test all(isfinite, series_spec.ydata)
        end
    end

    histogram=only(Cmp.histograms(result)).R
    @test Cmp._mc_histogram_cdf(histogram, first(histogram.edges) - 1) == 0.0
    @test Cmp._mc_histogram_cdf(histogram, last(histogram.edges) + 1) == 1.0
    @test Cmp._mc_histogram_cdf(histogram, 2.0) ≈ 0.25
    @test Cmp._mc_histogram_quantile(histogram, 0.0) == first(histogram.edges)
    @test Cmp._mc_histogram_quantile(histogram, 1.0) == last(histogram.edges)
    @test Cmp._mc_histogram_quantile(histogram, 0.5) ≈ 3.0
    @test_throws DomainError Cmp._mc_histogram_quantile(histogram, -eps())

    samples_only=Cmp.MonteCarloResult(
        result.formulation,
        result.values,
        result.stats,
        result.sample_values,
        nothing,
        result.root_seed,
        result.point_seeds,
        result.trial_counts
    )
    derived=PB.make_render(Spec, samples_only; mode = :pdf, nbins = 2)
    derived_series=only(only(only(derived.figures).views).series)
    @test derived_series.kind === :stairs
    @test length(derived_series.xdata) == 3

    samples_recipe=PB.resolve_input(
        Spec,
        PB.parse_kwargs(Spec, samples_only; mode = :hist, data = :samples)
    )
    @test samples_recipe.input.histogram === nothing
    @test first(samples_recipe.input.bins) == minimum(samples_recipe.input.values)
    @test last(samples_recipe.input.bins) == maximum(samples_recipe.input.values)
    @test length(samples_recipe.input.bins) == 3
    cdf_grid=Cmp._mc_cdf_grid(samples_recipe)
    @test length(cdf_grid) == 500
    @test first(cdf_grid) < minimum(samples_recipe.input.values)
    @test last(cdf_grid) > maximum(samples_recipe.input.values)

    series_contracts=(
        (Val(:samples), :histogram, "samples"),
        (Val(:histogram_pdf), :stairs, "model PDF"),
        (Val(:histogram_cdf), :line, "model CDF"),
        (Val(:empirical_cdf), :line, "empirical CDF"),
        (Val(:quantiles), :scatter, "quantiles"),
        (Val(:reference), :line, "perfect fit")
    )
    for (key, kind, label) in series_contracts
        @test PB.plot_kind(Spec, Val(:hist), samples_recipe, nothing, nothing, key) ===
              kind
        @test PB.legend_label(Spec, Val(:hist), samples_recipe, nothing, nothing, key) ==
              label
    end
    @test Cmp._mc_quantity_label(
        LineCableModels.Units.QuantityTag{:dimensionless}(),
        LineCableModels.Units.UnitExpr()
    ) == "Dimensionless"

    histogram_only=Cmp.MonteCarloResult(
        result.formulation,
        result.values,
        result.stats,
        nothing,
        result.histogram_values,
        result.root_seed,
        result.point_seeds,
        result.trial_counts
    )
    @test PB.make_render(Spec, histogram_only; mode = :ecdf, data = :pdf) isa
          PB.PlotRecipe
    histogram_recipe=PB.resolve_input(
        Spec,
        PB.parse_kwargs(Spec, histogram_only; mode = :ecdf, data = :pdf)
    )
    @test_throws ArgumentError Cmp._mc_values(histogram_recipe)
    @test_throws ArgumentError PB.make_render(
        Spec,
        histogram_only;
        mode = :ecdf,
        data = :samples
    )
    @test_throws ArgumentError PB.make_render(Spec, histogram_only; mode = :qq)
end
