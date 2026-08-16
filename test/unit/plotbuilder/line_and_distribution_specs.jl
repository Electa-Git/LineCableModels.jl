@testitem "Engine / plot specification / line selection and physical units" tags=[:unit] setup=[
    PlotBuilderTestSupport,
    UsePlotBuilderSupport,
    TestFixtures
] begin
    const E=LineCableModels.Engine
    const PB=LineCableModels.PlotBuilder

    parameters=TestFixtures.two_conductor_results()
    series=Z(parameters)
    shunt=Y(parameters)

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
        E.LineParameterPlotSpec,
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
        E.LineParameterPlotSpec,
        parameters;
        quantities = [R]
    )
    @test_throws ArgumentError PB.make_render(
        E.LineParameterPlotSpec,
        parameters;
        quantities = (Z, real)
    )
    @test_throws DimensionMismatch PB.make_render(
        E.LineParameterPlotSpec,
        series;
        frequencies = [50.0]
    )
    @test_throws DomainError PB.make_render(
        E.LineParameterPlotSpec,
        series;
        frequencies = [0.0, 50.0],
        quantities = (L,)
    )

    residual_conductance=fill(1.0e-17, size(shunt))
    lossless=LineParameters(
        copy(Z(parameters).values),
        complex.(residual_conductance, imag.(Y(parameters).values)),
        frequencies(parameters)
    )
    lossless_render=PB.make_render(
        E.LineParameterPlotSpec,
        lossless;
        quantities = (G,)
    )
    lossless_series=only(only(lossless_render.figures).views).series
    @test all(series_spec -> all(iszero, series_spec.ydata), lossless_series)
    @test all(==(1.0e-17), G(lossless))

    small_conductance=fill(1.0e-12, size(shunt))
    lossy=LineParameters(
        copy(Z(parameters).values),
        complex.(small_conductance, imag.(Y(parameters).values)),
        frequencies(parameters)
    )
    lossy_render=PB.make_render(E.LineParameterPlotSpec, lossy; quantities = (G,))
    lossy_series=only(only(lossy_render.figures).views).series
    @test all(
        series_spec -> all(==(1.0e-9), series_spec.ydata),
        lossy_series
    )
    @test all(==(1.0e-12), G(lossy))
end

@testitem "Computation / plot specification / empirical and model distributions" tags=[:unit] setup=[
    PlotBuilderTestSupport,
    UsePlotBuilderSupport,
    TestFixtures
] begin
    const Cmp=LineCableModels.Computation
    const PB=LineCableModels.PlotBuilder
    const Spec=Cmp.MCDistributionPlotSpec

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

    histogram=result.histograms.R
    @test Cmp._mc_histogram_cdf(histogram, first(histogram.edges) - 1) == 0.0
    @test Cmp._mc_histogram_cdf(histogram, last(histogram.edges) + 1) == 1.0
    @test Cmp._mc_histogram_cdf(histogram, 2.0) ≈ 0.25
    @test Cmp._mc_histogram_quantile(histogram, 0.0) == first(histogram.edges)
    @test Cmp._mc_histogram_quantile(histogram, 1.0) == last(histogram.edges)
    @test Cmp._mc_histogram_quantile(histogram, 0.5) ≈ 3.0
    @test_throws DomainError Cmp._mc_histogram_quantile(histogram, -eps())

    samples_only=Cmp.MonteCarloResult(
        result.representation,
        result.statistics,
        result.samples,
        nothing,
        nothing,
        result.trials,
        result.confidence,
        result.cdf_tol,
        result.distribution,
        result.seed,
        result.manifest
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
        (Cmp.MCSeriesKey{:samples}(), :histogram, "samples"),
        (Cmp.MCSeriesKey{:histogram_pdf}(), :stairs, "model PDF"),
        (Cmp.MCSeriesKey{:histogram_cdf}(), :line, "model CDF"),
        (Cmp.MCSeriesKey{:empirical_cdf}(), :line, "empirical CDF"),
        (Cmp.MCSeriesKey{:quantiles}(), :scatter, "quantiles"),
        (Cmp.MCSeriesKey{:reference}(), :line, "perfect fit")
    )
    for (key, kind, label) in series_contracts
        @test PB.plot_kind(Spec, Val(:hist), samples_recipe, nothing, nothing, key) ===
              kind
        @test PB.legend_label(Spec, Val(:hist), samples_recipe, nothing, nothing, key) ==
              label
    end
    @test Cmp._mc_quantity_label(
        LineCableModels.UnitHandler.QuantityTag{:dimensionless}(),
        LineCableModels.UnitHandler.Units()
    ) == "dimensionless"

    histogram_only=Cmp.MonteCarloResult(
        result.representation,
        result.statistics,
        nothing,
        result.histograms,
        nothing,
        result.trials,
        result.confidence,
        result.cdf_tol,
        result.distribution,
        result.seed,
        result.manifest
    )
    @test PB.make_render(Spec, histogram_only; mode = :ecdf, data = :pdf) isa
          PB.RenderSpec
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
