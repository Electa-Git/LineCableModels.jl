@testitem "Engine / detached line pages / selection and physical units" tags=[:unit] setup=[
    PlotBuilderTestSupport,
    UsePlotBuilderSupport,
    TestFixtures
] begin
    const E=LineCableModels.Engine
    const PB=LineCableModels.PlotBuilder
    const U=LineCableModels.Units

    parameters=TestFixtures.two_conductor_results()
    series=parameters.Z
    shunt=parameters.Y

    @test U.family(U.quantity(R)) === Val(:series)
    @test U.family(U.quantity(L)) === Val(:series)
    @test U.family(U.quantity(Z, abs)) === Val(:series)
    @test U.family(U.quantity(G)) === Val(:shunt)
    @test U.family(U.quantity(C)) === Val(:shunt)
    @test U.family(U.quantity(Y, angle)) === Val(:shunt)
    @test Base.ispublic(U, :family)
    @test !isdefined(E, :line_requests)
    @test !isdefined(E, :line_parent)

    requests=(
        @observe(R[1:2, [1, 2], :]),
        @observe(L[1:2, [1, 2], :]),
        @observe(G[1:2, [1, 2], :]),
        @observe(C[1:2, [1, 2], :])
    )
    render=PB.make_render(
        E.LineParameterPlotDefinition,
        parameters;
        requests
    )
    @test length(render.pages) == 2
    @test sort([page.title for page in render.pages]) ==
          ["Series impedance", "Shunt admittance"]
    @test sum(length(rows) * length(columns)
    for page in render.pages
    for (rows, columns, _) in page.payload.coordinates) == 16
    @test all(
        scales -> scales == (:linear, :log10),
        (scales for page in render.pages for scales in page.payload.xscales)
    )
    @test all(
        keys(page.payload.frequency) == (:values, :quantity, :unit) &&
            all(observation -> keys(observation) == (:values, :quantity, :unit),
                page.payload.observations)
    for page in render.pages
    )
    @test Tuple(request for page in render.pages for request in page.payload.requests) ==
          requests

    @test_throws ArgumentError PB.make_render(
        E.LineParameterPlotDefinition,
        parameters;
        requests = (@observe(R[:, :, :]), @observe(R[1, 1, :]))
    )
    @test_throws ArgumentError PB.make_render(E.LineParameterPlotDefinition, parameters)
    @test_throws DimensionMismatch PB.make_render(
        E.LineParameterPlotDefinition,
        series;
        frequencies = [50.0],
        requests = (@observe(R[:, :, :]),)
    )
    standalone_frequencies=frequencies(parameters)
    @test_throws DomainError PB.make_render(
        E.LineParameterPlotDefinition,
        series;
        frequencies = [0.0, 50.0],
        requests = (@observe(L[:, :, :]),)
    )
    standalone_shunt=PB.make_render(
        E.LineParameterPlotDefinition,
        shunt;
        frequencies = standalone_frequencies,
        requests = (@observe(G[:, :, :]), @observe(C[:, :, :]))
    )
    @test length(standalone_shunt.pages) == 1
    @test length(only(standalone_shunt.pages).payload.positions) == 2

    residual_conductance=fill(1.0e-17, size(shunt))
    lossless=LineParameters(
        copy(Z(parameters)),
        complex.(residual_conductance, imag.(Y(parameters))),
        frequencies(parameters)
    )
    lossless_render=PB.make_render(
        E.LineParameterPlotDefinition,
        lossless;
        requests = (@observe(G[:, :, :]),)
    )
    lossless_values=only(only(lossless_render.pages).payload.observations).values
    @test all(==(1.0e-14), lossless_values)
    @test all(==(1.0e-17), G(lossless))

    small_conductance=fill(1.0e-12, size(shunt))
    lossy=LineParameters(
        copy(Z(parameters)),
        complex.(small_conductance, imag.(Y(parameters))),
        frequencies(parameters)
    )
    lossy_render=PB.make_render(
        E.LineParameterPlotDefinition,
        lossy;
        requests = (@observe(G[:, :, :]),)
    )
    lossy_values=only(only(lossy_render.pages).payload.observations).values
    @test all(==(1.0e-9), lossy_values)
    @test all(==(1.0e-12), G(lossy))

    milli_resistance=PB.make_render(
        E.LineParameterPlotDefinition,
        parameters;
        requests = (@observe(R[:, :, :]),),
        quantity_units = Dict(R=>:milli)
    )
    expected_target=only(LineCableModels.Grammar.unit_targets(
        (R,),
        basis(parameters);
        length_prefix = :kilo,
        overrides = Dict(R=>:milli)
    ))
    @test only(only(milli_resistance.pages).payload.observations).unit ==
          expected_target

    presented=PB.make_render(
        E.LineParameterPlotDefinition,
        parameters;
        requests = (@observe(R[:, :, :]), @observe(X[:, :, :])),
        title = "Custom title",
        labels = ("Measured resistance", "Measured reactance"),
        legend = ("self 1", "mutual 1–2", "mutual 2–1", "self 2")
    )
    @test only(presented.pages).title == "Custom title"
    @test only(presented.pages).payload.titles ==
          ("Measured resistance", "Measured reactance")
    @test only(presented.pages).payload.legend_labels ==
          ("self 1", "mutual 1–2", "mutual 2–1", "self 2")
end

@testitem "Engine / detached comparison pages / matrix grid" tags=[:unit] setup=[
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
            result_basis = :pul)
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
        requests = (@observe(R[:, :, :]), @observe(C[:, :, :])),
        xscale = :log10
    )

    @test length(render.pages) == 2
    @test all(page -> page.size == (1200, 800), render.pages)
    @test all(page -> length(page.payload.positions) == 9, render.pages)
    @test first(first(render.pages).payload.titles) ==
          "R[1,1] · Series resistance"
    @test Set(
        first(render.pages).payload.positions
    ) == Set((row, column) for row in 1:3 for column in 1:3)
    @test first(render.pages).payload.legend_labels == labels
    first_colors=first(render.pages).payload.colors
    @test length(unique(first_colors)) ==
          length(parameters)
    @test all(
        page -> page.payload.colors == first_colors,
        render.pages
    )

    page=first(render.pages)
    @test page.legend.interactive
    @test page.legend.overflow === :ellipsis
    @test all(
        name -> !hasproperty(page.payload, name),
        (:title, :key, :colorbars, :widgets, :export_definition)
    )

    default_render=PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        parameters[1:2];
        legend = labels[1:2],
        requests = (
            @observe(R[:, :, :]), @observe(X[:, :, :]),
            @observe(G[:, :, :]), @observe(B[:, :, :])
        )
    )
    @test Tuple(page.key.request for page in default_render.pages) ==
          (
        @observe(R[:, :, :]), @observe(X[:, :, :]),
        @observe(G[:, :, :]), @observe(B[:, :, :])
    )

    @test_throws ArgumentError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        (first(parameters),);
        legend = ("one",),
        requests = (@observe(R[:, :, :]),)
    )
    @test_throws ArgumentError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        parameters
    )
    @test_throws ArgumentError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        parameters;
        legend = collect(labels),
        requests = (@observe(R[:, :, :]),)
    )
    @test_throws DimensionMismatch PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        parameters;
        legend = labels[1:2],
        requests = (@observe(R[:, :, :]),)
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
        legend = ("three conductors", "two conductors"),
        requests = (@observe(R[:, :, :]),)
    )
    @test_throws ArgumentError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        (parameters[1], result(1.0; result_frequency = [10.0, 100.0, 2_000.0]));
        legend = ("first", "frequency mismatch"),
        requests = (@observe(R[:, :, :]),)
    )
    @test_throws ArgumentError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        (parameters[1], result(1.0; result_domain = ModalDomain));
        legend = ("phase", "modal"),
        requests = (@observe(R[:, :, :]),)
    )
    @test_throws ArgumentError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        (parameters[1], result(1.0; result_basis = :total));
        legend = ("per length", "total"),
        requests = (@observe(R[:, :, :]),)
    )
    @test_throws DomainError PB.make_render(
        E.LineParametersBenchmarkPlotDefinition,
        (
            result(1.0; result_frequency = [0.0, 10.0, 100.0]),
            result(1.1; result_frequency = [0.0, 10.0, 100.0])
        );
        legend = ("first", "second"),
        requests = (@observe(L[:, :, :]),)
    )
end

@testitem "UQ / detached distribution page / empirical and model data" tags=[:unit] setup=[
    PlotBuilderTestSupport,
    UsePlotBuilderSupport,
    TestFixtures
] begin
    const Cmp=LineCableModels.UQ
    const PB=LineCableModels.PlotBuilder
    const Spec=Cmp.MCDistributionPlotDefinition

    result=TestFixtures.cable_monte_carlo_result()
    @test_throws ArgumentError PB.make_render(Spec, result; selector = :R)
    @test_throws ArgumentError PB.make_render(Spec, result; quantity = :R)
    @test PB.make_render(Spec, result; selector = L) isa PB.PlotRecipe
    detached_page=only(PB.make_render(Spec, result).pages)
    @test all(
        name -> !hasproperty(detached_page.payload, name),
        (:title, :key, :legend, :colorbars, :widgets, :export_definition)
    )
    kind_value(::Val{kind}) where {kind} = kind
    expected_kinds=Dict(
        (:hist, :samples)=>(:histogram,),
        (:hist, :pdf)=>(:stairs,),
        (:hist, :both)=>(:histogram, :stairs),
        (:pdf, :samples)=>(:stairs,),
        (:ecdf, :samples)=>(:line,),
        (:ecdf, :pdf)=>(:line,),
        (:ecdf, :both)=>(:line, :line),
        (:qq, :samples)=>(:scatter, :line)
    )
    for ((mode, data), kinds) in expected_kinds
        render=PB.make_render(Spec, result; mode, data)
        page=only(render.pages).payload
        @test page isa Cmp.MCDistributionPagePayload
        @test kind_value.(getfield.(page.layers, :kind)) == kinds
        @test keys(page.x_observation) == (:values, :quantity, :unit)
        @test keys(page.y_observation) == (:values, :quantity, :unit)
        for layer in page.layers
            layer.x===nothing||@test all(isfinite, layer.x)
            layer.y===nothing||@test all(isfinite, layer.y)
        end
    end

    pdf_page=only(PB.make_render(Spec, result; mode = :pdf).pages).payload
    @test pdf_page.y_observation.quantity isa
          LineCableModels.Units.Quantity{:probability_density}
    @test pdf_page.y_observation.unit == inv(pdf_page.x_observation.unit)
    @test LineCableModels.Units.label(
        pdf_page.y_observation.quantity,
        pdf_page.y_observation.unit
    ) == "Probability density [km/Ω]"
    milli_page=only(PB.make_render(
        Spec,
        result;
        selector = R,
        quantity_units = :milli
    ).pages).payload
    expected_target=only(LineCableModels.Grammar.unit_targets(
        (R,),
        basis(only(Cmp.statistics(result)));
        length_prefix = :kilo,
        overrides = :milli
    ))
    @test milli_page.x_observation.unit == expected_target
    @test !isdefined(Cmp, :_mc_plot_exponent)

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
    derived_layer=only(only(derived.pages).payload.layers)
    @test derived_layer.kind == Val(:stairs)
    @test length(derived_layer.x) == 3

    parsed=PB.parse(Spec, samples_only; mode = :hist, data = :samples)
    resolved=PB.resolve(Spec, samples_only, parsed)
    sample_pages=PB.fetch(Spec, samples_only, resolved)
    samples_layer=only(only(sample_pages).payload.layers)
    @test samples_layer.kind == Val(:histogram)
    @test first(samples_layer.style.bins) == minimum(samples_layer.x)
    @test last(samples_layer.style.bins) == maximum(samples_layer.x)
    @test length(samples_layer.style.bins) == 3
    @test !hasproperty(only(sample_pages).payload, :histogram)
    @test !hasproperty(only(sample_pages).payload, :values)
    cdf_grid=Cmp._mc_cdf_grid(samples_layer.x, nothing)
    @test length(cdf_grid) == 500
    @test first(cdf_grid) < minimum(samples_layer.x)
    @test last(cdf_grid) > maximum(samples_layer.x)

    combined=only(PB.make_render(
        Spec,
        result;
        mode = :ecdf,
        data = :both
    ).pages).payload
    @test getfield.(combined.layers, :group) == (:empirical_cdf, :histogram_cdf)
    @test getfield.(combined.layers, :label) == ("empirical CDF", "model CDF")
    @test LineCableModels.Units.label(
        LineCableModels.Units.Quantity{:dimensionless}(),
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
    histogram_recipe=PB.make_render(
        Spec,
        histogram_only;
        mode = :ecdf,
        data = :pdf
    )
    histogram_page=only(histogram_recipe.pages).payload
    @test only(histogram_page.layers).kind == Val(:line)
    @test only(histogram_page.layers).group === :histogram_cdf
    @test_throws ArgumentError Cmp._mc_values(nothing)
    @test_throws ArgumentError PB.make_render(
        Spec,
        histogram_only;
        mode = :ecdf,
        data = :samples
    )
    @test_throws ArgumentError PB.make_render(Spec, histogram_only; mode = :qq)
end
