@testitem "PlotBuilder / grammar / dispatch and render construction" tags=[:unit] setup=[
    PlotBuilderTestSupport, UsePlotBuilderSupport, TestNumerics] begin
    const PB=LineCableModels.PlotBuilder
    const UH=LineCableModels.UnitHandler

    @test length(methods(PB.make_render)) == 1
    if get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false")!="true"
        loaded_module_names=Set(nameof(module_value)
        for module_value in values(Base.loaded_modules))
        @test :Makie ∉ loaded_module_names
        @test :CairoMakie ∉ loaded_module_names
        @test :GLMakie ∉ loaded_module_names
        @test :WGLMakie ∉ loaded_module_names
    end

    struct ProfileResult
        frequency::Vector{Float64}
        response::Matrix{Float64}
    end

    struct ProfilePlotSpec<:PB.AbstractPlotSpec end

    PB.dispatch_on(::Type{ProfilePlotSpec}) = ProfileResult
    PB.input_kwargs(::Type{ProfilePlotSpec}) = (:grouping, :color)
    PB.renderer_kwargs(::Type{ProfilePlotSpec}) = (:size,)
    PB.input_defaults(::Type{ProfilePlotSpec}, ::ProfileResult) = (;
        grouping = :overlay, color = :steelblue)
    PB.renderer_defaults(::Type{ProfilePlotSpec}, ::ProfileResult) = (; size = (800, 400))
    PB.recipe_mode(::Type{ProfilePlotSpec}, recipe::PB.PlotRecipe) = Val(:profile)
    PB.grouping_mode(::Type{ProfilePlotSpec}, ::Val{:profile},
        recipe::PB.PlotRecipe) = Val(recipe.input.grouping)
    PB.group_facets(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe, page_key) = axes(
        recipe.object.response, 2)

    PB.axis_quantity(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x}, recipe::PB.PlotRecipe,
        page_key, view_key) = UH.QuantityTag{:freq}()
    PB.axis_quantity(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y}, recipe::PB.PlotRecipe,
        page_key, view_key) = UH.QuantityTag{:dimensionless}()
    PB.axis_unit(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x},
        quantity::UH.QuantityTag, recipe::PB.PlotRecipe, page_key, view_key) = UH.units(:base, :hertz)
    PB.axis_label(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y},
        quantity::UH.QuantityTag, unit::UH.Units, recipe::PB.PlotRecipe,
        page_key, view_key) = "Response"

    PB.series_data(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x}, recipe::PB.PlotRecipe,
        page_key, view_key, series_key::Int) = recipe.object.frequency
    PB.series_data(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y}, recipe::PB.PlotRecipe,
        page_key, view_key, series_key::Int) = recipe.object.response[:, series_key]
    PB.legend_label(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe,
        page_key, view_key, series_key::Int) = "response $series_key"
    profile_linestyle(::Val) = :solid
    profile_linestyle(::Val{2}) = :dash
    PB.series_attributes(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe,
        page_key, view_key, series_key::Int) = (;
        color = recipe.input.color,
        linewidth = 2,
        linestyle = profile_linestyle(Val(series_key))
    )
    PB.series_group(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe,
        page_key, view_key, series_key::Int) = Symbol("response_$series_key")

    profile_title(::Val{:overlay}, page_key, view_key) = "Frequency responses"
    profile_title(::Val{:panels}, page_key, view_key::Int) = "Response $view_key"
    profile_title(::Val{:panels}, page_key, ::Nothing) = "Frequency responses"
    profile_title(::Val{:pages}, page_key::Int, view_key) = "Response $page_key"
    PB.default_title(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe,
        page_key, view_key) = profile_title(Val(recipe.input.grouping), page_key, view_key)

    profile_layout(::Val{:overlay}) = :single
    profile_layout(::Val{:panels}) = :grid
    profile_layout(::Val{:pages}) = :single
    PB.layout_spec(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe,
        page_key) = profile_layout(Val(recipe.input.grouping))
    PB.default_figsize(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe, page_key) = recipe.renderer.size

    result=ProfileResult(
        [50.0, 100.0, 500.0],
        [1.0 1.2; 1.5 1.6; 2.0 2.1]
    )

    overlay=PB.make_render(ProfilePlotSpec, result; color = :navy)
    panels=PB.make_render(ProfilePlotSpec, result; grouping = :panels)
    pages=PB.make_render(ProfilePlotSpec, result; grouping = :pages)

    @test length(only(only(overlay.figures).views).series) == 2
    @test only(only(overlay.figures).views).series[1].attributes.color === :navy
    @test length(only(panels.figures).views) == 2
    @test all(view -> length(view.series) == 1, only(panels.figures).views)
    @test only(panels.figures).layout.name === :grid
    @test length(pages.figures) == 2
    @test all(page -> length(only(page.views).series) == 1, pages.figures)
    overlay_page=only(overlay.figures)
    overlay_view=only(overlay_page.views)
    root_grid=only(grid for grid in overlay_page.layout.grids if grid.parent===nothing)
    side_slots=Dict(slot.name=>slot for slot in overlay_page.layout.slots)
    @test root_grid.columngap == 12
    @test side_slots[:legend].halign === :left
    @test side_slots[:colorbars].halign === :left
    compact_objects=Any[
        PB.parse_kwargs(ProfilePlotSpec, result),
        PB.FixedTrack(36),
        PB.RelativeTrack(),
        PB.ContentTrack(),
        PB.GridArea(1, 1:2),
        first(overlay_page.layout.grids),
        first(overlay_page.layout.slots),
        overlay_page.layout,
        overlay_view.placement,
        overlay_page.controls,
        overlay_page.legend,
        PB.ColorbarSpec(
            "Response scale",
            :viridis,
            (0.0, 2.5),
            ([0.0, 2.5], ["0", "2.5"])
        ),
        overlay_page.status,
        overlay_page.export_spec,
        overlay_view.xaxis,
        first(overlay_view.series),
        overlay_view,
        overlay_page,
        overlay
    ]
    for object in compact_objects
        representation=sprint(show, MIME"text/plain"(), object)
        @test count(==('\n'), representation) == 0
        @test ncodeunits(representation) < 180
        @test !occursin("[50.0, 100.0, 500.0]", representation)
    end
    @test sprint(show, MIME"text/plain"(), overlay) ==
          "RenderSpec(spec=:ProfilePlotSpec, pages=1)"
    @test_throws ArgumentError PB.make_render(
        ProfilePlotSpec, result; grouping = :unsupported)
    @test_throws ArgumentError PB.make_render(ProfilePlotSpec, [1.0, 2.0])
    @test basename(String(which(PB.make_render, (Type{ProfilePlotSpec}, ProfileResult)).file)) ==
          "grammar.jl"

    maintained_recipes=(
        (LineCableModels.Engine.LineParameterPlotSpec,
            LineCableModels.Engine.LineParameters),
        (LineCableModels.Computation.MCDistributionPlotSpec,
            LineCableModels.MonteCarloResult),
        (LineCableModels.DataModel.CablePreviewPlotSpec,
            LineCableModels.DataModel.CableDesign),
        (LineCableModels.DataModel.SystemPreviewPlotSpec,
            LineCableModels.DataModel.LineCableSystem),
        (LineCableModels.DataModel.MaterialScalePlotSpec, Nothing)
    )
    for (specification, object_type) in maintained_recipes
        method=which(PB.make_render, (Type{specification}, object_type))
        @test basename(String(method.file)) == "grammar.jl"
    end

    summary=SampleSummary([1.0, 2.0, 3.0, 4.0])
    histogram=HistogramPDF([1.0, 3.0, 5.0], [0.25, 0.25])
    mc_result=MonteCarloResult(
        CableConstants(2.5, 2.5, 2.5),
        CableConstants(summary, summary, summary),
        CableConstants(
            [1.0, 2.0, 3.0, 4.0],
            [1.0, 2.0, 3.0, 4.0],
            [1.0, 2.0, 3.0, 4.0]
        ),
        CableConstants(histogram, histogram, histogram),
        nothing,
        4,
        0.95,
        0.02,
        :normal,
        UInt64(1),
        (hash = "plot-fixture",)
    )
    for mode in (:hist, :pdf, :ecdf, :qq)
        rendered=PB.make_render(
            LineCableModels.Computation.MCDistributionPlotSpec,
            mc_result;
            mode,
            data = :both
        )
        page=only(rendered.figures)
        @test page.key.mode === mode
        @test page.key.quantity === :R
        @test only(page.views).key.selection === nothing
    end
    samples_only=MonteCarloResult(
        mc_result.representation,
        mc_result.statistics,
        mc_result.samples,
        nothing,
        nothing,
        mc_result.trials,
        mc_result.confidence,
        mc_result.cdf_tol,
        mc_result.distribution,
        mc_result.seed,
        mc_result.manifest
    )
    @test PB.make_render(
        LineCableModels.Computation.MCDistributionPlotSpec,
        samples_only;
        mode = :pdf
    ) isa PB.RenderSpec

    custom_layout=PB.LayoutSpec(
        :profile_dashboard,
        [
            PB.GridSpec(
                :root;
                rows = PB.AbstractTrackSize[
                    PB.FixedTrack(36), PB.RelativeTrack(), PB.FixedTrack(20)],
                columns = PB.AbstractTrackSize[PB.ContentTrack(), PB.ContentTrack()],
                rowgap = 6,
                columngap = 6,
                padding = (20, 20, 28, 28)
            ),
            PB.GridSpec(
                :plots;
                parent = :root,
                area = PB.GridArea(2, 1),
                rows = PB.AbstractTrackSize[PB.RelativeTrack()],
                columns = PB.AbstractTrackSize[PB.RelativeTrack()]
            )
        ],
        [
            PB.SlotSpec(:toolbar, :root, PB.GridArea(1, 1:2)),
            PB.SlotSpec(:canvas, :plots, PB.GridArea(1, 1)),
            PB.SlotSpec(:legend, :root, PB.GridArea(2, 2)),
            PB.SlotSpec(:status, :root, PB.GridArea(3, 1:2))
        ]
    )
    custom=PB.make_render(
        ProfilePlotSpec,
        result;
        grouping = :panels,
        layout = custom_layout
    )
    @test only(custom.figures).layout === custom_layout
    @test all(view -> view.placement.slot === :canvas, only(custom.figures).views)
    @test only(PB.make_render(ProfilePlotSpec, result; layout = :grid).figures).layout.name ===
          :grid
    @test_throws ArgumentError PB.make_render(ProfilePlotSpec, result; layout = :unknown)
    @test_throws ArgumentError PB.make_render(ProfilePlotSpec, result; typo = true)

    struct BadDefaultsPlotSpec<:PB.AbstractPlotSpec end
    PB.dispatch_on(::Type{BadDefaultsPlotSpec}) = ProfileResult
    PB.input_kwargs(::Type{BadDefaultsPlotSpec}) = (:mode,)
    PB.input_defaults(::Type{BadDefaultsPlotSpec}, ::ProfileResult) = (;)
    @test_throws ArgumentError PB.make_render(BadDefaultsPlotSpec, result)

    struct DuplicateOptionPlotSpec<:PB.AbstractPlotSpec end
    PB.dispatch_on(::Type{DuplicateOptionPlotSpec}) = ProfileResult
    PB.input_kwargs(::Type{DuplicateOptionPlotSpec}) = (:mode,)
    PB.renderer_kwargs(::Type{DuplicateOptionPlotSpec}) = (:mode,)
    PB.input_defaults(::Type{DuplicateOptionPlotSpec}, ::ProfileResult) = (; mode = :a)
    PB.renderer_defaults(::Type{DuplicateOptionPlotSpec}, ::ProfileResult) = (; mode = :b)
    @test_throws ArgumentError PB.make_render(DuplicateOptionPlotSpec, result)

    struct CommonOptionPlotSpec<:PB.AbstractPlotSpec end
    PB.dispatch_on(::Type{CommonOptionPlotSpec}) = ProfileResult
    PB.input_kwargs(::Type{CommonOptionPlotSpec}) = (:layout,)
    PB.input_defaults(::Type{CommonOptionPlotSpec}, ::ProfileResult) = (; layout = :single)
    @test_throws ArgumentError PB.make_render(CommonOptionPlotSpec, result)

    @test_throws ArgumentError PB.AxisSpec(
        :x,
        UH.QuantityTag{:dimensionless}(),
        UH.Units(),
        "x",
        :log10;
        allowed_scales = (:linear,)
    )
    @test_throws ArgumentError PB.AxisSpec(
        :x,
        UH.QuantityTag{:dimensionless}(),
        UH.Units(),
        "x";
        attributes = (; xlabel = "hidden semantic override")
    )
    @test_throws ArgumentError PB.SeriesSpec(
        :line,
        [1.0],
        [2.0],
        nothing,
        "reserved";
        attributes = (; group = :invalid)
    )
    @test_throws ArgumentError PB.SeriesSpec(
        :line,
        [1.0],
        [2.0],
        nothing,
        "reserved";
        attributes = (; label = "hidden semantic override")
    )
    @test_throws ArgumentError PB.ViewSpec(
        nothing,
        nothing,
        nothing,
        "reserved",
        PB.SeriesSpec[],
        (;);
        attributes = (; limits = ((0.0, 1.0), (0.0, 1.0)))
    )
    @test_throws DimensionMismatch PB.ColorbarSpec(
        "invalid ticks",
        :viridis,
        (0.0, 1.0),
        ([0.0, 1.0], ["zero"])
    )
    @test_throws ArgumentError PB.ColorbarSpec(
        "out-of-range ticks",
        :viridis,
        (0.0, 1.0),
        ([2.0], ["two"])
    )
    @test_throws ArgumentError PB.GridArea(2:1, 1:1)
    @test_throws ArgumentError PB.FixedTrack(-1)
    @test_throws ArgumentError PB.RelativeTrack(0)

    bad_parent=PB.GridSpec(
        :child;
        parent = :missing,
        area = PB.GridArea(1, 1)
    )
    @test_throws ArgumentError PB.LayoutSpec(
        :missing_parent,
        [PB.GridSpec(:root), bad_parent],
        PB.SlotSpec[]
    )
    @test_throws ArgumentError PB.LayoutSpec(
        :overlap,
        [PB.GridSpec(
            :root;
            rows = PB.AbstractTrackSize[PB.RelativeTrack()],
            columns = PB.AbstractTrackSize[PB.RelativeTrack()]
        )],
        [
            PB.SlotSpec(:first, :root, PB.GridArea(1, 1)),
            PB.SlotSpec(:second, :root, PB.GridArea(1, 1))
        ]
    )

    cycle_first=PB.GridSpec(
        :first;
        parent = :second,
        area = PB.GridArea(1, 1)
    )
    cycle_second=PB.GridSpec(
        :second;
        parent = :first,
        area = PB.GridArea(1, 1)
    )
    @test_throws ArgumentError PB.LayoutSpec(
        :cycle,
        [cycle_first, cycle_second],
        PB.SlotSpec[]
    )

    axis=PB.AxisSpec(
        :x,
        UH.QuantityTag{:dimensionless}(),
        UH.Units(),
        "x"
    )
    yaxis=PB.AxisSpec(
        :y,
        UH.QuantityTag{:dimensionless}(),
        UH.Units(),
        "y"
    )
    line=PB.SeriesSpec(:line, [1.0], [2.0], nothing, "line")
    automatic=PB.ViewSpec(axis, yaxis, nothing, "automatic", [line], (;))
    explicit=PB.ViewSpec(
        axis,
        yaxis,
        nothing,
        "explicit",
        [line],
        (; explicit = true);
        placement = PB.PlacementSpec(:canvas, PB.GridArea(1, 1))
    )
    @test_throws ArgumentError PB.PageSpec(
        "mixed placements",
        (400, 300),
        (; kind = :mixed),
        PB.layout_preset(Val(:single), 2),
        [automatic, explicit]
    )

    duplicate_key=PB.ViewSpec(
        axis,
        yaxis,
        nothing,
        "duplicate",
        [line],
        (; explicit = true)
    )
    @test_throws ArgumentError PB.PageSpec(
        "duplicate view identity",
        (400, 300),
        (; kind = :duplicate),
        PB.layout_preset(Val(:single), 2),
        [explicit, duplicate_key]
    )

    negative_axis=PB.AxisSpec(
        :y,
        UH.QuantityTag{:dimensionless}(),
        UH.Units(),
        "negative",
        :linear;
        allowed_scales = (:linear, :log10)
    )
    negative_view=PB.ViewSpec(
        axis,
        negative_axis,
        nothing,
        "negative logarithmic data",
        [PB.SeriesSpec(:line, [1.0], [-1.0], nothing, "negative")],
        (; kind = :negative)
    )
    @test_throws DomainError PB.PageSpec(
        "negative logarithmic data",
        (400, 300),
        (; kind = :negative),
        PB.layout_preset(Val(:single), 1),
        [negative_view]
    )

    reversed_limits_view=PB.ViewSpec(
        axis,
        yaxis,
        nothing,
        "reversed limits",
        [line],
        (; kind = :reversed_limits);
        limits = ((2.0, 1.0), (0.0, 1.0))
    )
    @test_throws ArgumentError PB.PageSpec(
        "reversed limits",
        (400, 300),
        (; kind = :reversed_limits),
        PB.layout_preset(Val(:single), 1),
        [reversed_limits_view]
    )

    invalid_aspect_view=PB.ViewSpec(
        axis,
        yaxis,
        nothing,
        "invalid aspect",
        [line],
        (; kind = :invalid_aspect);
        aspect = :square
    )
    @test_throws ArgumentError PB.PageSpec(
        "invalid aspect",
        (400, 300),
        (; kind = :invalid_aspect),
        PB.layout_preset(Val(:single), 1),
        [invalid_aspect_view]
    )

    missing_slot_layout=PB.LayoutSpec(
        :missing_slot,
        [PB.GridSpec(:root)],
        [PB.SlotSpec(:canvas, :root, PB.GridArea(1, 1))]
    )
    @test_throws ArgumentError PB.PageSpec(
        "missing typed component slots",
        (400, 300),
        (; kind = :missing_slot),
        missing_slot_layout,
        [automatic]
    )
end

@testitem "PlotBuilder / grammar / default hook contracts" tags=[:unit] setup=[
    PlotBuilderTestSupport,
    UsePlotBuilderSupport
] begin
    const PB=LineCableModels.PlotBuilder
    const UH=LineCableModels.UnitHandler

    struct DefaultHookSpec<:PB.AbstractPlotSpec end
    recipe=PB.parse_kwargs(DefaultHookSpec, :payload)
    mode=PB.recipe_mode(DefaultHookSpec, recipe)

    @test PB.dispatch_on(DefaultHookSpec) === Any
    @test PB.resolve_input(DefaultHookSpec, recipe) === recipe
    @test mode === Val(:default)
    @test PB.grouping_mode(DefaultHookSpec, mode, recipe) === Val(:overlay)
    @test PB.page_facets(DefaultHookSpec, mode, recipe) == (nothing,)
    @test PB.group_facets(DefaultHookSpec, mode, recipe, nothing) == (nothing,)
    @test PB.page_keys(DefaultHookSpec, mode, Val(:overlay), recipe) == (nothing,)
    @test PB.page_keys(DefaultHookSpec, mode, Val(:panels), recipe) == (nothing,)
    @test PB.page_keys(DefaultHookSpec, mode, Val(:pages), recipe) == (nothing,)
    @test PB.page_keys(DefaultHookSpec, mode, Val(:faceted_pages), recipe) == (nothing,)
    @test PB.page_keys(DefaultHookSpec, mode, Val(:empty), recipe) == (nothing,)
    @test_throws ArgumentError PB.page_keys(
        DefaultHookSpec,
        mode,
        Val(:unsupported),
        recipe
    )

    @test PB.view_keys(DefaultHookSpec, mode, Val(:faceted_pages), recipe, nothing) ==
          (nothing,)
    @test PB.series_keys(
        DefaultHookSpec,
        mode,
        Val(:faceted_pages),
        recipe,
        nothing,
        nothing
    ) == (nothing,)
    @test PB.geom_axes(DefaultHookSpec, mode, recipe, nothing, nothing) == (:x, :y)

    quantity=PB.axis_quantity(DefaultHookSpec, Val(:x), recipe)
    @test quantity isa UH.QuantityTag{:unknown}
    @test PB.axis_quantity(
        DefaultHookSpec,
        mode,
        Val(:x),
        recipe,
        nothing,
        nothing
    ) === quantity
    unit=PB.axis_unit(DefaultHookSpec, Val(:x), quantity, recipe)
    @test PB.axis_label(DefaultHookSpec, Val(:x), quantity, unit, recipe) == "unknown"
    @test PB.axis_scale(DefaultHookSpec, Val(:x), recipe) === :linear
    @test PB.axis_scales(DefaultHookSpec, Val(:x), recipe, PB.SeriesSpec[]) == (:linear,)
    @test PB.axis_exponent(DefaultHookSpec, Val(:x), recipe, PB.SeriesSpec[]) == 0

    @test PB.plot_kind(DefaultHookSpec, recipe, nothing) === :line
    @test PB.series_data(DefaultHookSpec, Val(:x), recipe, nothing) === nothing
    @test PB.legend_label(DefaultHookSpec, recipe, nothing) === nothing
    @test PB.series_group(DefaultHookSpec, recipe, nothing) === nothing
    @test PB.series_visible(DefaultHookSpec, recipe, nothing)
    @test PB.series_attributes(DefaultHookSpec, recipe, nothing) == NamedTuple()
    @test PB.default_title(DefaultHookSpec, recipe) == ""
    @test PB.default_figsize(DefaultHookSpec, recipe) == (800, 400)
    @test PB.layout_spec(DefaultHookSpec, recipe) === :single
    @test PB.page_identity(DefaultHookSpec, recipe, nothing) == NamedTuple()
    @test PB.page_identity(DefaultHookSpec, recipe, (quantity = :R,)) == (quantity = :R,)
    @test PB.page_identity(DefaultHookSpec, recipe, :phase) == (facet = :phase,)

    @test_throws ArgumentError PB.make_render(DefaultHookSpec, :payload)
end
