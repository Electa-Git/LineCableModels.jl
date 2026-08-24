@testitem "PlotBuilder / grammar / dispatch and render construction" tags=[:unit] setup=[
    PlotBuilderTestSupport, UsePlotBuilderSupport, TestNumerics, TestFixtures] begin
    const PB=LineCableModels.PlotBuilder
    const UH=LineCableModels.QuantityUnits

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

    struct ProfilePlotDefinition<:PB.AbstractPlotDefinition end

    PB.dispatch_on(::Type{ProfilePlotDefinition}) = ProfileResult
    PB.input_kwargs(::Type{ProfilePlotDefinition}) = (:grouping, :color)
    PB.renderer_kwargs(::Type{ProfilePlotDefinition}) = (:size,)
    PB.input_defaults(::Type{ProfilePlotDefinition}, ::ProfileResult) = (;
        grouping = :overlay, color = :steelblue)
    PB.renderer_defaults(::Type{ProfilePlotDefinition}, ::ProfileResult) = (;
        size = (800, 400))
    PB._recipe_variant(::Type{ProfilePlotDefinition}, recipe::PB.PlotRecipe) = Val(:profile)
    PB._composition(::Type{ProfilePlotDefinition}, ::Val{:profile},
        recipe::PB.PlotRecipe) = Val(recipe.input.grouping)
    PB._series_items(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, recipe::PB.PlotRecipe, page_key) = axes(
        recipe.object.response, 2)

    PB.axis_quantity(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, ::Val{:x}, recipe::PB.PlotRecipe,
        page_key, view_key) = UH.QuantityTag{:frequency}()
    PB.axis_quantity(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, ::Val{:y}, recipe::PB.PlotRecipe,
        page_key, view_key) = UH.QuantityTag{:dimensionless}()
    PB.axis_unit(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, ::Val{:x},
        quantity::UH.QuantityTag, recipe::PB.PlotRecipe, page_key, view_key) = UH.units(:base, :hertz)
    PB.axis_label(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, ::Val{:y},
        quantity::UH.QuantityTag, unit::UH.Units, recipe::PB.PlotRecipe,
        page_key, view_key) = "Response"

    PB.series_data(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, ::Val{:x}, recipe::PB.PlotRecipe,
        page_key, view_key, series_key::Int) = recipe.object.frequency
    PB.series_data(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, ::Val{:y}, recipe::PB.PlotRecipe,
        page_key, view_key, series_key::Int) = recipe.object.response[:, series_key]
    PB.legend_label(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, recipe::PB.PlotRecipe,
        page_key, view_key, series_key::Int) = "response $series_key"
    profile_linestyle(::Val) = :solid
    profile_linestyle(::Val{2}) = :dash
    PB.series_attributes(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, recipe::PB.PlotRecipe,
        page_key, view_key, series_key::Int) = (;
        color = recipe.input.color,
        linewidth = 2,
        linestyle = profile_linestyle(Val(series_key))
    )
    PB.series_group(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, recipe::PB.PlotRecipe,
        page_key, view_key, series_key::Int) = Symbol("response_$series_key")

    profile_title(::Val{:overlay}, page_key, view_key) = "Frequency responses"
    profile_title(::Val{:panels}, page_key, view_key::Int) = "Response $view_key"
    profile_title(::Val{:panels}, page_key, ::Nothing) = "Frequency responses"
    profile_title(::Val{:pages}, page_key::Int, view_key) = "Response $page_key"
    PB.default_title(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, recipe::PB.PlotRecipe,
        page_key, view_key) = profile_title(Val(recipe.input.grouping), page_key, view_key)

    profile_layout(::Val{:overlay}) = :single
    profile_layout(::Val{:panels}) = :grid
    profile_layout(::Val{:pages}) = :single
    PB.layout_spec(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, recipe::PB.PlotRecipe,
        page_key) = profile_layout(Val(recipe.input.grouping))
    PB.default_figsize(
        ::Type{ProfilePlotDefinition}, ::Val{:profile}, recipe::PB.PlotRecipe, page_key) = recipe.renderer.size

    result=ProfileResult(
        [50.0, 100.0, 500.0],
        [1.0 1.2; 1.5 1.6; 2.0 2.1]
    )

    overlay=PB.make_render(ProfilePlotDefinition, result; color = :navy)
    panels=PB.make_render(ProfilePlotDefinition, result; grouping = :panels)
    pages=PB.make_render(ProfilePlotDefinition, result; grouping = :pages)

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
    @test root_grid.columns[1] isa PB.RelativeTrack
    @test side_slots[:legend].halign === :left
    @test side_slots[:colorbars].halign === :left
    side_grid=only(grid for grid in overlay_page.layout.grids if grid.name===:side)
    @test side_grid.rows[1] isa PB.RelativeTrack
    @test side_grid.rows[2] isa PB.ContentTrack
    @test overlay_page.legend.overflow === :ellipsis
    @test PB.LegendSpec(overflow = :show_all).overflow === :show_all
    @test_throws ArgumentError PB.LegendSpec(overflow = :truncate)

    content_legend_grids=map(overlay_page.layout.grids) do grid
        grid.name===:side||return grid
        return PB.GridSpec(
            grid.name;
            parent = grid.parent,
            area = grid.area,
            rows = PB.AbstractTrackSize[PB.ContentTrack(), PB.ContentTrack()],
            columns = grid.columns,
            rowgap = grid.rowgap,
            columngap = grid.columngap,
            padding = grid.padding
        )
    end
    content_legend_layout=PB.LayoutSpec(
        :content_legend,
        content_legend_grids,
        overlay_page.layout.slots
    )
    @test_throws ArgumentError PB.PageSpec(
        overlay_page.title,
        overlay_page.size,
        overlay_page.key,
        content_legend_layout,
        overlay_page.views;
        controls = overlay_page.controls,
        legend = PB.LegendSpec(overflow = :ellipsis),
        colorbars = overlay_page.colorbars,
        status = overlay_page.status,
        export_spec = overlay_page.export_spec
    )
    @test PB.PageSpec(
        overlay_page.title,
        overlay_page.size,
        overlay_page.key,
        content_legend_layout,
        overlay_page.views;
        controls = overlay_page.controls,
        legend = PB.LegendSpec(overflow = :show_all),
        colorbars = overlay_page.colorbars,
        status = overlay_page.status,
        export_spec = overlay_page.export_spec
    ).legend.overflow === :show_all
    compact_objects=Any[
        PB.parse_kwargs(ProfilePlotDefinition, result),
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
          "PlotRecipe(spec=:ProfilePlotDefinition, pages=1, object=:ProfileResult)"
    @test_throws ArgumentError PB.make_render(
        ProfilePlotDefinition, result; grouping = :unsupported)
    @test_throws ArgumentError PB.make_render(ProfilePlotDefinition, [1.0, 2.0])
    @test basename(String(which(PB.make_render, (
        Type{ProfilePlotDefinition}, ProfileResult)).file)) ==
          "render.jl"

    maintained_recipes=(
        (LineCableModels.Engine.LineParameterPlotDefinition,
            LineCableModels.Engine.LineParameters),
        (LineCableModels.Engine.LineParametersBenchmarkPlotDefinition,
            Tuple{
                LineCableModels.Engine.LineParameters,
                LineCableModels.Engine.LineParameters
            }),
        (LineCableModels.UQ.MCDistributionPlotDefinition,
            LineCableModels.MonteCarloResult),
        (LineCableModels.DataModel.CablePreviewPlotDefinition,
            LineCableModels.DataModel.CableDesign),
        (LineCableModels.DataModel.SystemPreviewPlotDefinition,
            LineCableModels.DataModel.LineCableSystem),
        (LineCableModels.DataModel.MaterialScalePlotDefinition, Nothing)
    )
    for (specification, object_type) in maintained_recipes
        method=which(PB.make_render, (Type{specification}, object_type))
        @test basename(String(method.file)) == "render.jl"
    end

    mc_result=TestFixtures.cable_monte_carlo_result()
    for mode in (:hist, :pdf, :ecdf, :qq)
        rendered=PB.make_render(
            LineCableModels.UQ.MCDistributionPlotDefinition,
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
        mc_result.formulation,
        mc_result.values,
        mc_result.stats,
        mc_result.sample_values,
        nothing,
        mc_result.root_seed,
        mc_result.point_seeds,
        mc_result.trial_counts
    )
    @test PB.make_render(
        LineCableModels.UQ.MCDistributionPlotDefinition,
        samples_only;
        mode = :pdf
    ) isa PB.PlotRecipe

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
        ProfilePlotDefinition,
        result;
        grouping = :panels,
        layout = custom_layout
    )
    @test only(custom.figures).layout === custom_layout
    @test all(view -> view.placement.slot === :canvas, only(custom.figures).views)
    @test only(PB.make_render(ProfilePlotDefinition, result; layout = :grid).figures).layout.name ===
          :grid
    @test_throws ArgumentError PB.make_render(ProfilePlotDefinition, result; layout = :unknown)
    @test_throws ArgumentError PB.make_render(ProfilePlotDefinition, result; typo = true)

    struct BadDefaultsPlotDefinition<:PB.AbstractPlotDefinition end
    PB.dispatch_on(::Type{BadDefaultsPlotDefinition}) = ProfileResult
    PB.input_kwargs(::Type{BadDefaultsPlotDefinition}) = (:mode,)
    PB.input_defaults(::Type{BadDefaultsPlotDefinition}, ::ProfileResult) = (;)
    @test_throws ArgumentError PB.make_render(BadDefaultsPlotDefinition, result)

    struct DuplicateOptionPlotDefinition<:PB.AbstractPlotDefinition end
    PB.dispatch_on(::Type{DuplicateOptionPlotDefinition}) = ProfileResult
    PB.input_kwargs(::Type{DuplicateOptionPlotDefinition}) = (:mode,)
    PB.renderer_kwargs(::Type{DuplicateOptionPlotDefinition}) = (:mode,)
    PB.input_defaults(::Type{DuplicateOptionPlotDefinition}, ::ProfileResult) = (;
        mode = :a)
    PB.renderer_defaults(::Type{DuplicateOptionPlotDefinition}, ::ProfileResult) = (;
        mode = :b)
    @test_throws ArgumentError PB.make_render(DuplicateOptionPlotDefinition, result)

    struct CommonOptionPlotDefinition<:PB.AbstractPlotDefinition end
    PB.dispatch_on(::Type{CommonOptionPlotDefinition}) = ProfileResult
    PB.input_kwargs(::Type{CommonOptionPlotDefinition}) = (:layout,)
    PB.input_defaults(::Type{CommonOptionPlotDefinition}, ::ProfileResult) = (;
        layout = :single)
    @test_throws ArgumentError PB.make_render(CommonOptionPlotDefinition, result)

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

    switchable_axis=PB.AxisSpec(
        :y,
        UH.QuantityTag{:dimensionless}(),
        UH.Units(),
        "switchable",
        :linear;
        allowed_scales = (:linear, :log10)
    )
    switchable_view=PB.ViewSpec(
        axis,
        switchable_axis,
        nothing,
        "linear limits",
        [line],
        (; kind = :linear_limits);
        limits = ((1.0, 2.0), (-1.0, 1.0))
    )
    @test PB.PageSpec(
        "linear limits",
        (400, 300),
        (; kind = :linear_limits),
        PB.layout_preset(Val(:single), 1),
        [switchable_view]
    ) isa PB.PageSpec

    active_log_axis=PB.AxisSpec(
        :y,
        UH.QuantityTag{:dimensionless}(),
        UH.Units(),
        "active log",
        :log10;
        allowed_scales = (:linear, :log10)
    )
    active_log_view=PB.ViewSpec(
        axis,
        active_log_axis,
        nothing,
        "invalid logarithmic limits",
        [line],
        (; kind = :invalid_log_limits);
        limits = ((1.0, 2.0), (-1.0, 1.0))
    )
    @test_throws DomainError PB.PageSpec(
        "invalid logarithmic limits",
        (400, 300),
        (; kind = :invalid_log_limits),
        PB.layout_preset(Val(:single), 1),
        [active_log_view]
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
    const UH=LineCableModels.QuantityUnits

    struct DefaultHookSpec<:PB.AbstractPlotDefinition end
    recipe=PB.parse_kwargs(DefaultHookSpec, :payload)
    variant=PB._recipe_variant(DefaultHookSpec, recipe)

    @test PB.dispatch_on(DefaultHookSpec) === Any
    @test PB.resolve_input(DefaultHookSpec, recipe) === recipe
    @test PB.observe(DefaultHookSpec, recipe) === recipe
    @test variant === Val(:default)
    @test !isdefined(PB, :recipe_mode)
    @test !isdefined(PB, :grouping_mode)
    @test !isdefined(PB, :page_facets)
    @test !isdefined(PB, :group_facets)
    @test !isdefined(PB, :page_keys)
    @test !isdefined(PB, :view_keys)
    @test !isdefined(PB, :series_keys)
    @test !isdefined(PB, :RenderSpec)
    @test !isdefined(PB, :LineFamilyKey)
    @test !isdefined(PB, :MCSeriesKey)
    @test all(name -> name ∉ names(PB),
        (
            :AxisSpec, :SeriesSpec, :ViewSpec, :PageSpec, :LayoutSpec,
            :ControlSpec, :LegendSpec, :ColorbarSpec, :StatusSpec, :ExportSpec
        ))
    @test PB.geom_axes(DefaultHookSpec, variant, recipe, nothing, nothing) == (:x, :y)

    quantity=PB.axis_quantity(DefaultHookSpec, Val(:x), recipe)
    @test quantity isa UH.QuantityTag{:unknown}
    @test PB.axis_quantity(
        DefaultHookSpec,
        variant,
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
