@testitem "PlotBuilder / final recipe grammar and Template Method" tags=[:unit] setup=[
    PlotBuilderTestSupport, UsePlotBuilderSupport, TestFixtures
] begin
    const PB=LineCableModels.PlotBuilder
    const U=LineCableModels.Units

    if get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false")!="true"
        loaded=Set(nameof(module_value) for module_value in values(Base.loaded_modules))
        @test :Makie ∉ loaded
        @test :CairoMakie ∉ loaded
        @test :GLMakie ∉ loaded
        @test :WGLMakie ∉ loaded
    end

    @test length(methods(PB.make_render)) == 1
    @test fieldnames(PB.PlotPage) == (:title, :size, :key, :payload)
    @test fieldnames(PB.PlotRecipe) == (:definition, :pages)

    struct ProfileResult
        frequency::Vector{Float64}
        response::Matrix{Float64}
    end
    struct ProfilePlotDefinition<:PB.AbstractPlotDefinition end

    profile_frequency(result::ProfileResult) = result.frequency
    profile_response(result::ProfileResult) = result.response

    LineCableModels.basis(::ProfileResult) = :total
    function LineCableModels.observe(
            result::ProfileResult,
            ::typeof(profile_frequency),
            indices...
    )
        return isempty(indices) ? result.frequency : getindex(result.frequency, indices...)
    end
    function LineCableModels.observe(
            result::ProfileResult,
            ::typeof(profile_response),
            indices...
    )
        return isempty(indices) ? result.response : getindex(result.response, indices...)
    end
    LineCableModels.observables(::Type{<:ProfileResult}) = (
        profile_frequency, profile_response)
    U.quantity(::typeof(profile_frequency)) = U.Quantity{:frequency}()
    U.quantity(::typeof(profile_response)) = U.Quantity{:test_response}()
    U.native_unit(::U.Quantity{:test_response}) = U.UnitExpr()
    U.display_unit(::U.Quantity{:test_response}) = U.UnitExpr()
    U.label(::U.Quantity{:test_response}) = "Response"
    U.symbol(::U.Quantity{:test_response}) = "u"

    const profile_stage_calls=Symbol[]
    PB.dispatch_on(::Type{ProfilePlotDefinition}) = ProfileResult
    function PB.entitle(::Type{ProfilePlotDefinition}, source::ProfileResult)
        push!(profile_stage_calls, :entitle)
        return source
    end
    PB.input_kwargs(::Type{ProfilePlotDefinition}) = (:grouping, :color)
    PB.renderer_kwargs(::Type{ProfilePlotDefinition}) = (:size,)
    function PB.input_defaults(::Type{ProfilePlotDefinition}, ::ProfileResult)
        push!(profile_stage_calls, :parse)
        return (; grouping = :overlay, color = :steelblue)
    end
    PB.renderer_defaults(::Type{ProfilePlotDefinition}, ::ProfileResult) = (;
        size = (800, 400))
    function PB.resolve(
            ::Type{ProfilePlotDefinition},
            ::ProfileResult,
            request::NamedTuple
    )
        push!(profile_stage_calls, :resolve)
        request.input.grouping in (:overlay, :pages)||throw(ArgumentError(
            "grouping must be :overlay or :pages",
        ))
        return request
    end
    function PB.fetch(
            ::Type{ProfilePlotDefinition},
            source::ProfileResult,
            request::NamedTuple
    )
        push!(profile_stage_calls, :fetch)
        published=observables(
            source,
            (frequency = profile_frequency, response = profile_response)
        )
        legend=PB.LegendDefinition()
        export_definition=PB.ExportDefinition(
            theme = request.renderer.export_theme,
            name = "profile",
            open_file = request.renderer.open_export
        )
        count=request.input.grouping===:overlay ? 1 : size(source.response, 2)
        return Tuple(
            PB.PlotPage(
                "Profile $index",
                request.renderer.size,
                (; page = index),
                (;
                    published,
                    color = request.input.color,
                    legend,
                    colorbars = (),
                    export_definition
                )
            )
        for index in 1:count)
    end
    function PB.finish(
            ::Type{ProfilePlotDefinition},
            pages::P
    ) where {P <: Tuple}
        push!(profile_stage_calls, :finish)
        return PB.PlotRecipe(ProfilePlotDefinition, pages)
    end

    source=ProfileResult(
        [50.0, 100.0, 500.0],
        [1.0 1.2; 1.5 1.6; 2.0 2.1]
    )
    empty!(profile_stage_calls)
    recipe=PB.make_render(ProfilePlotDefinition, source; color = :navy)
    @test profile_stage_calls == [:entitle, :parse, :resolve, :fetch, :finish]
    @test recipe.definition === ProfilePlotDefinition
    @test recipe.pages isa Tuple
    @test length(recipe.pages) == 1
    @test only(recipe.pages).payload.color === :navy
    @test only(recipe.pages).payload.published.response.values == source.response
    @test !hasproperty(recipe, :source)
    @test !hasproperty(recipe, :object)
    @test !hasproperty(recipe, :input)
    @test !hasproperty(recipe, :renderer)
    @test sprint(show, MIME"text/plain"(), recipe) ==
          "PlotRecipe(definition=:ProfilePlotDefinition, pages=1)"

    paged=PB.make_render(ProfilePlotDefinition, source; grouping = :pages)
    @test length(paged.pages) == 2
    @test getfield.(getfield.(paged.pages, :key), :page) == (1, 2)
    @test_throws ArgumentError PB.make_render(
        ProfilePlotDefinition,
        source;
        grouping = :unsupported
    )
    @test_throws ArgumentError PB.make_render(ProfilePlotDefinition, source; layout = :grid)
    @test_throws ArgumentError PB.make_render(ProfilePlotDefinition, source; typo = true)

    empty!(profile_stage_calls)
    @test_throws ArgumentError PB.make_render(ProfilePlotDefinition, [1.0, 2.0])
    @test isempty(profile_stage_calls)

    struct MissingDispatchDefinition<:PB.AbstractPlotDefinition end
    @test_throws MethodError PB.make_render(MissingDispatchDefinition, source)

    struct MissingResolveDefinition<:PB.AbstractPlotDefinition end
    PB.dispatch_on(::Type{MissingResolveDefinition}) = ProfileResult
    @test_throws MethodError PB.make_render(MissingResolveDefinition, source)

    struct MissingFetchDefinition<:PB.AbstractPlotDefinition end
    PB.dispatch_on(::Type{MissingFetchDefinition}) = ProfileResult
    PB.resolve(::Type{MissingFetchDefinition}, ::ProfileResult, request::NamedTuple) = request
    @test_throws MethodError PB.make_render(MissingFetchDefinition, source)

    struct BadDefaultsDefinition<:PB.AbstractPlotDefinition end
    PB.dispatch_on(::Type{BadDefaultsDefinition}) = ProfileResult
    PB.input_kwargs(::Type{BadDefaultsDefinition}) = (:mode,)
    PB.input_defaults(::Type{BadDefaultsDefinition}, ::ProfileResult) = (;)
    @test_throws ArgumentError PB.make_render(BadDefaultsDefinition, source)

    struct DuplicateOptionDefinition<:PB.AbstractPlotDefinition end
    PB.dispatch_on(::Type{DuplicateOptionDefinition}) = ProfileResult
    PB.input_kwargs(::Type{DuplicateOptionDefinition}) = (:mode,)
    PB.renderer_kwargs(::Type{DuplicateOptionDefinition}) = (:mode,)
    PB.input_defaults(::Type{DuplicateOptionDefinition}, ::ProfileResult) = (; mode = :a)
    PB.renderer_defaults(::Type{DuplicateOptionDefinition}, ::ProfileResult) = (; mode = :b)
    @test_throws ArgumentError PB.make_render(DuplicateOptionDefinition, source)

    page=only(recipe.pages)
    @test_throws ArgumentError PB.PlotPage("invalid", (0, 400), (; page = 0), (;))
    @test_throws ArgumentError PB.PlotRecipe(ProfilePlotDefinition, (page, page))
    @test_throws ArgumentError PB.PlotRecipe(ProfilePlotDefinition, (page, :not_a_page))

    legend=PB.LegendDefinition(overflow = :show_all)
    @test legend.overflow === :show_all
    @test_throws ArgumentError PB.LegendDefinition(overflow = :truncate)
    @test_throws DimensionMismatch PB.ColorbarDefinition(
        "invalid ticks",
        :viridis,
        (0.0, 1.0),
        ([0.0, 1.0], ["zero"])
    )
    @test_throws ArgumentError PB.ColorbarDefinition(
        "out of range",
        :viridis,
        (0.0, 1.0),
        ([2.0], ["two"])
    )
    @test PB.ExportDefinition(theme = :publication).theme === :publication
    @test_throws ArgumentError PB.ExportDefinition(theme = :unknown)

    compact=Any[
        legend,
        PB.ColorbarDefinition("scale", :viridis, (0.0, 1.0), ([0.0], ["0"])),
        PB.ExportDefinition(),
        page,
        recipe
    ]
    for value in compact
        representation=sprint(show, MIME"text/plain"(), value)
        @test count(==('\n'), representation) == 0
        @test ncodeunits(representation) < 180
        @test !occursin("[50.0, 100.0, 500.0]", representation)
    end

    maintained=(
        (LineCableModels.Engine.LineParameterPlotDefinition,
            LineCableModels.Engine.LineParameters),
        (LineCableModels.Engine.LineParametersBenchmarkPlotDefinition,
            Tuple{LineCableModels.Engine.LineParameters,
                LineCableModels.Engine.LineParameters}),
        (LineCableModels.UQ.MCDistributionPlotDefinition,
            LineCableModels.MonteCarloResult),
        (LineCableModels.DataModel.CablePreviewPlotDefinition,
            LineCableModels.DataModel.CableDesign),
        (LineCableModels.DataModel.SystemPreviewPlotDefinition,
            LineCableModels.DataModel.LineCableSystem),
        (LineCableModels.DataModel.MaterialScalePlotDefinition, Nothing)
    )
    for (definition, source_type) in maintained
        method=which(PB.make_render, (Type{definition}, source_type))
        @test basename(String(method.file)) == "render.jl"
    end

    line_recipe=@inferred PB.make_render(
        LineCableModels.Engine.LineParameterPlotDefinition,
        TestFixtures.two_conductor_results()
    )
    mc_recipe=@inferred PB.make_render(
        LineCableModels.UQ.MCDistributionPlotDefinition,
        TestFixtures.cable_monte_carlo_result()
    )
    preview_recipe=@inferred PB.make_render(
        LineCableModels.DataModel.CablePreviewPlotDefinition,
        TestFixtures.mv_cable_design()
    )
    @test line_recipe isa PB.PlotRecipe
    @test mc_recipe isa PB.PlotRecipe
    @test preview_recipe isa PB.PlotRecipe

    compiler_types=(
        :AbstractTrackSize, :FixedTrack, :RelativeTrack, :ContentTrack,
        :GridArea, :GridSpec, :SlotSpec, :LayoutSpec, :PlacementSpec,
        :ControlSpec, :LegendSpec, :ColorbarSpec, :StatusSpec, :ExportSpec,
        :AxisSpec, :SeriesSpec, :ViewSpec, :PageSpec
    )
    compiler_hooks=(
        :layout_preset, :default_title, :default_figsize, :layout_spec,
        :page_identity, :control_spec, :legend_spec, :colorbar_specs,
        :status_spec, :export_spec, :axis_payload, :series_values,
        :series_group, :series_attributes, :plot_kind, :geom_axes,
        :make_axes, :make_series, :make_views, :make_pages,
        :_recipe_variant, :_composition, :_page_keys, :_view_keys,
        :_series_items
    )
    @test all(name -> !isdefined(PB, name), (compiler_types..., compiler_hooks...))
    @test all(
        name -> !isfile(joinpath(pkgdir(LineCableModels), "src", "plotbuilder", name)),
        ("axes.jl", "series.jl", "views.jl", "pages.jl", "composition.jl",
            "validation.jl"))

    shell_source=read(
        joinpath(
            pkgdir(LineCableModels),
            "ext",
            "LineCableModelsMakieExt",
            "shell.jl"
        ), String)
    build_page_start=first(findfirst("function build_page(", shell_source))
    build_page_stop=first(findnext("\nend\n\nfunction build(", shell_source,
        build_page_start))
    build_page_source=shell_source[build_page_start:build_page_stop]
    ui_stages=(
        "build_context(",
        "build_shell(",
        "draw!(",
        "format_axes!(",
        "place_legend!(",
        "place_colorbars!(",
        "build_widgets!(",
        "assemble(",
        "display!("
    )
    positions=map(stage->first(findfirst(stage, build_page_source)), ui_stages)
    @test issorted(positions)

    @test all(name -> name in names(PB), (:plotwindow, :axis!, :register!))
    @test all(name -> name ∉ names(LineCableModels), (:plotwindow, :axis!, :register!))
end
