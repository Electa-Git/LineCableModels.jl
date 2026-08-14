@testitem "PlotBuilder: dispatch recipe grammar" setup = [defaults] begin
    const PB = LineCableModels.PlotBuilder
    const UH = LineCableModels.UnitHandler

    struct ProfileResult
        frequency::Vector{Float64}
        response::Matrix{Float64}
    end

    struct ProfilePlotSpec <: PB.AbstractPlotSpec end

    PB.dispatch_on(::Type{ProfilePlotSpec}) = ProfileResult
    PB.input_kwargs(::Type{ProfilePlotSpec}) = (:grouping, :color)
    PB.renderer_kwargs(::Type{ProfilePlotSpec}) = (:size,)
    PB.input_defaults(::Type{ProfilePlotSpec}, ::ProfileResult) = (;
        grouping = :overlay, color = :steelblue)
    PB.renderer_defaults(::Type{ProfilePlotSpec}, ::ProfileResult) = (; size = (800, 400))
    PB.recipe_mode(::Type{ProfilePlotSpec}, recipe::NamedTuple) = Val(:profile)
    PB.grouping_mode(::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple) = Val(recipe.input.grouping)
    PB.group_facets(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple, page_key) = axes(
        recipe.object.response, 2)

    PB.axis_quantity(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x}, recipe::NamedTuple,
        page_key, view_key) = UH.QuantityTag{:freq}()
    PB.axis_quantity(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y}, recipe::NamedTuple,
        page_key, view_key) = UH.QuantityTag{:dimensionless}()
    PB.axis_unit(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x},
        quantity::UH.QuantityTag, recipe::NamedTuple, page_key, view_key) = UH.units(:base, :hertz)
    PB.axis_label(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y},
        quantity::UH.QuantityTag, unit::UH.Units, recipe::NamedTuple,
        page_key, view_key) = "Response"

    PB.series_data(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x}, recipe::NamedTuple,
        page_key, view_key, series_key::Int) = recipe.object.frequency
    PB.series_data(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y}, recipe::NamedTuple,
        page_key, view_key, series_key::Int) = recipe.object.response[:, series_key]
    PB.legend_label(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple,
        page_key, view_key, series_key::Int) = "response $series_key"
    profile_linestyle(::Val) = :solid
    profile_linestyle(::Val{2}) = :dash
    PB.series_attributes(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple,
        page_key, view_key, series_key::Int) = (;
        color = recipe.input.color,
        linewidth = 2,
        linestyle = profile_linestyle(Val(series_key)),
        group = Symbol("response_$series_key")
    )

    profile_title(::Val{:overlay}, page_key, view_key) = "Frequency responses"
    profile_title(::Val{:panels}, page_key, view_key::Int) = "Response $view_key"
    profile_title(::Val{:panels}, page_key, ::Nothing) = "Frequency responses"
    profile_title(::Val{:pages}, page_key::Int, view_key) = "Response $page_key"
    PB.default_title(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple,
        page_key, view_key) = profile_title(Val(recipe.input.grouping), page_key, view_key)

    profile_layout(::Val{:overlay}) = :single
    profile_layout(::Val{:panels}) = :grid
    profile_layout(::Val{:pages}) = :single
    PB.figure_layout(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple,
        page_key) = profile_layout(Val(recipe.input.grouping))
    PB.default_figsize(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple, page_key) = recipe.renderer.size

    result = ProfileResult(
        [50.0, 100.0, 500.0],
        [1.0 1.2; 1.5 1.6; 2.0 2.1]
    )

    overlay = PB.make_render(ProfilePlotSpec, result; color = :navy)
    panels = PB.make_render(ProfilePlotSpec, result; grouping = :panels)
    pages = PB.make_render(ProfilePlotSpec, result; grouping = :pages)

    @test length(only(only(overlay.figures).views).series) == 2
    @test only(only(overlay.figures).views).series[1].attributes.color === :navy
    @test length(only(panels.figures).views) == 2
    @test all(view -> length(view.series) == 1, only(panels.figures).views)
    @test only(panels.figures).layout === :grid
    @test length(pages.figures) == 2
    @test all(page -> length(only(page.views).series) == 1, pages.figures)
    @test_throws ArgumentError PB.make_render(
        ProfilePlotSpec, result; grouping = :unsupported)
    @test_throws ArgumentError PB.make_render(ProfilePlotSpec, [1.0, 2.0])
    @test basename(String(which(PB.make_render, (Type{ProfilePlotSpec}, ProfileResult)).file)) ==
          "grammar.jl"

    maintained_recipes = (
        (LineCableModels.Engine.LineParameterPlotSpec,
            LineCableModels.Engine.LineParameters),
        (LineCableModels.UQ.MCDistributionPlotSpec, LineCableModels.UQ.CableConstantsMC),
        (LineCableModels.DataModel.CablePreviewPlotSpec,
            LineCableModels.DataModel.CableDesign),
        (LineCableModels.DataModel.SystemPreviewPlotSpec,
            LineCableModels.DataModel.LineCableSystem),
        (LineCableModels.DataModel.MaterialScalePlotSpec, Nothing)
    )
    for (specification, object_type) in maintained_recipes
        method = which(PB.make_render, (Type{specification}, object_type))
        @test basename(String(method.file)) == "grammar.jl"
    end
end
