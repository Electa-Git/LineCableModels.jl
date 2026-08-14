@testitem "PlotBuilder: cable and system preview RenderSpec semantics" setup = [defaults] begin
    library = CablesLibrary()
    load!(library; file_name = joinpath(pkgdir(LineCableModels), "test", "cable_test.json"))
    design = first(values(library.data))
    plot_builder = LineCableModels.PlotBuilder

    cable_render = plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotSpec,
        design
    )
    @test length(cable_render.figures) == 1
    cable_page = only(cable_render.figures)
    cable_view = only(cable_page.views)
    @test cable_page.layout === :preview
    @test cable_view.xaxis.label == "y [m]"
    @test cable_view.yaxis.label == "z [m]"
    @test cable_view.attributes.aspect === :data
    @test all(series -> series.kind === :polygon, cable_view.series)
    @test all(series -> haskey(series.attributes, :group), cable_view.series)
    @test all(series -> haskey(series.attributes, :color), cable_view.series)
    @test any(
        series -> hasproperty(series.zdata, :interiors) &&
                  !isempty(series.zdata.interiors),
        cable_view.series
    )
    @test length(cable_page.kwargs.colorbars) == 3
    @test cable_page.kwargs.colorbars[2].ticks == ([0.5], ["1"])
    @test all(
        descriptor -> length(unique(descriptor.ticks[2])) == length(descriptor.ticks[2]),
        cable_page.kwargs.colorbars)
    @test cable_page.kwargs.display_legend
    @test cable_page.kwargs.controls.reset
    @test cable_page.kwargs.controls.export_svg
    @test cable_page.kwargs.export_theme === :default
    @test cable_page.kwargs.open_export
    @test !cable_page.kwargs.controls.xlog
    @test !cable_page.kwargs.controls.ylog
    @test cable_page.kwargs.configuration.x_offset == 0.0
    @test cable_page.kwargs.configuration.y_offset == 0.0

    cable_without_chrome = plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotSpec,
        design;
        display_legend = false,
        display_colorbars = false
    )
    @test !only(cable_without_chrome.figures).kwargs.display_legend
    @test isempty(only(cable_without_chrome.figures).kwargs.colorbars)

    publication_cable = plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotSpec,
        design;
        export_theme = :publication,
        open_export = false
    )
    @test only(publication_cable.figures).kwargs.export_theme === :publication
    @test !only(publication_cable.figures).kwargs.open_export
    @test_throws ArgumentError plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotSpec,
        design;
        export_theme = :unsupported
    )

    position = CablePosition(
        design,
        0.0,
        -0.20,
        Dict(component.id => (index == 1 ? 1 : 0)
        for (index, component) in enumerate(design.components))
    )
    system = LineCableSystem("render-spec-system", 1000.0, position)
    earth = EarthModel([50.0, 100.0], 100.0, 10.0, 1.0)
    system_render = plot_builder.make_render(
        LineCableModels.DataModel.SystemPreviewPlotSpec,
        system;
        earth_model = earth
    )
    system_page = only(system_render.figures)
    system_view = only(system_page.views)
    @test system_view.attributes.aspect === :data
    @test any(series -> series.kind === :hline, system_view.series)
    @test any(series -> series.kind === :polygon, system_view.series)
    earth_reference = only(series
    for series in system_view.series if series.kind === :hline)
    @test earth_reference.ydata == [0.0]
    @test earth_reference.attributes.color === :black
    @test earth_reference.attributes.linewidth == 1.5
    @test length(system_page.kwargs.colorbars) == 3
    @test all(
        descriptor -> length(descriptor.ticks[1]) == 1 &&
                      length(descriptor.ticks[2]) == 1,
        system_page.kwargs.colorbars
    )
    @test getproperty.(system_page.kwargs.colorbars, :ticks) ==
          (([0.5], ["100"]), ([0.5], ["1"]), ([0.5], ["10"]))
    @test system_page.kwargs.configuration.zoom_factor === nothing

    zoomed_render = plot_builder.make_render(
        LineCableModels.DataModel.SystemPreviewPlotSpec,
        system;
        earth_model = earth,
        zoom_factor = 0.5
    )
    default_limits = system_view.attributes.limits
    zoomed_limits = only(only(zoomed_render.figures).views).attributes.limits
    @test zoomed_limits[1][2] - zoomed_limits[1][1] <
          default_limits[1][2] - default_limits[1][1]
    @test zoomed_limits[2][2] - zoomed_limits[2][1] <
          default_limits[2][2] - default_limits[2][1]
    @test_throws ArgumentError plot_builder.make_render(
        LineCableModels.DataModel.SystemPreviewPlotSpec,
        system;
        zoom_factor = 0.0
    )

    scale_render = plot_builder.make_render(
        LineCableModels.DataModel.MaterialScalePlotSpec,
        nothing
    )
    scale_page = only(scale_render.figures)
    @test scale_page.layout === :material_scale
    @test isempty(scale_page.views)
    @test length(scale_page.kwargs.colorbars) == 3
    @test !scale_page.kwargs.controls.reset
    @test scale_page.kwargs.controls.export_svg
    @test scale_page.kwargs.export_theme === :default
    @test scale_page.kwargs.open_export
    @test !scale_page.kwargs.controls.legend
end
