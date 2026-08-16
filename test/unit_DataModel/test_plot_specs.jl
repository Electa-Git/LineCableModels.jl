@testitem "PlotBuilder: cable and system preview RenderSpec semantics" setup = [defaults] begin
    @test :show_material_scale ∉ names(LineCableModels)
    @test :show_material_scale ∉ names(LineCableModels.DataModel)

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
    @test cable_page.layout.name === :preview
    @test cable_view.xaxis.label == "y [m]"
    @test cable_view.yaxis.label == "z [m]"
    @test cable_view.aspect === :data
    @test all(series -> series.kind === :polygon, cable_view.series)
    @test all(series -> series.group isa Symbol, cable_view.series)
    @test all(series -> haskey(series.attributes, :color), cable_view.series)
    legend_labels = String[
        series.label for series in cable_view.series if series.label !== nothing
    ]
    @test !isempty(legend_labels)
    @test all(label -> !occursin("ρ=", label), legend_labels)
    @test all(label -> occursin(": ", label), legend_labels)
    @test length(unique(legend_labels)) == length(legend_labels)
    @test any(
        series -> hasproperty(series.zdata, :interiors) &&
                  !isempty(series.zdata.interiors),
        cable_view.series
    )
    @test length(cable_page.colorbars) == 3
    @test cable_page.colorbars[2].ticks == ([0.5], ["1"])
    @test all(
        descriptor -> length(unique(descriptor.ticks[2])) == length(descriptor.ticks[2]),
        cable_page.colorbars)
    @test cable_page.legend.enabled
    @test cable_page.controls.reset
    @test cable_page.controls.export_svg
    @test cable_page.export_spec.theme === :default
    @test cable_page.export_spec.open_file
    @test cable_page.key == (; kind = :cable, id = design.cable_id)

    cable_without_chrome = plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotSpec,
        design;
        display_legend = false,
        display_colorbars = false
    )
    @test !only(cable_without_chrome.figures).legend.enabled
    @test isempty(only(cable_without_chrome.figures).colorbars)

    publication_cable = plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotSpec,
        design;
        export_theme = :publication,
        open_export = false
    )
    @test only(publication_cable.figures).export_spec.theme === :publication
    @test !only(publication_cable.figures).export_spec.open_file
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
    @test system_view.aspect === :data
    @test any(series -> series.kind === :hline, system_view.series)
    @test any(series -> series.kind === :polygon, system_view.series)
    earth_reference = only(series
    for series in system_view.series if series.kind === :hline)
    @test earth_reference.ydata == [0.0]
    @test earth_reference.attributes.color === :black
    @test earth_reference.attributes.linewidth == 1.5
    @test length(system_page.colorbars) == 3
    @test all(
        descriptor -> length(descriptor.ticks[1]) == 1 &&
                      length(descriptor.ticks[2]) == 1,
        system_page.colorbars
    )
    @test getproperty.(system_page.colorbars, :ticks) ==
          [([0.5], ["100"]), ([0.5], ["1"]), ([0.5], ["10"])]

    zoomed_render = plot_builder.make_render(
        LineCableModels.DataModel.SystemPreviewPlotSpec,
        system;
        earth_model = earth,
        zoom_factor = 0.5
    )
    default_limits = system_view.limits
    zoomed_limits = only(only(zoomed_render.figures).views).limits
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
    @test scale_page.layout.name === :material_scale
    @test isempty(scale_page.views)
    @test length(scale_page.colorbars) == 3
    @test !scale_page.controls.reset
    @test scale_page.controls.export_svg
    @test scale_page.export_spec.theme === :default
    @test scale_page.export_spec.open_file
    @test !scale_page.legend.enabled
end
