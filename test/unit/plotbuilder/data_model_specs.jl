@testitem "PlotBuilder / data models / preview RenderSpec semantics" tags=[:unit] setup=[
    PlotBuilderTestSupport, UsePlotBuilderSupport, TestNumerics] begin
    @test :show_material_scale ∉ names(LineCableModels)
    @test :show_material_scale ∉ names(LineCableModels.DataModel)

    library=CablesLibrary()
    load!(library;
        file_name = joinpath(
            pkgdir(LineCableModels),
            "test",
            "fixtures",
            "data",
            "mv_cable_design.json"
        ))
    design=first(values(library.data))
    plot_builder=LineCableModels.PlotBuilder

    cable_render=plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotSpec,
        design
    )
    @test length(cable_render.figures) == 1
    cable_page=only(cable_render.figures)
    cable_view=only(cable_page.views)
    @test cable_page.layout.name === :preview
    @test cable_view.xaxis.label == "y [m]"
    @test cable_view.yaxis.label == "z [m]"
    @test cable_view.aspect === :data
    @test all(series -> series.kind === :polygon, cable_view.series)
    @test all(series -> series.group isa Symbol, cable_view.series)
    @test all(series -> haskey(series.attributes, :color), cable_view.series)
    legend_labels=String[series.label
                         for series in cable_view.series if series.label!==nothing]
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

    cable_without_chrome=plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotSpec,
        design;
        display_legend = false,
        display_colorbars = false
    )
    @test !only(cable_without_chrome.figures).legend.enabled
    @test isempty(only(cable_without_chrome.figures).colorbars)

    publication_cable=plot_builder.make_render(
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

    position=CablePosition(
        design,
        0.0,
        -0.20,
        Dict(component.id=>(index==1 ? 1 : 0)
        for (index, component) in enumerate(design.components))
    )
    system=LineCableSystem("render-spec-system", 1000.0, position)
    earth=EarthModel(100.0, 10.0, 1.0)
    system_render=plot_builder.make_render(
        LineCableModels.DataModel.SystemPreviewPlotSpec,
        system;
        earth_model = earth
    )
    system_page=only(system_render.figures)
    system_view=only(system_page.views)
    @test system_view.aspect === :data
    @test any(series -> series.kind === :hline, system_view.series)
    @test any(series -> series.kind === :polygon, system_view.series)
    earth_reference=only(series
    for series in system_view.series if series.kind===:hline)
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

    zoomed_render=plot_builder.make_render(
        LineCableModels.DataModel.SystemPreviewPlotSpec,
        system;
        earth_model = earth,
        zoom_factor = 0.5
    )
    default_limits=system_view.limits
    zoomed_limits=only(only(zoomed_render.figures).views).limits
    @test zoomed_limits[1][2] - zoomed_limits[1][1] <
          default_limits[1][2] - default_limits[1][1]
    @test zoomed_limits[2][2] - zoomed_limits[2][1] <
          default_limits[2][2] - default_limits[2][1]
    @test_throws ArgumentError plot_builder.make_render(
        LineCableModels.DataModel.SystemPreviewPlotSpec,
        system;
        zoom_factor = 0.0
    )

    scale_render=plot_builder.make_render(
        LineCableModels.DataModel.MaterialScalePlotSpec,
        nothing
    )
    scale_page=only(scale_render.figures)
    @test scale_page.layout.name === :material_scale
    @test isempty(scale_page.views)
    @test length(scale_page.colorbars) == 3
    @test !scale_page.controls.reset
    @test scale_page.controls.export_svg
    @test scale_page.export_spec.theme === :default
    @test scale_page.export_spec.open_file
    @test !scale_page.legend.enabled
end

@testitem "PlotBuilder / cable geometry / sector strands and fallback identities" tags=[:unit] setup=[
    PlotBuilderTestSupport,
    UsePlotBuilderSupport,
    TestNumerics
] begin
    using GeometryBasics
    import LineCableModels.DataModel as DM

    conductor=LineCableModels.Materials.Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    dielectric=LineCableModels.Materials.Material(1.0e14, 2.3, 1.0, 20.0, 0.0)
    semiconductor=LineCableModels.Materials.Material(1000.0, 1000.0, 1.0, 20.0, 0.0)
    outer=sqrt(0.01^2+12*0.002*0.001/π)
    strands=RectStrands(0.01, outer, 0.001, 0.002, 12, 12.0, conductor)

    wedge=DM._radial_wedge(0.01, strands.r_ex, strands.width, 0.0, 1.0, -2.0)
    @test length(wedge) == 64
    @test all(point -> isfinite(point[1]) && isfinite(point[2]), wedge)
    collapsed=DM._radial_wedge(0.0, 0.01, 0.002, 0.0, 0.0, 0.0)
    @test length(unique(collapsed)) == 2

    series=LineCableModels.PlotBuilder.SeriesSpec[]
    @test DM._layer_series!(series, strands, "sector", :sector, 0.0, 0.0) === series
    @test length(series) == strands.num_wires
    @test count(item -> item.label == "sector", series) == 1
    @test all(item -> item.kind === :polygon, series)
    @test all(item -> length(item.zdata) == 64, series)

    group=ConductorGroup(strands)
    grouped_series=LineCableModels.PlotBuilder.SeriesSpec[]
    DM._layer_series!(grouped_series, group, "group", :group, 0.0, 0.0)
    @test length(grouped_series) == strands.num_wires
    @test count(item -> item.label == "group", grouped_series) == 1

    semicon=Semicon(strands.r_ex, strands.r_ex+0.001, semiconductor)
    insulator=Insulator(semicon.r_ex, semicon.r_ex+0.002, dielectric)
    @test DM._preview_layer_name(strands) == "sector strands"
    @test DM._preview_layer_name(semicon) == "semiconductor"
    @test DM._preview_layer_name(insulator) == "insulation"
    @test DM._preview_layer_name(group) == "conductor group"

    struct UnsupportedPreviewLayer
        material_props
    end
    unsupported=UnsupportedPreviewLayer(conductor)
    @test DM._preview_layer_name(unsupported) == "unsupportedpreviewlayer"
    @test_logs (:warn, r"unsupported cable-preview layer") begin
        @test isempty(DM._layer_series!(
            LineCableModels.PlotBuilder.SeriesSpec[],
            unsupported,
            "unsupported",
            :unsupported,
            0.0,
            0.0
        ))
    end

    @test DM._base_material_color(Inf) == DM._INSULATOR_COLORS[end]
    @test DM._cable_title(Val(false),
        CableDesign("id", CableComponent(
            "core",
            group,
            InsulatorGroup(semicon)
        ))) == "Cable design preview"
    @test DM._cable_title(Val(true),
        CableDesign("id", CableComponent(
            "core",
            group,
            InsulatorGroup(semicon)
        ))) == "Cable design preview: id"
    @test isempty(DM._cable_colorbars(Val(false), nothing))
    @test DM._system_title(Val(false),
        LineCableSystem(
            "system",
            100.0,
            CablePosition(
                CableDesign("id", CableComponent("core", group, InsulatorGroup(semicon))),
                0.0,
                -1.0,
                Dict("core" => 1)
            )
        )) == "Cable system cross-section"
    @test isempty(DM._system_colorbars(Val(false), nothing))
end
