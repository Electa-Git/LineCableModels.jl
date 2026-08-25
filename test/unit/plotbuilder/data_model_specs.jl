@testitem "PlotBuilder / data models / preview PlotRecipe semantics" tags=[:unit] setup=[
    PlotBuilderTestSupport, UsePlotBuilderSupport, TestNumerics] begin
    @test :show_material_scale ∉ names(LineCableModels)
    @test :show_material_scale ∉ names(LineCableModels.DataModel)
    @test all(name -> name in names(LineCableModels.PlotBuilder),
        (:plotwindow, :axis!, :register!))
    @test all(name -> name ∉ names(LineCableModels),
        (:plotwindow, :axis!, :register!))

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
        LineCableModels.DataModel.CablePreviewPlotDefinition,
        design
    )
    @test length(cable_render.figures) == 1
    cable_page=only(cable_render.figures)
    cable_payload=cable_render.input.payload
    @test cable_page.layout.name === :preview
    @test isempty(cable_page.views)
    @test all(polygon -> polygon.group isa Symbol, cable_payload.polygons)
    @test all(polygon -> polygon.color !== nothing, cable_payload.polygons)
    legend_labels=String[polygon.label
                         for polygon in cable_payload.polygons if polygon.label!==nothing]
    @test !isempty(legend_labels)
    @test all(label -> !occursin("ρ=", label), legend_labels)
    @test all(label -> occursin(": ", label), legend_labels)
    @test length(unique(legend_labels)) == length(legend_labels)
    @test any(
        polygon -> hasproperty(polygon.geometry, :interiors) &&
                   !isempty(polygon.geometry.interiors),
        cable_payload.polygons
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
        LineCableModels.DataModel.CablePreviewPlotDefinition,
        design;
        display_legend = false,
        display_colorbars = false
    )
    @test !only(cable_without_chrome.figures).legend.enabled
    @test isempty(only(cable_without_chrome.figures).colorbars)

    publication_cable=plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotDefinition,
        design;
        export_theme = :publication,
        open_export = false
    )
    @test only(publication_cable.figures).export_spec.theme === :publication
    @test !only(publication_cable.figures).export_spec.open_file
    @test_throws ArgumentError plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotDefinition,
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
        LineCableModels.DataModel.SystemPreviewPlotDefinition,
        system;
        earth_model = earth
    )
    system_page=only(system_render.figures)
    system_payload=system_render.input.payload
    @test isempty(system_page.views)
    @test !isempty(system_payload.polygons)
    earth_reference=only(system_payload.references)
    @test earth_reference.values == [0.0]
    @test earth_reference.color === :black
    @test earth_reference.width == 1.5
    @test length(system_page.colorbars) == 3
    @test all(
        descriptor -> length(descriptor.ticks[1]) == 1 &&
                      length(descriptor.ticks[2]) == 1,
        system_page.colorbars
    )
    @test getproperty.(system_page.colorbars, :ticks) ==
          [([0.5], ["100"]), ([0.5], ["1"]), ([0.5], ["10"])]

    zoomed_render=plot_builder.make_render(
        LineCableModels.DataModel.SystemPreviewPlotDefinition,
        system;
        earth_model = earth,
        zoom_factor = 0.5
    )
    default_limits=system_payload.limits
    zoomed_limits=zoomed_render.input.payload.limits
    @test zoomed_limits[1][2] - zoomed_limits[1][1] <
          default_limits[1][2] - default_limits[1][1]
    @test zoomed_limits[2][2] - zoomed_limits[2][1] <
          default_limits[2][2] - default_limits[2][1]
    @test_throws ArgumentError plot_builder.make_render(
        LineCableModels.DataModel.SystemPreviewPlotDefinition,
        system;
        zoom_factor = 0.0
    )

    scale_render=plot_builder.make_render(
        LineCableModels.DataModel.MaterialScalePlotDefinition,
        nothing
    )
    scale_page=only(scale_render.figures)
    @test scale_page.layout.name === :material_scale
    @test isempty(scale_page.views)
    @test length(scale_page.colorbars) == 3
    @test scale_page.controls.reset
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

    shapes=DM.PreviewPolygon[]
    @test DM._layer_shapes!(shapes, strands, "sector", :sector, 0.0, 0.0) === shapes
    @test length(shapes) == strands.num_wires
    @test count(item -> item.label == "sector", shapes) == 1
    @test all(item -> length(item.geometry) == 64, shapes)

    group=ConductorGroup(strands)
    grouped_shapes=DM.PreviewPolygon[]
    DM._layer_shapes!(grouped_shapes, group, "group", :group, 0.0, 0.0)
    @test length(grouped_shapes) == strands.num_wires
    @test count(item -> item.label == "group", grouped_shapes) == 1

    semicon=Semicon(strands.r_ex, strands.r_ex+0.001, semiconductor)
    insulator=Insulator(semicon.r_ex, semicon.r_ex+0.002, dielectric)
    @test DM._preview_layer_name(strands) == "rectangular strands"
    @test DM._preview_layer_name(semicon) == "semiconductor"
    @test DM._preview_layer_name(insulator) == "insulation"
    @test DM._preview_layer_name(group) == "conductor group"

    struct UnsupportedPreviewLayer
        material_props
    end
    unsupported=UnsupportedPreviewLayer(conductor)
    @test DM._preview_layer_name(unsupported) == "unsupportedpreviewlayer"
    @test_logs (:warn, r"unsupported cable-preview layer") begin
        @test isempty(DM._layer_shapes!(
            DM.PreviewPolygon[],
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
