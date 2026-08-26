@testitem "PlotBuilder / data models / detached preview pages" tags=[:unit] setup=[
    PlotBuilderTestSupport, UsePlotBuilderSupport, TestNumerics] begin
    using GeometryBasics: Point2f

    const DM=LineCableModels.DataModel
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
    @test Base.ispublic(DM, :preview_shapes)
    @test Base.ispublic(DM, :preview_materials)

    preview_context=(;
        label = "test layer",
        group = :test_layer,
        xcenter = 0.0,
        ycenter = 0.0,
        include_label = true
    )
    layers=Any[
        layer
        for component in design.components
        for group in (component.conductor_group, component.insulator_group)
        for layer in group.layers
    ]
    for layer in layers
        @test applicable(DM.preview_shapes, layer, preview_context)
        @test applicable(DM.preview_materials, layer)
        @test !isempty(DM.preview_shapes(layer, preview_context))
        @test !isempty(DM.preview_materials(layer))
    end
    for component in design.components
        @test !isempty(DM.preview_shapes(component.conductor_group, preview_context))
        @test !isempty(DM.preview_materials(component.conductor_group))
        @test !isempty(DM.preview_shapes(component.insulator_group, preview_context))
        @test !isempty(DM.preview_materials(component.insulator_group))
    end

    cable_render=plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotDefinition,
        design
    )
    @test length(cable_render.pages) == 1
    cable_page=only(cable_render.pages)
    cable_payload=cable_page.payload
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
    @test cable_page.export_definition.theme === :default
    @test cable_page.export_definition.open_file
    @test all(
        name -> !hasproperty(cable_payload, name),
        (:title, :key, :legend, :colorbars, :widgets, :export_definition)
    )
    @test cable_page.key == (; kind = :cable, id = design.cable_id)
    @test_throws ArgumentError DM.PreviewPolygon(
        Point2f[(0.0, 0.0), (Inf, 1.0)], nothing, :invalid, :red, :black, 1.0)
    @test_throws ArgumentError DM.PreviewPolygon(
        Point2f[(0.0, 0.0), (1.0, 1.0)], nothing, :invalid, :red, :black, -1.0)
    @test_throws ArgumentError DM.PreviewReferenceLine(
        [Inf], :invalid, :black, 1.0)
    @test_throws ArgumentError DM.PreviewPayload(
        cable_payload.polygons,
        (),
        ((1.0, 0.0), (0.0, 1.0)),
        nothing
    )
    @test_throws ArgumentError DM.PreviewPayload(
        cable_payload.polygons,
        (),
        nothing,
        (hidden_groups = (:absent,), current_limits = nothing)
    )
    @test_throws ArgumentError DM.PreviewPayload(1, (), nothing, nothing)
    @test_throws ArgumentError DM.PreviewPayload((), 1, nothing, nothing)

    cable_without_chrome=plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotDefinition,
        design;
        display_legend = false,
        display_colorbars = false
    )
    @test !only(cable_without_chrome.pages).legend.enabled
    @test isempty(only(cable_without_chrome.pages).colorbars)

    cable_designs=[
        design,
        equivalent(design; new_id = "test_cable_equivalent"),
        CableDesign(
            "test_cable_copy",
            design.components;
            nominal_data = design.nominal_data
        )
    ]
    collection_render=plot_builder.make_render(
        LineCableModels.DataModel.CableCollectionPreviewPlotDefinition,
        cable_designs
    )
    collection_page=only(collection_render.pages)
    @test collection_page.payload.layout == (2, 2)
    @test getproperty.(collection_page.payload.panels, :position) ==
          [(1, 1), (1, 2), (2, 1)]
    @test getproperty.(collection_page.payload.panels, :title) ==
          getproperty.(cable_designs, :cable_id)
    @test all(
        panel -> all(polygon -> polygon.label === nothing, panel.payload.polygons),
        collection_page.payload.panels
    )
    @test !collection_page.legend.enabled
    @test length(collection_page.colorbars) == 3
    @test collection_page.key == (;
        kind = :cable_collection,
        ids = Tuple(getproperty.(cable_designs, :cable_id)),
        layout = (2, 2)
    )
    five_designs=fill(design, 5)
    five_design_page=only(plot_builder.make_render(
        LineCableModels.DataModel.CableCollectionPreviewPlotDefinition,
        five_designs
    ).pages)
    @test five_design_page.payload.layout == (2, 3)

    explicit_collection=plot_builder.make_render(
        LineCableModels.DataModel.CableCollectionPreviewPlotDefinition,
        cable_designs;
        layout = (1, 3),
        display_colorbars = false,
        size = (1000, 400)
    )
    explicit_page=only(explicit_collection.pages)
    @test explicit_page.payload.layout == (1, 3)
    @test isempty(explicit_page.colorbars)
    @test explicit_page.size == (1000, 400)
    @test_throws ArgumentError plot_builder.make_render(
        LineCableModels.DataModel.CableCollectionPreviewPlotDefinition,
        CableDesign[]
    )
    @test_throws ArgumentError plot_builder.make_render(
        LineCableModels.DataModel.CableCollectionPreviewPlotDefinition,
        cable_designs;
        layout = (0, 3)
    )
    @test_throws DimensionMismatch plot_builder.make_render(
        LineCableModels.DataModel.CableCollectionPreviewPlotDefinition,
        cable_designs;
        layout = (1, 2)
    )

    publication_cable=plot_builder.make_render(
        LineCableModels.DataModel.CablePreviewPlotDefinition,
        design;
        export_theme = :publication,
        open_export = false
    )
    @test only(publication_cable.pages).export_definition.theme === :publication
    @test !only(publication_cable.pages).export_definition.open_file
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
    system_page=only(system_render.pages)
    system_payload=system_page.payload
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
          (([0.5], ["100"]), ([0.5], ["1"]), ([0.5], ["10"]))

    zoomed_render=plot_builder.make_render(
        LineCableModels.DataModel.SystemPreviewPlotDefinition,
        system;
        earth_model = earth,
        zoom_factor = 0.5
    )
    default_limits=system_payload.limits
    zoomed_limits=only(zoomed_render.pages).payload.limits
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
    scale_page=only(scale_render.pages)
    @test scale_page.payload === nothing
    @test length(scale_page.colorbars) == 3
    @test scale_page.export_definition.theme === :default
    @test scale_page.export_definition.open_file
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

    sector_context=(;
        label = "sector",
        group = :sector,
        xcenter = 0.0,
        ycenter = 0.0,
        include_label = true
    )
    shapes=DM.preview_shapes(strands, sector_context)
    @test length(shapes) == strands.num_wires
    @test count(item -> item.label == "sector", shapes) == 1
    @test all(item -> length(item.geometry) == 64, shapes)

    group=ConductorGroup(strands)
    grouped_shapes=DM.preview_shapes(
        group,
        merge(sector_context, (; label = "group", group = :group))
    )
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
    @test_throws MethodError DM.preview_shapes(unsupported, sector_context)
    @test_throws MethodError DM.preview_materials(unsupported)

    @test DM._base_material_color(Inf) == DM._INSULATOR_COLORS[end]

    broad_conductor=LineCableModels.Materials.Material(1.0e-6, 4.0, 10.0, 20.0, 0.0)
    broad_dielectric=LineCableModels.Materials.Material(1.0e14, 20.0, 2.0, 20.0, 0.0)
    broad_design=CableDesign(
        "broad-scale",
        ConductorGroup(Tubular(0.0, 0.01, broad_conductor)),
        InsulatorGroup(Insulator(0.01, 0.02, broad_dielectric))
    )
    narrow_design=CableDesign(
        "narrow-scale",
        ConductorGroup(Tubular(0.0, 0.01, conductor)),
        InsulatorGroup(Insulator(0.01, 0.02, dielectric))
    )
    @test DM._property_ranges([narrow_design, broad_design]) ==
          ((1.7241e-8, 1.0e14), (1.0, 10.0), (1.0, 20.0))
    global_scale=only(LineCableModels.PlotBuilder.make_render(
        DM.CableCollectionPreviewPlotDefinition,
        [narrow_design, broad_design]
    ).pages).colorbars
    @test getproperty.(global_scale, :ticks) ==
          (([0.0, 1.0], ["1.72e-08", "1e+14"]),
        ([0.0, 1.0], ["1", "10"]),
        ([0.0, 1.0], ["1", "20"]))
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
