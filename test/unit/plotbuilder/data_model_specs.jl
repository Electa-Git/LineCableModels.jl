@testitem "PlotBuilder / data models / v1 detached preview pages" tags=[:unit] setup=[
    PlotBuilderTestSupport, UsePlotBuilderSupport, TestNumerics, TestFixtures] begin
    using GeometryBasics: Point2f

    const DM=LineCableModels.DataModel
    const PB=LineCableModels.PlotBuilder
    design=TestFixtures.mv_cable_design()

    @test Base.ispublic(DM, :preview_shapes)
    @test Base.ispublic(DM, :preview_materials)
    @test :show_material_scale ∉ names(LineCableModels)
    @test :show_material_scale ∉ names(DM)

    context=(;
        label = "resolved region",
        group = :resolved_region,
        include_label = true
    )
    for region in design.geometry.regions
        @test applicable(DM.preview_shapes, region, context)
        @test applicable(DM.preview_materials, region)
        @test length(DM.preview_shapes(region, context)) == 1
        @test DM.preview_materials(region) == (region.source.material,)
    end

    rendered=PB.make_render(DM.CablePreviewPlotDefinition, design)
    @test length(rendered.pages) == 1
    page=only(rendered.pages)
    payload=page.payload
    @test length(payload.polygons) == length(design.geometry.regions)
    @test all(polygon -> polygon.group isa Symbol, payload.polygons)
    @test all(polygon -> polygon.color !== nothing, payload.polygons)
    labels=String[polygon.label for polygon in payload.polygons
                  if polygon.label!==nothing]
    @test !isempty(labels)
    @test length(unique(labels)) == length(labels)
    @test any(
        polygon -> hasproperty(polygon.geometry, :interiors) &&
                   !isempty(polygon.geometry.interiors),
        payload.polygons
    )
    @test length(page.colorbars) == 3
    @test page.legend.enabled
    @test page.key == (; kind = :cable, id = design.cable_id)
    @test page.export_definition.theme === :default

    @test_throws ArgumentError DM.PreviewPolygon(
        Point2f[(0.0, 0.0), (Inf, 1.0)], nothing, :invalid, :red, :black, 1.0
    )
    @test_throws ArgumentError DM.PreviewPolygon(
        Point2f[(0.0, 0.0), (1.0, 1.0)], nothing, :invalid, :red, :black, -1.0
    )
    @test_throws ArgumentError DM.PreviewReferenceLine([Inf], :invalid, :black, 1.0)
    @test_throws ArgumentError DM.PreviewPayload(
        payload.polygons, (), ((1.0, 0.0), (0.0, 1.0)), nothing
    )
    @test_throws ArgumentError DM.PreviewPayload(
        payload.polygons,
        (),
        nothing,
        (hidden_groups = (:absent,), current_limits = nothing)
    )

    without_chrome=PB.make_render(
        DM.CablePreviewPlotDefinition,
        design;
        display_legend = false,
        display_colorbars = false
    )
    @test !only(without_chrome.pages).legend.enabled
    @test isempty(only(without_chrome.pages).colorbars)

    designs=[
        design,
        build(
            CableDesign,
            "$(design.cable_id)-2",
            design.root
        ),
        build(
            CableDesign,
            "$(design.cable_id)-3",
            design.root
        )
    ]
    collection=PB.make_render(DM.CableCollectionPreviewPlotDefinition, designs)
    collection_page=only(collection.pages)
    @test collection_page.payload.layout == (2, 2)
    @test getproperty.(collection_page.payload.panels, :position) ==
          [(1, 1), (1, 2), (2, 1)]
    @test getproperty.(collection_page.payload.panels, :title) ==
          getproperty.(designs, :cable_id)
    @test all(
        panel -> all(polygon -> polygon.label === nothing, panel.payload.polygons),
        collection_page.payload.panels
    )
    @test !collection_page.legend.enabled
    @test length(collection_page.colorbars) == 3
    @test only(PB.make_render(
        DM.CableCollectionPreviewPlotDefinition,
        fill(design, 5)
    ).pages).payload.layout == (2, 3)

    explicit=PB.make_render(
        DM.CableCollectionPreviewPlotDefinition,
        designs;
        layout = (1, 3),
        display_colorbars = false,
        size = (1000, 400)
    )
    explicit_page=only(explicit.pages)
    @test explicit_page.payload.layout == (1, 3)
    @test isempty(explicit_page.colorbars)
    @test explicit_page.size == (1000, 400)
    @test_throws ArgumentError PB.make_render(
        DM.CableCollectionPreviewPlotDefinition, CableDesign[]
    )
    @test_throws DimensionMismatch PB.make_render(
        DM.CableCollectionPreviewPlotDefinition, designs; layout = (1, 2)
    )

    system=build(
        LineCableSystem,
        design,
        Pose2(0.0, -0.20, 0.0);
        connections = Dict(
            terminal=>(index==1 ? 1 : 0)
        for (index, terminal) in enumerate(design.terminal_order)
        ),
        system_id = "render-spec-system",
        line_length = 1000.0
    )
    earth=EarthModel(100.0, 10.0, 1.0)
    system_page=only(PB.make_render(
        DM.SystemPreviewPlotDefinition,
        system;
        earth_model = earth
    ).pages)
    @test length(system_page.payload.polygons) >= length(system.geometry)
    @test only(system_page.payload.references).values == [0.0]
    @test length(system_page.colorbars) == 3
    @test system_page.key == (; kind = :system, id = system.system_id)

    zoomed=only(PB.make_render(
        DM.SystemPreviewPlotDefinition,
        system;
        earth_model = earth,
        zoom_factor = 0.5
    ).pages).payload.limits
    ordinary=system_page.payload.limits
    @test zoomed[1][2] - zoomed[1][1] < ordinary[1][2] - ordinary[1][1]
    @test zoomed[2][2] - zoomed[2][1] < ordinary[2][2] - ordinary[2][1]
    @test_throws ArgumentError PB.make_render(
        DM.SystemPreviewPlotDefinition, system; zoom_factor = 0.0
    )

    scale=only(PB.make_render(DM.MaterialScalePlotDefinition, nothing).pages)
    @test scale.payload === nothing
    @test length(scale.colorbars) == 3
    @test !scale.legend.enabled
end

@testitem "PlotBuilder / cable geometry / v1 shape dispatch and global scales" tags=[:unit] setup=[
    PlotBuilderTestSupport, UsePlotBuilderSupport, TestNumerics] begin
    import LineCableModels.DataModel as DM

    conductor=Material(kind = :conductor, rho = 1.7241e-8)
    dielectric=Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)
    context=(; label = "shape", group = :shape, include_label = true)
    primitives=(
        DiskDefinition(0.01),
        AnnulusDefinition(0.005, 0.01),
        RectangleDefinition(0.02, 0.01),
        EllipseDefinition(0.02, 0.01),
        SectorDefinition(0.005, 0.01, 0.0, pi/2),
        PolygonDefinition(((-0.01, -0.01), (0.01, -0.01), (0.0, 0.01)))
    )
    for (index, primitive) in enumerate(primitives)
        source=Region(Symbol(:shape_, index), primitive, conductor)
        resolved=only(DM.resolve(DM.EmptyBoundary(), source).regions)
        shape=only(DM.preview_shapes(resolved, context))
        @test shape.label == "shape"
        @test shape.group === :shape
        @test !isempty(shape.geometry.exterior)
        primitive isa Rectangle&&@test length(shape.geometry.exterior) == 4
    end

    struct UnsupportedPreviewLayer end
    unsupported=UnsupportedPreviewLayer()
    @test_throws MethodError DM.preview_shapes(unsupported, context)
    @test_throws MethodError DM.preview_materials(unsupported)

    broad=Material(kind = :conductor, rho = 1.0e-6, mu_r = 10.0, eps_r = 4.0)
    broad_design=build(
        CableDesign,
        "broad-scale",
        Group(:core, Region(:broad, DiskDefinition(0.01), broad))
    )
    narrow_design=build(
        CableDesign,
        "narrow-scale",
        Group(:core, Region(:narrow, DiskDefinition(0.01), conductor))
    )
    @test DM._property_ranges([narrow_design, broad_design]) ==
          ((1.7241e-8, 1.0e-6), (1.0, 10.0), (1.0, 4.0))
    scales=only(LineCableModels.PlotBuilder.make_render(
        DM.CableCollectionPreviewPlotDefinition,
        [narrow_design, broad_design]
    ).pages).colorbars
    @test getproperty.(scales, :ticks) ==
          (([0.0, 1.0], ["1.72e-08", "1e-06"]),
        ([0.0, 1.0], ["1", "10"]),
        ([0.0, 1.0], ["1", "4"]))

    filled=build(
        CableDesign,
        "filled",
        Stack(
            Group(:core, Region(:metal, DiskDefinition(0.01), conductor)),
            Region(:insulation, ShellDefinition(0.005), dielectric)
        )
    )
    @test length(only(LineCableModels.PlotBuilder.make_render(
        DM.CablePreviewPlotDefinition,
        filled
    ).pages).payload.polygons) == 2
end
