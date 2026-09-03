@testitem "Plotting / core boundary / shell handles and domain geometry" tags=[:unit] setup=[
    TestFixtures
] begin
    using LineCableModels

    @test isdefined(LineCableModels, :PlotBuilder)
    @test parentmodule(LineCableModels.plot) === LineCableModels.PlotBuilder
    @test parentmodule(LineCableModels.preview) === LineCableModels.PlotBuilder
    @test parentmodule(LineCableModels.show_material_scale) ===
          LineCableModels.PlotBuilder
    @test !isdefined(LineCableModels, :set_backend!)
    @test !isdefined(LineCableModels.Engine, :plot)
    @test !isdefined(LineCableModels.DataModel, :preview)
    @test !isdefined(LineCableModels.DataModel, :show_material_scale)
    @test !isdefined(LineCableModels.Engine, :_prepare_line_observations)
    @test !isdefined(LineCableModels.DataModel, :NativePreviewPolygon)
    @test !isdefined(LineCableModels.DataModel, :NativePreviewReference)
    @test !isdefined(LineCableModels.DataModel, :material_color_scheme)

    @test fieldnames(Base.unwrap_unionall(UIPlot)) == (
        :figure, :title, :axes, :controls, :legend, :panel_legends, :colorbars,
        :addon_state, :export_name, :export_theme, :open_export
    )

    figure=Ref(:native_figure)
    axes=(:axis_1, :axis_2)
    controls=Dict{Symbol, Any}(:reset=>:native_button)
    legend=Ref(:native_legend)
    colorbars=(Ref(:native_colorbar),)
    addon_state=(semantic_groups = (:a, :b),)
    plot_handle=UIPlot(
        figure,
        axes;
        controls,
        legend,
        colorbars,
        addon_state,
        export_name = "native",
        open_export = false
    )
    @test plot_handle.figure === figure
    @test plot_handle.title === nothing
    @test plot_handle.axes === axes
    @test plot_handle.controls === controls
    @test plot_handle.legend === legend
    @test isempty(plot_handle.panel_legends)
    @test plot_handle.colorbars === colorbars
    @test plot_handle.addon_state === addon_state
    @test !hasproperty(plot_handle, :context)
    @test !hasproperty(plot_handle, :panels)
    @test !hasproperty(plot_handle, :legend_data)

    design=TestFixtures.mv_cable_design()
    placed=first(design.geometry.regions)
    shapes=LineCableModels.DataModel.preview_shapes(placed)
    @test !isempty(shapes)
    @test all(shape -> shape isa LineCableModels.DataModel.PreviewShape, shapes)
    @test all(shape -> shape.material === placed.source.material, shapes)
    @test all(shape -> shape.tag === placed.source.tag, shapes)

    copper=Material(kind = :conductor, rho = 1.7241e-8)
    matrix=Material(kind = :insulator, rho = Inf, eps_r = 2.3)
    packed=Enclosure(
        :preview_matrix,
        terminal(
            :preview_core,
            stranded(
                copper;
                shape = Disk(0.5),
                layers = 1,
                boundary = Disk(1.4)
            )
        );
        primitive = Disk(1.5),
        fill = matrix
    )
    fill=last(LineCableModels.DataModel.resolve(
        LineCableModels.DataModel.EmptyBoundary(), packed
    ).regions)
    perforated=only(LineCableModels.DataModel.preview_shapes(fill))
    @test length(perforated.geometry.interiors) == 7
    @test eltype(perforated.geometry.interiors) !== Any

    ranges=LineCableModels.DataModel.material_property_ranges(design)
    @test keys(ranges) == (:rho, :mu_r, :eps_r)
    @test all(range -> range isa Tuple{<:Real, <:Real}, values(ranges))

    root=pkgdir(LineCableModels)
    @test !ispath(joinpath(
        root, "src", "engine", "lineparameters", "plotdata.jl"
    ))
    @test !ispath(joinpath(
        root, "src", "engine", "lineparameters", "comparisondata.jl"
    ))
    @test isfile(joinpath(
        root, "ext", "LineCableModelsMakieExt", "recipes", "line_data.jl"
    ))
    @test isfile(joinpath(
        root, "ext", "LineCableModelsMakieExt", "recipes", "preview_data.jl"
    ))
end
