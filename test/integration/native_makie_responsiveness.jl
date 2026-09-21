@testitem "Makie addons / responsive native legend and visible-series limits" tags=[:visual] setup=[
    UseNativePlotSupport, TestFixtures
] begin
    get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false")=="true"||
    error("set LINECABLEMODELS_TEST_PLOTTING=true to run the visual contract")
    using CairoMakie

    frequency=10.0 .^ range(1, 4; length = 20)
    parameters=TestFixtures.two_conductor_results(; frequencies = frequency)
    compact=Makie.plot(
        parameters, parameters, parameters, parameters,
        @observe R[1, 1, :];
        series_labels = ("one", "two", "three", "four"),
        backend = :cairo,
        display_plot = false,
        fig_size = (420, 260),
        legend_position = :right,
        legend_overflow = :ellipsis
    )
    Makie.colorbuffer(compact.figure)
    labels=[entry.label[] for entry in last(first(compact.legend.entrygroups[]))]
    @test !isempty(labels)
    @test length(labels) <= 4
    @test last(labels) == "(...)" || length(labels) == 4

    complete=Makie.plot(
        parameters, parameters, parameters, parameters,
        @observe R[1, 1, :];
        series_labels = ("one", "two", "three", "four"),
        backend = :cairo,
        display_plot = false,
        fig_size = (420, 260),
        legend_position = :right,
        legend_overflow = :show_all
    )
    Makie.colorbuffer(complete.figure)
    complete_labels=[entry.label[] for entry in last(first(complete.legend.entrygroups[]))]
    @test length(complete_labels) == 4
    @test_throws ArgumentError Makie.plot(
        parameters,
        (R, 1, 1, Colon());
        backend = :cairo,
        display_plot = false,
        legend_overflow = :invalid
    )

    constant=Makie.plot(
        parameters,
        (R, 1, 1, Colon());
        series_labels = ("result",),
        backend = :cairo,
        display_plot = false
    )
    axis=only(constant.axes)
    limits=axis.finallimits[]
    @test all(isfinite, limits.origin)
    @test all(isfinite, limits.widths)
    first_entry=first(last(first(constant.legend.entrygroups[])))
    Makie.toggle_visibility!(first_entry)
    @test all(isfinite, axis.finallimits[].origin)
    Makie.toggle_visibility!(first_entry)
    @test all(isfinite, axis.finallimits[].origin)

    narrow=Makie.plot(
        parameters,
        @observe Z[1, 1, :];
        backend = :cairo,
        display_plot = false,
        fig_size = (600, 320)
    )
    @test length(narrow)==2
    foreach(p -> Makie.colorbuffer(p.figure),narrow)
    narrow_bounds=[only(p.axes).layoutobservables.computedbbox[] for p in narrow]
    @test all(bounds -> bounds.widths[1]>150,narrow_bounds)
    @test all(p -> only(p.axes).xlabelvisible[],narrow)

    tall=Makie.plot(
        parameters,
        @observe Z[1, 1, :];
        backend = :cairo,
        display_plot = false,
        fig_size = (600, 700)
    )
    @test length(tall)==2
    foreach(p -> Makie.colorbuffer(p.figure),tall)
    tall_bounds=[only(p.axes).layoutobservables.computedbbox[] for p in tall]
    @test all(p -> p.figure.scene.viewport[].widths[1]<p.figure.scene.viewport[].widths[2],tall)
    @test all(bounds -> bounds.widths[1]>150,tall_bounds)
    @test all(bounds -> bounds.widths[2]>220,tall_bounds)
    @test all(p -> only(p.axes).xlabelvisible[],tall)
end

@testitem "Makie addons / compact preview and material scheme geometry" tags=[:visual] setup=[
    UseNativePlotSupport, TestFixtures
] begin
    get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false")=="true"||
    error("set LINECABLEMODELS_TEST_PLOTTING=true to run the visual contract")
    using CairoMakie

    function inside(viewport, object)
        bounds=object.layoutobservables.computedbbox[]
        lower=bounds.origin
        upper=bounds.origin .+ bounds.widths
        viewport_lower=viewport.origin
        viewport_upper=viewport.origin .+ viewport.widths
        return all(lower .>= viewport_lower)&&all(upper .<= viewport_upper)
    end

    design=TestFixtures.coaxial_design()
    compact=preview(
        design;
        backend = :cairo,
        display_plot = false,
        open_export = false,
        size = (900, 350)
    )
    Makie.colorbuffer(compact.figure)
    viewport=compact.figure.scene.viewport[]
    axis=only(compact.axes)
    axis_bounds=axis.layoutobservables.computedbbox[]
    @test axis_bounds.widths[1] >= 160
    @test axis_bounds.widths[2] >= 160
    axis_viewport=axis.scene.viewport[]
    @test isapprox(axis_viewport.widths[1], axis_viewport.widths[2]; atol = 1)
    @test inside(viewport, axis)
    @test inside(viewport, compact.legend)
    @test length(compact.colorbars) == 3
    @test all(colorbar -> inside(viewport, colorbar), compact.colorbars)
    @test all(
        colorbar -> colorbar.layoutobservables.computedbbox[].widths[1] >= 135,
        compact.colorbars
    )

    collection=preview(
        fill(design, 4);
        layout = (2, 2),
        backend = :cairo,
        display_plot = false,
        controls = false,
        open_export = false,
        size = (1000, 850),
        colorbar_position = :bottom,
        colorbar_attributes = (; vertical = false)
    )
    Makie.colorbuffer(collection.figure)
    collection_viewport=collection.figure.scene.viewport[]
    collection_bounds=[axis.layoutobservables.computedbbox[] for axis in collection.axes]
    @test length(collection_bounds) == 4
    @test all(bounds -> bounds.widths[1] > 220, collection_bounds)
    @test all(bounds -> bounds.widths[2] > 200, collection_bounds)
    @test all(axis -> inside(collection_viewport, axis), collection.axes)
    @test all(
        colorbar -> inside(collection_viewport, colorbar), collection.colorbars
    )
    @test maximum(
        bounds.origin[2] + bounds.widths[2] for bounds in collection_bounds
    ) > 0.75collection_viewport.widths[2]

    wide=preview(
        design;
        backend = :cairo,
        display_plot = false,
        open_export = false,
        size = (1670, 965)
    )
    Makie.colorbuffer(wide.figure)
    wide_viewport=wide.figure.scene.viewport[]
    wide_axis=only(wide.axes)
    wide_axis_bounds=wide_axis.layoutobservables.computedbbox[]
    wide_legend_bounds=wide.legend.layoutobservables.computedbbox[]
    @test inside(wide_viewport, wide_axis)
    @test inside(wide_viewport, wide.legend)
    wide_viewport=only(wide.axes).scene.viewport[]
    @test isapprox(wide_viewport.widths[1], wide_viewport.widths[2]; atol = 1)
    @test 0 <
          wide_legend_bounds.origin[1] -
          (wide_axis_bounds.origin[1] + wide_axis_bounds.widths[1]) < 100

    reference=LineCableModels.show_material_scale(
        backend = :cairo,
        display_plot = false,
        open_export = false
    )
    Makie.colorbuffer(reference.figure)
    reference_viewport=reference.figure.scene.viewport[]
    @test length(reference.colorbars) == 3
    @test all(colorbar -> inside(reference_viewport, colorbar), reference.colorbars)
    @test all(
        colorbar -> colorbar.layoutobservables.computedbbox[].widths[1] > 500,
        reference.colorbars
    )
    vertical_positions=[colorbar.layoutobservables.computedbbox[].origin[2]
                        for colorbar in reference.colorbars]
    @test length(unique(vertical_positions)) == 3
end

@testitem "Makie / local native edits and flow identities preserve frames" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    options=(backend=:cairo,display_plot=false,open_export=false)
    raw=TestFixtures.two_conductor_results()
    p=LineCableModels.plot(raw;options...,ydata=(R,),layout=(1,1),series_labels=("one",)) |> first
    frame=only(p.axes).layoutobservables.computedbbox[].widths
    view=only(p.axes).targetlimits[]
    before=p.figure.scene.viewport[].widths
    p.legend.labelsize[]=30
    @test all(isapprox.(only(p.axes).layoutobservables.computedbbox[].widths,frame;atol=1))
    @test only(p.axes).targetlimits[]==view
    @test p.figure.scene.viewport[].widths[2]>before[2]
    widget=addwidget!((p,cell) -> Button(cell;label=repeat("Wide control ",30)),p,:wide)
    @test p.figure.scene.viewport[].widths[1]>=widget.layoutobservables.computedbbox[].widths[1]
    @test all(isapprox.(only(p.axes).layoutobservables.computedbbox[].widths,frame;atol=1))
    removewidget!(p,:wide)
    @test only(p.axes).targetlimits[]==view

    design=TestFixtures.coaxial_design()
    pages=preview(fill(design,4);options...,controls=false,layout=(1,2),
        panel_titles=Dict(1=>"First",4=>"Fourth"))
    @test first(pages).axes[1].title[]=="First"
    @test last(pages).axes[2].title[]=="Fourth"
    @test first(pages).axes[2].title[]==design.cable_id
    auto=preview(fill(design,3);options...,controls=false,panel_titles=i -> "Cable $i")
    identities=auto.addon_state.panel_page.coordinates
    resize!(auto.figure,1800,450)
    @test auto.addon_state.panel_page.coordinates==identities
    @test auto.addon_state.panel_page.dimensions==(1,3)
    resize!(auto.figure,450,1800)
    @test auto.addon_state.panel_page.dimensions==(3,1)
    @test [auto.addon_state.panel_data[i].axis.title[] for i in identities]==["Cable $i" for i in identities]

    canvas=LineCableModels.plotwindow(;title="Nested",options...,controls=false) do grid
        nested=GridLayout(grid[1,1])
        lines!(Axis(nested[1,1]),[1.,2.],[2.,4.];label="A")
        lines!(Axis(nested[2,1]),[1.,2.],[3.,6.];label="B")
    end
    plots=[copy(axis.scene.plots) for axis in canvas.axes]
    parent=[Makie.GridLayoutBase.gridcontent(axis).parent for axis in canvas.axes]
    for (id,label) in ((1,"A"),(2,"B"))
        legend=panellegend!(canvas,id;position=:right)
        @test only(last(only(legend.entrygroups[]))).label[]==label
    end
    @test [axis.scene.plots for axis in canvas.axes]==plots
    @test [Makie.GridLayoutBase.gridcontent(axis).parent for axis in canvas.axes]==parent
end
