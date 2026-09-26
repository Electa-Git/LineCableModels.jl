@testitem "Makie addons / responsive native legend and visible-series limits" tags=[:visual] setup=[
    UseNativePlotSupport, TestFixtures
] begin
    get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false")=="true"||
    error("set LINECABLEMODELS_TEST_PLOTTING=true to run the visual test")
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
        legend_cap = 0.5
    )
    Makie.colorbuffer(compact.figure)
    labels=[entry.label[] for entry in last(first(compact.legend.entrygroups[]))]
    @test !isempty(labels)
    @test length(labels) <= 4
    @test last(labels) == "(...)" || length(labels) == 4
    figurelegend!(compact; title = "Updated legend")
    resize!(compact.figure, 420, 140)
    @test length(last(only(compact.legend.entrygroups[])))<4
    @test first(only(compact.legend.entrygroups[]))=="Updated legend"
    # The edited title must fit the new width cap as well as the height.
    resize!(compact.figure, 680, 400)
    @test length(last(only(compact.legend.entrygroups[])))==4
    @test first(only(compact.legend.entrygroups[]))=="Updated legend"

    complete=Makie.plot(
        parameters, parameters, parameters, parameters,
        @observe R[1, 1, :];
        series_labels = ("one", "two", "three", "four"),
        backend = :cairo,
        display_plot = false,
        fig_size = (420, 260),
        legend_position = :right,
        legend_cap = 0.5
    )
    Makie.colorbuffer(complete.figure)
    complete_labels=[entry.label[] for entry in last(first(complete.legend.entrygroups[]))]
    @test length(complete_labels) == 4
    @test_throws ArgumentError Makie.plot(
        parameters,
        (R, 1, 1, Colon());
        backend = :cairo,
        display_plot = false,
        legend_cap = :invalid
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
    foreach(p->Makie.colorbuffer(p.figure), narrow)
    narrow_bounds=[only(p.axes).layoutobservables.computedbbox[] for p in narrow]
    @test all(bounds -> bounds.widths[1]>150, narrow_bounds)
    @test all(p -> only(p.axes).xlabelvisible[], narrow)

    tall=Makie.plot(
        parameters,
        @observe Z[1, 1, :];
        backend = :cairo,
        display_plot = false,
        fig_size = (600, 700)
    )
    @test length(tall)==2
    foreach(p->Makie.colorbuffer(p.figure), tall)
    tall_bounds=[only(p.axes).layoutobservables.computedbbox[] for p in tall]
    @test all(p -> p.figure.scene.viewport[].widths[1]<p.figure.scene.viewport[].widths[2], tall)
    @test all(bounds -> bounds.widths[1]>150, tall_bounds)
    @test all(bounds -> bounds.widths[2]>220, tall_bounds)
    @test all(p -> only(p.axes).xlabelvisible[], tall)
end

@testitem "Makie addons / numeric legend limits preserve data and native resizing" tags=[:visual] begin
    using CairoMakie
    f=[1.0, 10.0, 100.0]
    z=ComplexF64[complex(i+j+k, 1) for i in 1:6, j in 1:6, k in eachindex(f)]
    source=LineParameters(z, z .* 1e-6, f)
    options=(ydata = (R,), overlay = :coordinates, backend = :cairo, display_plot = false,
        open_export = false, fig_size = (680, 440))
    for position in (:bottom, :top, :left, :right, :inside)
        p=LineCableModels.plot(source; options..., legend_position = position)
        axis=only(p.axes)
        curves=filter(q->q isa Makie.Lines, axis.scene.plots)
        original=[copy(q[1][]) for q in curves]
        @test length(curves)==36
        entries=last(only(p.legend.entrygroups[]))
        @test 1<length(entries)<36
        @test last(entries).label[]=="(...)"
        dimension=position in (:left, :right) ? 1 : 2
        @test p.legend.layoutobservables.autosize[][dimension]<=0.5axis.scene.viewport[].widths[dimension]+1
        @test !isempty(Makie.colorbuffer(p.figure))
        views=[axis.targetlimits[] for axis in p.axes]
        resize!(p.figure, 1000, 900)
        @test Tuple(p.figure.scene.viewport[].widths)==(1000, 900)
        @test [axis.targetlimits[] for axis in p.axes]==views
        @test p.legend.layoutobservables.autosize[][dimension]<=0.5axis.scene.viewport[].widths[dimension]+1
        default_count=length(last(only(p.legend.entrygroups[])))
        figurelegend!(p; max_fraction = 1)
        @test length(last(only(p.legend.entrygroups[])))>=default_count
        @test p.legend.layoutobservables.autosize[][dimension]<=axis.scene.viewport[].widths[dimension]+1
        @test [q[1][] for q in curves]==original
        for invalid in (0, -0.1, 1.01, Inf, NaN, true, :show_all, :ellipsis)
            @test_throws ArgumentError figurelegend!(p; max_fraction = invalid)
        end
        @test p.addon_state.guides[(:legend, nothing)].max_fraction[]==1
        figurelegend!(p; max_fraction = 0.001)
        @test !p.legend.blockscene.visible[]
        figurelegend!(p; position = nothing)
        figurelegend!(p; max_fraction = 1)
        figurelegend!(p; position)
        @test p.legend.blockscene.visible[]
        empty!(p.figure)
    end
    for invalid in (0, -0.1, 1.01, Inf, NaN, true, :show_all, :ellipsis)
        @test_throws ArgumentError LineCableModels.plot(source; options..., legend_cap = invalid)
    end
    for sizing in ((fig_size = (680, 440),), (
        fig_size = (400, 300), figure = (size = (680, 440),)))
        p=LineCableModels.plot(source; options..., sizing...)
        @test p.addon_state.shell.reference_size==(680, 440)
        @test !isempty(Makie.colorbuffer(p.figure))
        empty!(p.figure)
    end
    fixed=LineCableModels.plot(source; options..., legend_attributes = (
        width = 1200, height = 800))
    axis=only(fixed.axes)
    @test fixed.legend.layoutobservables.computedbbox[].widths[1]<=axis.scene.viewport[].widths[1]+1
    @test fixed.legend.layoutobservables.computedbbox[].widths[2]<=0.5axis.scene.viewport[].widths[2]+1
    resize!(fixed.figure, 1000, 900)
    @test fixed.legend.layoutobservables.computedbbox[].widths[2]<=0.5axis.scene.viewport[].widths[2]+1
    @test !isempty(Makie.colorbuffer(fixed.figure))
    empty!(fixed.figure)
    panel=LineCableModels.plot(source; options..., legend_position = nothing,
        panel_legends = (1=>(position = :bottom, max_fraction = 0.5),))
    @test panel.panel_legends[1].layoutobservables.autosize[][2]<=0.5only(panel.axes).scene.viewport[].widths[2]+1
    @test last(last(only(panel.panel_legends[1].entrygroups[]))).label[]=="(...)"
    @test !isempty(Makie.colorbuffer(panel.figure))
    empty!(panel.figure)
end

@testitem "Makie addons / legend fitting ignores unchanged bounds and retains native edits" tags=[:visual] begin
    using CairoMakie
    ext=Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
    f=[1.0, 10.0, 100.0]
    z=ComplexF64[complex(i+j+k, 1) for i in 1:2, j in 1:2, k in eachindex(f)]
    p=LineCableModels.plot(
        LineParameters(z, z .* 1e-6, f); ydata = (R,), overlay = :coordinates,
        backend = :cairo, display_plot = false, controls = false, open_export = false)
    probe=only(filter(q->q isa Makie.Text && !q.visible[], p.legend.blockscene.plots))
    measurements=Ref(0)
    subscription=on(_->(measurements[]+=1), probe.text)
    publications=Ref(0)
    entry_subscription=on(_->(publications[]+=1), p.legend.entrygroups)
    bounds=p.addon_state.inside_bbox
    initial=bounds[]
    frames=[axis.scene.viewport[] for axis in p.axes]
    listeners=length(p.figure.scene.viewport.listeners)
    for _ in 1:3
        notify(bounds)
        ext._addon_edit_presentation!(() -> nothing, p)
    end
    @test measurements[]==0
    @test publications[]==0
    @test [axis.scene.viewport[] for axis in p.axes]==frames
    @test length(p.figure.scene.viewport.listeners)==listeners
    bounds[]=Makie.Rect2f(initial.origin .+ (1, 2), initial.widths)
    @test measurements[]==0
    bounds[]=Makie.Rect2f(initial.origin, initial.widths .* (0.6, 1))
    # Reflow reuses the cached measurements of unchanged labels and fonts.
    @test measurements[]==0
    before=measurements[]
    p.legend.labelsize[]+=2
    @test measurements[]>before
    before=measurements[]
    first(last(only(p.legend.entrygroups[]))).label[]="A longer coefficient label"
    @test measurements[]>before
    p.legend.nbanks[]=2
    resize!(p.figure, 900, 600)
    @test p.legend.nbanks[]==2
    @test !isempty(Makie.colorbuffer(p.figure))
    off(subscription)
    off(entry_subscription)
    empty!(p.figure)
end

@testitem "Makie addons / compact preview and material scheme geometry" tags=[:visual] setup=[
    UseNativePlotSupport, TestFixtures
] begin
    get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false")=="true"||
    error("set LINECABLEMODELS_TEST_PLOTTING=true to run the visual test")
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
        # Keep the compact right-column case explicit now that previews default
        # to a bottom strip.
        colorbar_position = :right,
        size = (900, 350)
    )
    Makie.colorbuffer(compact.figure)
    viewport=compact.figure.scene.viewport[]
    # size is the reference allocation; the physical frame and complete
    # right-side guides determine the fitted width.
    @test 0<viewport.widths[1]<900
    @test viewport.widths[2]>=350
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

    strip=preview(design; backend = :cairo, display_plot = false,
        open_export = false, size = (900, 350))
    Makie.colorbuffer(strip.figure)
    strip_viewport=strip.figure.scene.viewport[]
    @test all(colorbar -> inside(strip_viewport, colorbar), strip.colorbars)
    @test inside(strip_viewport, only(strip.axes))
    @test strip.addon_state.guides[(:colorbars, nothing)].position[]===:bottom

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

    pages=preview(fill(design, 3); backend = :cairo, display_plot = false,
        controls = false, layout = (1, 2))
    foreach(p->Makie.colorbuffer(p.figure), pages)
    frames=[axis.layoutobservables.computedbbox[].widths for p in pages for axis in p.axes]
    @test all(frame -> isapprox(frame, first(frames); atol = 2), frames)
    for p in pages
        viewport=p.figure.scene.viewport[]
        @test all(object -> inside(viewport, object), (p.axes..., p.colorbars...))
    end

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
    options=(backend = :cairo, display_plot = false, open_export = false)
    raw=TestFixtures.two_conductor_results()
    p=LineCableModels.plot(
        raw; options..., ydata = (R,), layout = (1, 1), series_labels = ("one",))|>first
    frame=only(p.axes).layoutobservables.computedbbox[].widths
    view=only(p.axes).targetlimits[]
    before=p.figure.scene.viewport[].widths
    p.legend.labelsize[]=30
    @test all(isapprox.(only(p.axes).layoutobservables.computedbbox[].widths, frame; atol = 1))
    @test only(p.axes).targetlimits[]==view
    @test p.figure.scene.viewport[].widths[2]>before[2]
    widget=addwidget!((p, cell)->Button(cell; label = repeat("Wide control ", 30)), p, :wide)
    @test p.figure.scene.viewport[].widths[1]>=widget.layoutobservables.computedbbox[].widths[1]
    @test all(isapprox.(only(p.axes).layoutobservables.computedbbox[].widths, frame; atol = 1))
    removewidget!(p, :wide)
    @test only(p.axes).targetlimits[]==view

    design=TestFixtures.coaxial_design()
    pages=preview(fill(design, 4); options..., controls = false, layout = (1, 2),
        panel_titles = Dict(1=>"First", 4=>"Fourth"))
    @test first(pages).axes[1].title[]=="First"
    @test last(pages).axes[2].title[]=="Fourth"
    @test first(pages).axes[2].title[]==design.cable_id
    auto=preview(fill(design, 3); options..., controls = false, panel_titles = i->"Cable $i")
    identities=auto.addon_state.panel_page.coordinates
    resize!(auto.figure, 1800, 450)
    @test auto.addon_state.panel_page.coordinates==identities
    @test auto.addon_state.panel_page.dimensions==(1, 3)
    resize!(auto.figure, 450, 1800)
    @test auto.addon_state.panel_page.dimensions==(3, 1)
    @test [auto.addon_state.panel_data[i].axis.title[] for i in identities]==["Cable $i"
                                                                              for i in identities]

    canvas=LineCableModels.plotwindow(; title = "Nested", options..., controls = false) do grid
        nested=GridLayout(grid[1, 1])
        lines!(Axis(nested[1, 1]), [1.0, 2.0], [2.0, 4.0]; label = "A")
        lines!(Axis(nested[2, 1]), [1.0, 2.0], [3.0, 6.0]; label = "B")
    end
    plots=[copy(axis.scene.plots) for axis in canvas.axes]
    parent=[Makie.GridLayoutBase.gridcontent(axis).parent for axis in canvas.axes]
    for (id, label) in ((1, "A"), (2, "B"))
        legend=panellegend!(canvas, id; position = :right)
        @test only(last(only(legend.entrygroups[]))).label[]==label
    end
    @test [axis.scene.plots for axis in canvas.axes]==plots
    @test [Makie.GridLayoutBase.gridcontent(axis).parent for axis in canvas.axes]==parent
end

@testitem "Makie addons / fitting preserves physical frames and native resize allocations" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    ext=Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
    design=TestFixtures.coaxial_design()
    p=preview(fill(design, 2); controls = false, display_plot = false,
        backend = :cairo, size = (1200, 900))
    @test p.figure.scene.viewport[].widths[2]<800
    @test p.addon_state.guide_gap==(8.0, 8.0, 24.0, 8.0)
    frames=[axis.scene.viewport[] for axis in p.axes]
    @test all(frame -> isapprox(frame.widths[1], frame.widths[2]; atol = 1), frames)
    @test all(all(isapprox.(axis.layoutobservables.computedbbox[].widths,
                  axis.scene.viewport[].widths; atol = 1)) for axis in p.axes)
    views=[axis.targetlimits[] for axis in p.axes]
    counts=[length(axis.scene.plots) for axis in p.axes]
    listeners=length(p.figure.scene.viewport.listeners)
    dimensions=Tuple(p.figure.scene.viewport[].widths)
    for _ in 1:5
        ext._addon_edit_presentation!(()->nothing, p)
        @test Tuple(p.figure.scene.viewport[].widths)==dimensions
        @test all(all(isapprox.(axis.scene.viewport[].widths, frame.widths; atol = 1))
        for (axis, frame) in zip(p.axes, frames))
    end
    @test [axis.targetlimits[] for axis in p.axes]==views
    @test [length(axis.scene.plots) for axis in p.axes]==counts
    @test length(p.figure.scene.viewport.listeners)==listeners
    resize!(p.figure, 1000, 800)
    @test Tuple(p.figure.scene.viewport[].widths)==(1000, 800)
    resized=[axis.scene.viewport[] for axis in p.axes]
    figuretitle!(p, "A title\nwith two lines")
    @test all(all(isapprox.(axis.scene.viewport[].widths, frame.widths; atol = 1))
    for (axis, frame) in zip(p.axes, resized))
    @test [axis.targetlimits[] for axis in p.axes]==views
    dimensions=Tuple(p.figure.scene.viewport[].widths)
    mktempdir() do directory
        for theme in (:default, :publication)
            @test isfile(export_svg(p; path = joinpath(directory, "$theme.svg"), theme, open_file = false))
            @test Tuple(p.figure.scene.viewport[].widths)==dimensions
            @test all(all(isapprox.(axis.scene.viewport[].widths, frame.widths; atol = 1))
            for (axis, frame) in zip(p.axes, resized))
        end
    end
    @test_throws ErrorException ext._addon_export_presentation!(p, :publication) do
        error("failing native writer")
    end
    @test Tuple(p.figure.scene.viewport[].widths)==dimensions
    @test [axis.targetlimits[] for axis in p.axes]==views
    @test all(all(isapprox.(axis.scene.viewport[].widths, frame.widths; atol = 1))
    for (axis, frame) in zip(p.axes, resized))
end

@testitem "Makie addons / native axis sizing constraints survive content fitting" tags=[:visual] begin
    using CairoMakie
    using LineCableModels
    ext=Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
    f=[1.0, 10.0, 100.0]
    z=reshape(complex.([1.0, 2.0, 3.0], [2.0, 3.0, 4.0]), 1, 1, :)
    raw=LineParameters(z, z .* 1e-6, f)
    for attributes in ((width = 700.0, height = 250.0),
        (width = Relative(0.5), height = Relative(0.5)), (alignmode = Outside(10),))
        p=LineCableModels.plot(raw; ydata = (R,), axis = attributes, layout = (1, 2),
            fig_size = (800, 500), controls = false, display_plot = false, backend = :cairo)
        axis=only(p.axes)
        initial=axis.scene.viewport[]
        size=Tuple(p.figure.scene.viewport[].widths)
        @test all(initial.origin .>= 0)
        @test all(initial.origin+initial.widths .<= size)
        haskey(attributes, :width) && attributes.width isa Real &&
            @test Tuple(initial.widths)==(700, 250)
        for _ in 1:3
            ext._addon_edit_presentation!(() -> nothing, p)
            @test all(isapprox.(axis.scene.viewport[].widths, initial.widths; atol = 1))
            @test Tuple(p.figure.scene.viewport[].widths)==size
        end
    end
end
