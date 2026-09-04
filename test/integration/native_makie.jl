@testitem "Makie addons / native figures, layouts, docks, widgets, and export" tags=[:visual] setup=[
    NativePlotTestSupport, UseNativePlotSupport, TestFixtures
] begin
    get(ENV,
        "LINECABLEMODELS_TEST_PLOTTING",
        "false")=="true"||
    error("set LINECABLEMODELS_TEST_PLOTTING=true to run the visual contract")
    using CairoMakie

    legend_labels(legend) = [entry.label[] for entry in last(first(legend.entrygroups[]))]

    frequency=[50.0, 100.0, 500.0]
    parameters=TestFixtures.two_conductor_results(; frequencies = frequency)

    # With no selection, Z is processed before Y and their physically meaningful
    # components are rendered as four matrix dashboards: R, X, G, then B.
    automatic=Makie.plot(
        parameters;
        backend = :cairo,
        display_plot = false,
        open_export = false
    )
    @test automatic isa Vector{UIPlot}
    @test length(automatic) == 4
    @test all(page -> length(page.axes) == 4, automatic)
    @test all(page -> page.legend === nothing, automatic)
    @test [first(page.axes).title[] for page in automatic] == [
        "Self series resistance — conductor 1",
        "Self series reactance — conductor 1",
        "Self shunt conductance — conductor 1",
        "Self shunt susceptance — conductor 1"
    ]
    automatic_page=first(automatic)
    @test automatic_page.figure isa Makie.Figure
    @test !hasproperty(automatic_page, :context)
    @test !hasproperty(automatic_page, :panels)
    @test :legend ∉ keys(automatic_page.controls)

    self_impedance_request=@observe Z[1, 1, :]
    side_by_side=Makie.plot(
        parameters,
        self_impedance_request;
        backend = :cairo,
        display_plot = false,
        controls = false
    )
    @test side_by_side isa UIPlot
    @test [axis.title[] for axis in side_by_side.axes] == [
        "Self series resistance — conductor 1",
        "Self series reactance — conductor 1"
    ]
    @test side_by_side.legend === nothing

    stacked_complex=Makie.plot(
        parameters,
        self_impedance_request;
        backend = :cairo,
        display_plot = false,
        controls = false,
        layout = (2, 1),
        panel_titles = ("Self resistance", "Self reactance")
    )
    @test stacked_complex isa UIPlot
    @test [axis.title[] for axis in stacked_complex.axes] ==
          ["Self resistance", "Self reactance"]
    @test Set(keys(stacked_complex.addon_state.panel_data)) ==
          Set(((1, 1), (2, 1)))

    individual=Makie.plot(
        parameters,
        self_impedance_request;
        backend = :cairo,
        display_plot = false,
        controls = false,
        layout = (1, 1)
    )
    @test individual isa Vector{UIPlot}
    @test length(individual) == 2
    @test all(page -> length(page.axes) == 1, individual)

    paired_all=Makie.plot(
        parameters;
        backend = :cairo,
        display_plot = false,
        controls = false,
        layout = (1, 2)
    )
    @test length(paired_all) == 8
    @test all(page -> length(page.axes) == 2, paired_all)

    stacked_all=Makie.plot(
        parameters;
        backend = :cairo,
        display_plot = false,
        controls = false,
        layout = (2, 1)
    )
    @test length(stacked_all) == 8
    @test all(
        page -> Set(keys(page.addon_state.panel_data)) ==
                Set(((1, 1), (2, 1))), stacked_all)

    one_quantity=Makie.plot(
        parameters,
        (R,);
        backend = :cairo,
        display_plot = false,
        controls = false
    )
    @test one_quantity isa UIPlot
    @test length(one_quantity.axes) == 4
    @test isempty(one_quantity.controls)
    @test one_quantity.legend === nothing

    initially_logarithmic=Makie.plot(
        parameters,
        (R, 1, 1, Colon());
        backend = :cairo,
        display_plot = false,
        controls = false,
        xscale = :log10
    )
    Makie.colorbuffer(initially_logarithmic.figure)
    @test only(initially_logarithmic.axes).xscale[] === Makie.log10
    @test only(initially_logarithmic.axes).finallimits[].origin[1] > 0

    compared_self=Makie.plot(
        parameters,
        parameters,
        self_impedance_request;
        series_labels = ("reference", "candidate"),
        backend = :cairo,
        display_plot = false,
        controls = false
    )
    @test compared_self isa UIPlot
    @test length(compared_self.axes) == 2
    @test legend_labels(compared_self.legend) == ["reference", "candidate"]
    @test Set(keys(compared_self.addon_state.groups)) == Set((:result_1, :result_2))

    three_sources=Makie.plot(
        parameters,
        parameters,
        parameters,
        self_impedance_request;
        series_labels = ("reference", "candidate A", "candidate B"),
        backend = :cairo,
        display_plot = false,
        controls = false
    )
    @test legend_labels(three_sources.legend) ==
          ["reference", "candidate A", "candidate B"]

    different_frequency=TestFixtures.two_conductor_results(
        ;
        frequencies = [40.0, 80.0, 160.0, 320.0])
    asynchronous=Makie.plot(
        parameters,
        different_frequency,
        self_impedance_request;
        series_labels = ("three samples", "four samples"),
        backend = :cairo,
        display_plot = false,
        controls = false
    )
    @test asynchronous isa UIPlot
    @test legend_labels(asynchronous.legend) == ["three samples", "four samples"]

    combined=Makie.plot(
        (; reference = parameters, candidate = parameters),
        (R,);
        backend = :cairo,
        display_plot = false,
        layout = (2, 2),
        legend_position = :bottom,
        legend_attributes = (; orientation = :horizontal),
        legend_overflow = :show_all
    )
    @test combined isa UIPlot
    @test length(combined.axes) == 4
    @test combined.legend.orientation[] == :horizontal
    @test legend_labels(combined.legend) == ["reference", "candidate"]

    for (position,
        orientation) in (
        (:left, :vertical), (:right, :vertical),
        (:top, :horizontal), (:bottom, :horizontal),
        ((1, 3), :horizontal)
    )
        placed=Makie.plot(
            parameters,
            parameters,
            (R,);
            series_labels = ("reference", "candidate"),
            backend = :cairo,
            display_plot = false,
            legend_position = position,
            legend_attributes = (; orientation),
            legend_overflow = :show_all
        )
        @test placed.legend.orientation[] == orientation
    end

    titled_inside=Makie.plot(
        parameters,
        parameters,
        self_impedance_request;
        series_labels = ("reference", "candidate"),
        backend = :cairo,
        display_plot = false,
        controls = false,
        figure_title = "Cable impedance",
        panel_titles = ("Resistance", "Reactance"),
        legend_position = :inside,
        legend_anchor = :lt,
        legend_title = "Result set",
        legend_overflow = :show_all
    )
    Makie.colorbuffer(titled_inside.figure)
    @test titled_inside.title.text[] == "Cable impedance"
    @test [axis.title[] for axis in titled_inside.axes] ==
          ["Resistance", "Reactance"]
    @test first(only(titled_inside.legend.entrygroups[])) == "Result set"
    @test titled_inside.legend.halign[] == :left
    @test titled_inside.legend.valign[] == :top
    @test legend_labels(titled_inside.legend) == ["reference", "candidate"]
    canvas_bounds=titled_inside.addon_state.shell.canvas.layoutobservables.computedbbox[]
    legend_bounds=titled_inside.legend.layoutobservables.computedbbox[]
    @test legend_bounds.origin[1] >= canvas_bounds.origin[1]
    @test legend_bounds.origin[2] >= canvas_bounds.origin[2]
    @test legend_bounds.origin[1] + legend_bounds.widths[1] <=
          canvas_bounds.origin[1] + canvas_bounds.widths[1]
    @test legend_bounds.origin[2] + legend_bounds.widths[2] <=
          canvas_bounds.origin[2] + canvas_bounds.widths[2]

    figuretitle!(titled_inside, "Updated impedance"; fontsize = 20)
    paneltitle!(titled_inside, (1, 2), "Updated reactance")
    @test titled_inside.title.text[] == "Updated impedance"
    @test titled_inside.title.fontsize[] == 20
    @test titled_inside.axes[2].title[] == "Updated reactance"

    figurelegend!(
        titled_inside;
        position = :right,
        title = "Renamed results",
        legend_labels = Dict(:result_1=>"baseline"),
        overflow = :show_all
    )
    @test first(only(titled_inside.legend.entrygroups[])) == "Renamed results"
    @test legend_labels(titled_inside.legend) == ["baseline", "candidate"]
    figurelegend!(
        titled_inside;
        position = :inside,
        anchor = :rb,
        overflow = :show_all
    )
    @test titled_inside.legend.halign[] == :right
    @test titled_inside.legend.valign[] == :bottom

    panellegend!(
        titled_inside,
        (1, 1);
        position = :inside,
        anchor = :lb,
        title = "Resistance results",
        legend_labels = ("base R", "candidate R"),
        overflow = :show_all
    )
    Makie.colorbuffer(titled_inside.figure)
    panel_inside=titled_inside.panel_legends[(1, 1)]
    @test first(only(panel_inside.entrygroups[])) == "Resistance results"
    @test panel_inside.halign[] == :left
    @test panel_inside.valign[] == :bottom
    @test legend_labels(panel_inside) == ["base R", "candidate R"]
    first_panel_entry=first(last(first(panel_inside.entrygroups[])))
    Makie.toggle_visibility!(first_panel_entry)
    @test all(plot -> !plot.visible[],
        titled_inside.addon_state.panel_data[(1, 1)].groups[:result_1])
    @test titled_inside.addon_state.shell.status[] ==
          "Axis limits fitted to visible series"
    Makie.toggle_visibility!(first_panel_entry)
    panellegend!(titled_inside, (1, 1); position = :right, overflow = :show_all)
    @test titled_inside.panel_legends[(1, 1)].orientation[] == :vertical

    paneltitle!(titled_inside, (1, 2), nothing)
    @test titled_inside.axes[2].title[] == ""
    figuretitle!(titled_inside, nothing)
    @test titled_inside.title === nothing

    per_page_titles=Makie.plot(
        parameters,
        (R, X);
        backend = :cairo,
        display_plot = false,
        controls = false,
        figure_title = ("Resistance dashboard", "Reactance dashboard")
    )
    @test [page.title.text[] for page in per_page_titles] ==
          ["Resistance dashboard", "Reactance dashboard"]

    @test_throws DimensionMismatch Makie.plot(
        parameters,
        self_impedance_request;
        backend = :cairo,
        display_plot = false,
        panel_titles = ("only one",)
    )
    @test_throws ArgumentError Makie.plot(
        parameters,
        (R,);
        backend = :cairo,
        display_plot = false,
        layout = (1, 3)
    )

    compared=Makie.plot(
        (; reference = parameters, candidate = parameters),
        (R, X, G, B);
        backend = :cairo,
        display_plot = false,
        legend_position = :top,
        legend_attributes = (; orientation = :horizontal)
    )
    @test length(compared) == 4
    @test all(plot -> length(plot.axes) == 4, compared)
    @test all(plot -> legend_labels(plot.legend) == ["reference", "candidate"], compared)

    first_axis=first(automatic_page.axes)
    first_axis.title[]="caller-owned title"
    extra=lines!(first_axis, frequency, fill(1.0e-4, 3); color = :red)
    @test first_axis.title[] == "caller-owned title"
    @test extra in first_axis.scene.plots
    automatic_page.controls[:xlog].active[]=true
    @test all(axis -> axis.xscale[] === Makie.log10, automatic_page.axes)
    automatic_page.controls[:xlog].active[]=false
    @test all(axis -> axis.xscale[] === Makie.identity, automatic_page.axes)
    @test occursin("10", sprint(show, first_axis.ylabel[]))

    library=CablesLibrary()
    load!(library;
        file_name = joinpath(
            pkgdir(LineCableModels), "test", "fixtures", "data", "mv_cable_design.json"
        )
    )
    design=first(values(library.data))
    cable=preview(
        design;
        backend = :cairo,
        display_plot = false,
        open_export = false
    )
    @test cable isa UIPlot
    @test length(cable.axes) == 1
    @test cable.legend !== nothing
    @test length(cable.colorbars) == 3
    @test all(colorbar -> !colorbar.vertical[], cable.colorbars)
    @test Set(keys(cable.controls)) == Set((:reset, :export_svg, :legend))

    grouped_cable=preview(
        design;
        backend = :cairo,
        display_plot = false,
        controls = false,
        display_colorbars = false,
        legend_group = region->startswith(
            String(region.source.tag), "core_wire_") ?
                               :stranded_core : region.source.tag,
        legend_labels = Dict(:stranded_core=>"Stranded core")
    )
    @test count(==("Stranded core"), legend_labels(grouped_cable.legend)) == 1

    late_legend_cable=preview(
        design;
        backend = :cairo,
        display_plot = false,
        controls = false,
        display_legend = false
    )
    @test late_legend_cable.legend === nothing
    figurelegend!(late_legend_cable; position = :right, overflow = :show_all)
    Makie.colorbuffer(late_legend_cable.figure)
    late_legend_bounds=late_legend_cable.legend.layoutobservables.computedbbox[]
    late_colorbar_bounds=[colorbar.layoutobservables.computedbbox[]
                          for colorbar in late_legend_cable.colorbars]
    @test late_legend_bounds.origin[2] > maximum(
        bounds.origin[2] + bounds.widths[2] for bounds in late_colorbar_bounds
    )
    figurelegend!(
        late_legend_cable;
        position = :top,
        orientation = :horizontal,
        overflow = :show_all
    )
    Makie.colorbuffer(late_legend_cable.figure)
    @test late_legend_cable.legend.orientation[] == :horizontal
    @test length(late_legend_cable.colorbars) == 3

    ranges=LineCableModels.DataModel.material_property_ranges(design)
    rho_scheme=LineCableModels.materialcolors(:rho, ranges.rho)
    native_scale_figure=Figure(; size = (500, 300))
    native_scale_axis=Axis(native_scale_figure[1, 1])
    lines!(native_scale_axis, 1:3, 1:3)
    native_scale=materialscale!(
        native_scale_figure[2, 1],
        rho_scheme;
        vertical = false
    )
    Makie.colorbuffer(native_scale_figure)
    @test !native_scale.vertical[]
    @test native_scale.layoutobservables.computedbbox[].widths[1] > 300

    moved_docks=preview(
        design;
        backend = :cairo,
        display_plot = false,
        legend_position = :top,
        legend_attributes = (; orientation = :horizontal),
        colorbar_position = :left,
        colorbar_attributes = (; vertical = true)
    )
    @test moved_docks.legend.orientation[] == :horizontal
    @test all(colorbar -> colorbar.vertical[], moved_docks.colorbars)

    collection=preview(
        [design, design, design, design];
        layout = (2, 2),
        backend = :cairo,
        display_plot = false,
        figure_title = "Design family",
        panel_titles = ("A", "B", "C", "D")
    )
    @test length(collection.axes) == 4
    @test length(collection.colorbars) == 3
    @test collection.title.text[] == "Design family"
    @test [axis.title[] for axis in collection.axes] == ["A", "B", "C", "D"]

    system=TestFixtures.three_phase_system()
    earth=EarthModel(100.0, 10.0, 1.0)
    system_plot=preview(
        system;
        earth_model = earth,
        backend = :cairo,
        display_plot = false
    )
    @test length(system_plot.axes) == 1
    @test system_plot.legend !== nothing
    @test length(system_plot.colorbars) == 3

    scale=LineCableModels.show_material_scale(
        backend = :cairo,
        display_plot = false,
        figure_title = "Material schemes",
        colorbar_position = :top,
        colorbar_attributes = (; vertical = false)
    )
    @test isempty(scale.axes)
    @test length(scale.colorbars) == 3
    @test scale.title.text[] == "Material schemes"

    native_grid=LineCableModels.plotwindow(
        ; title = "native 2x2", size = (700, 500), layout = (2, 2),
        figure_title = "Native dashboard",
        backend = :cairo, display_plot = false, open_export = false
    ) do grid
        for row in 1:2, column in 1:2

            axis=Axis(grid[row, column]; title = "$row,$column")
            lines!(axis, 1:3, (1:3) .* (row+column))
        end
    end
    @test length(native_grid.axes) == 4
    @test Set(keys(native_grid.controls)) == Set((:reset, :export_svg))
    @test native_grid.title.text[] == "Native dashboard"

    report_plot=report(
        TableReportDefinition(
            ((R, Colon(), Colon(), Colon()),);
            illustration = true,
            plot_options = (;
                backend = :cairo,
                display_plot = false,
                open_export = false
            )
        ),
        parameters
    ).illustration
    @test report_plot isa UIPlot
    @test length(report_plot.axes) == 1
    @test occursin("Frequency", sprint(show, only(report_plot.axes).xlabel[]))
    @test length(report_plot.addon_state.groups) == 4

    extension=Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
    publication_snapshot=Pair{Any, Any}[]
    publication_layout=Ref{Any}(nothing)
    root=combined.figure.layout
    row_sizes_before=copy(root.rowsizes)
    row_gap_before=root.default_rowgap
    background_before=combined.figure.scene.backgroundcolor[]
    font_before=combined.figure.scene.theme[:fonts][:regular][]
    reset_visible_before=combined.controls[:reset].blockscene.visible[]
    try
        publication_layout[]=extension._native_publication_snapshot!(
            publication_snapshot,
            combined,
            :publication
        )
        @test root.rowsizes[1] == Makie.Fixed(0)
        @test root.rowsizes[3] == Makie.Fixed(0)
        @test !combined.controls[:reset].blockscene.visible[]
        publication_font=Makie.to_font(
            combined.figure.scene.theme[:fonts],
            :regular
        )
        @test occursin("NewComputerModern", sprint(show, publication_font))
    finally
        extension._native_restore_interactive_chrome!(publication_layout[])
        extension._native_restore_snapshot!(publication_snapshot)
    end
    @test root.rowsizes == row_sizes_before
    @test root.default_rowgap == row_gap_before
    @test combined.figure.scene.backgroundcolor[] == background_before
    @test combined.figure.scene.theme[:fonts][:regular][] === font_before
    @test combined.controls[:reset].blockscene.visible[] == reset_visible_before

    mktempdir() do directory
        path=export_svg(
            combined;
            path = joinpath(directory, "current.svg"),
            theme = :publication,
            open_file = false
        )
        @test isfile(path)
        @test filesize(path) > 0
        @test combined.legend.orientation[] == :horizontal
        @test extension.current_backend_symbol() == :cairo
        @test root.rowsizes == row_sizes_before
        @test root.default_rowgap == row_gap_before
        @test combined.controls[:reset].blockscene.visible[] == reset_visible_before
    end
end
