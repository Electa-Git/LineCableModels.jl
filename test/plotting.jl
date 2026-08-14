@testitem "PlotBuilder Cairo: golden renders, callbacks, and SVG export" setup = [defaults] begin
    if get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false") != "true"
        @test_skip "Cairo plotting gate is disabled"
    else
        using CairoMakie
        using Measurements: measurement

        function pixel_error(current, reference)
            size(current) == size(reference) || return Inf
            channel_error = abs.(Makie.red.(current) .- Makie.red.(reference)) .+
                            abs.(Makie.green.(current) .- Makie.green.(reference)) .+
                            abs.(Makie.blue.(current) .- Makie.blue.(reference))
            return sum(channel_error) / (3length(channel_error))
        end

        function test_golden(handle::UIPlot, name; tolerance = 0.015)
            reference_path = joinpath(pkgdir(LineCableModels), "test", "reference", "$name.png")
            @test isfile(reference_path)
            reference = CairoMakie.FileIO.load(reference_path)
            current = Makie.colorbuffer(handle.figure)
            alternatives = (
                reference,
                reverse(reference; dims = 1),
                reverse(reference; dims = 2),
                reverse(reverse(reference; dims = 1); dims = 2)
            )
            @test minimum(pixel_error(current, candidate) for candidate in alternatives) <
                  tolerance
            return nothing
        end

        frequency = [50.0, 100.0, 500.0]
        omega = reshape(2π .* frequency, 1, 1, :)
        resistance_values = reshape([1.0, 0.2, 0.2, 2.0], 2, 2, 1) .*
                            ones(1, 1, length(frequency)) .* 1.0e-4
        inductance_values = fill(2.0e-7, 2, 2, length(frequency))
        conductance_values = fill(3.0e-9, 2, 2, length(frequency))
        capacitance_values = fill(4.0e-10, 2, 2, length(frequency))
        parameters = LineParameters(
            complex.(resistance_values, inductance_values .* omega),
            complex.(conductance_values, capacitance_values .* omega),
            frequency
        )

        rlcg = Makie.plot(
            parameters;
            mode = :RLCG,
            backend = :cairo,
            display_plot = false,
            open_export = false
        )
        cartesian = Makie.plot(
            parameters;
            mode = :ZY,
            coord = :cart,
            backend = :cairo,
            display_plot = false
        )
        polar = Makie.plot(
            parameters;
            mode = :ZY,
            coord = :polar,
            backend = :cairo,
            display_plot = false
        )
        @test rlcg isa Vector{UIPlot}
        @test cartesian isa Vector{UIPlot}
        @test polar isa Vector{UIPlot}
        @test length(rlcg) == 4
        @test length(cartesian) == 4
        @test length(polar) == 4
        @test rlcg[1].context !== rlcg[2].context
        @test rlcg[1].context.status !== rlcg[2].context.status
        series_plots = Makie.plot(
            parameters.Z,
            frequency;
            mode = :ZY,
            backend = :cairo,
            display_plot = false
        )
        shunt_plots = Makie.plot(
            parameters.Y,
            frequency;
            mode = :RLCG,
            backend = :cairo,
            display_plot = false
        )
        @test series_plots isa Vector{UIPlot}
        @test shunt_plots isa Vector{UIPlot}
        @test length(series_plots) == 2
        @test length(shunt_plots) == 2
        test_golden(first(rlcg), "line_rlcg")
        test_golden(first(cartesian), "line_zy_cartesian")
        test_golden(first(polar), "line_zy_polar")

        measurement_parameters = LineParameters(
            complex.(measurement.(resistance_values, resistance_values .* 0.05),
                measurement.(inductance_values .* omega, inductance_values .* omega .*
                                                         0.05)),
            complex.(measurement.(conductance_values, conductance_values .* 0.05),
                measurement.(capacitance_values .* omega, capacitance_values .* omega .*
                                                          0.05)),
            frequency
        )
        measurement_plots = Makie.plot(
            measurement_parameters;
            mode = :RLCG,
            backend = :cairo,
            display_plot = false
        )
        measurement_plot = first(measurement_plots)
        @test length(only(measurement_plot.panels).plots) >
              length(only(measurement_plot.page.views).series)
        test_golden(measurement_plot, "line_measurements")

        measurement_conductance = measurement_plots[3]
        measurement_conductance_axis = only(measurement_conductance.panels).axis
        measurement_limits = measurement_conductance_axis.finallimits[]
        measurement_ymin = measurement_limits.origin[2]
        measurement_ymax = measurement_ymin + measurement_limits.widths[2]
        @test measurement_ymin < 2.85e-6
        @test measurement_ymax > 3.15e-6
        measurement_conductance.controls[:ylog].active[] = true
        log_measurement_limits = measurement_conductance_axis.finallimits[]
        log_measurement_ymin = log_measurement_limits.origin[2]
        log_measurement_ymax = log_measurement_ymin + log_measurement_limits.widths[2]
        @test log_measurement_ymin < 2.85e-6
        @test log_measurement_ymax > 3.15e-6

        handle = first(rlcg)
        @test handle.figure.scene.backgroundcolor[] == Makie.to_color(:grey90)
        padding = handle.figure.layout.alignmode[].padding
        @test (padding.left, padding.right, padding.bottom, padding.top) ==
              (20.0, 20.0, 28.0, 28.0)
        @test handle.figure.layout.colsizes[2] == Makie.Auto(true)
        @test handle.figure.layout.addedrowgaps == fill(Makie.Fixed(6), 2)
        @test handle.figure.layout.addedcolgaps == [Makie.Fixed(6)]
        @test handle.controls[:reset].buttoncolor[] == Makie.RGBf(0.94, 0.94, 0.94)
        @test occursin("\\ue5d5", sprint(show, handle.controls[:reset].label[]))
        @test occursin("\\ue161", sprint(show, handle.controls[:export_svg].label[]))
        line_axis = only(handle.panels).axis
        initial_line_limits = line_axis.finallimits[]
        initial_xlabel = line_axis.xlabel[]
        handle.controls[:xlog].active[] = true
        @test line_axis.xscale[] === Makie.log10
        @test line_axis.xlabel[] == "Frequency [Hz]"
        log_line_limits = line_axis.finallimits[]
        @test log_line_limits.origin[1] < minimum(frequency)
        @test log_line_limits.origin[1] + log_line_limits.widths[1] > maximum(frequency)
        handle.controls[:xlog].active[] = false
        @test line_axis.xscale[] === Makie.identity
        @test line_axis.xlabel[] == initial_xlabel
        @test line_axis.finallimits[] == initial_line_limits
        handle.controls[:xlog].active[] = true
        handle.controls[:ylog].active[] = true
        @test line_axis.yscale[] === Makie.log10
        handle.controls[:reset].clicks[] += 1
        @test handle.context.status[] == "Axis limits reset"
        legend = handle.controls[:legend]
        first_entry = first(last(first(legend.entrygroups[])))
        Makie.toggle_visibility!(first_entry)
        @test any(plot_object -> !plot_object.visible[], only(handle.panels).plots)
        Makie.xlims!(only(handle.panels).axis, 100.0, 300.0)
        ui_components = Base.get_extension(
            LineCableModels,
            :LineCableModelsMakieExt
        ).UIComponents
        default_export_theme = ui_components._theme(
            export_mode = true,
            export_theme = :default
        )
        publication_export_theme = ui_components._theme(
            export_mode = true,
            export_theme = :publication
        )
        @test !haskey(default_export_theme[:fonts], :regular)
        @test haskey(publication_export_theme[:fonts], :regular)
        @test haskey(publication_export_theme[:fonts], :italic)
        @test publication_export_theme[:Axis][:titlesize][] == 15
        @test publication_export_theme[:Axis][:xticklabelsize][] == 14
        current_page = ui_components._current_page(handle)
        current_view = only(current_page.views)
        @test current_view.xaxis.scale === :log10
        @test current_view.yaxis.scale === :log10
        @test any(series -> !series.attributes.visible, current_view.series)
        @test haskey(current_view.attributes, :limits)
        @test collect(current_view.attributes.limits[1]) ≈ [100.0, 300.0]
        Makie.toggle_visibility!(first_entry)
        @test all(plot_object -> plot_object.visible[], only(handle.panels).plots)

        susceptance_handle = last(cartesian)
        susceptance_axis = only(susceptance_handle.panels).axis
        @test susceptance_axis.ylabel[] isa Makie.RichText
        @test occursin("−3", sprint(show, susceptance_axis.ylabel[]))
        initial_susceptance_limits = susceptance_axis.finallimits[]
        susceptance_handle.controls[:ylog].active[] = true
        @test susceptance_axis.yscale[] === Makie.log10
        @test susceptance_axis.ytickformat[] === Makie.automatic
        @test susceptance_axis.ylabel[] == "Capacitive susceptance [S/km]"
        limits = susceptance_axis.finallimits[]
        ymin = limits.origin[2]
        ymax = ymin + limits.widths[2]
        tick_values, tick_labels = Makie.get_ticks(
            susceptance_axis.yticks[],
            susceptance_axis.yscale[],
            susceptance_axis.ytickformat[],
            ymin,
            ymax
        )
        @test length(tick_values) in 1:4
        @test all(isinteger, log10.(tick_values))
        @test all(isone, round.(diff(log10.(tick_values)); digits = 8))
        @test all(label -> label isa Makie.RichText, tick_labels)
        susceptance_handle.controls[:ylog].active[] = false
        @test susceptance_axis.yscale[] === Makie.identity
        @test susceptance_axis.ylabel[] isa Makie.RichText
        @test susceptance_axis.finallimits[] == initial_susceptance_limits

        conductance_handle = rlcg[3]
        conductance_axis = only(conductance_handle.panels).axis
        conductance_handle.controls[:ylog].active[] = true
        conductance_limits = conductance_axis.finallimits[]
        @test conductance_limits.origin[2] == 1.0e-6
        @test conductance_limits.origin[2] + conductance_limits.widths[2] == 1.0e-5
        conductance_ticks, conductance_labels = Makie.get_ticks(
            conductance_axis.yticks[],
            conductance_axis.yscale[],
            conductance_axis.ytickformat[],
            1.0e-6,
            1.0e-5
        )
        @test conductance_ticks == [1.0e-6, 1.0e-5]
        @test all(label -> label isa Makie.RichText, conductance_labels)

        capacitance_axis = only(rlcg[4].panels).axis
        capacitance_limits = capacitance_axis.finallimits[]
        @test capacitance_limits.origin[2] == 0.38
        @test capacitance_limits.origin[2] + capacitance_limits.widths[2] ≈ 0.42

        mktempdir() do directory
            svg_path = joinpath(directory, "line parameters.svg")
            exported = export_svg(handle; path = svg_path, open_file = false)
            @test exported == abspath(svg_path)
            @test handle.context.status[] == "Saved SVG to $exported"
            @test filesize(exported) > 100
            @test occursin("<svg", read(exported, String))
            @test !occursin("Export SVG", read(exported, String))
            @test_throws ArgumentError export_svg(handle; path = svg_path, open_file = false)
            @test_throws ArgumentError export_svg(
                handle;
                path = joinpath(directory, "not-an-svg.png"),
                open_file = false
            )
            @test_throws ArgumentError export_svg(
                handle;
                path = joinpath(directory, "unsupported-theme.svg"),
                theme = :unsupported,
                open_file = false
            )
            @test nameof(Makie.current_backend()) === :CairoMakie

            cd(directory) do
                first_default = export_svg(handle; open_file = false)
                second_default = export_svg(handle; open_file = false)
                @test first_default != second_default
                @test dirname(first_default) == directory
                @test occursin(
                    r"^series_resistance_\d{8}_\d{6}(?:_\d+)?\.svg$",
                    basename(first_default)
                )
                @test occursin("rgb(100%, 100%, 100%)", read(first_default, String))
                @test filesize(second_default) > 100
            end

            publication_path = joinpath(directory, "publication.svg")
            publication = export_svg(
                handle;
                path = publication_path,
                theme = :publication,
                open_file = false
            )
            @test publication == publication_path
            @test filesize(publication) > 100
            @test read(publication) != read(exported)
        end

        fallback_export = cd(pkgdir(LineCableModels)) do
            export_svg(handle; open_file = false)
        end
        @test dirname(fallback_export) ==
              joinpath(tempdir(), "linecablemodels-exports")
        @test filesize(fallback_export) > 100
        rm(fallback_export)

        summary = SampleSummary([1.0, 2.0, 3.0, 4.0])
        distribution_model = HistogramPDF([1.0, 3.0, 5.0], [0.25, 0.25])
        mc_result = CableConstantsMC(
            CableConstants(summary, summary, summary),
            CableConstants(
                [1.0, 2.0, 3.0, 4.0],
                [1.0, 2.0, 3.0, 4.0],
                [1.0, 2.0, 3.0, 4.0]
            ),
            CableConstants(distribution_model, distribution_model, distribution_model),
            CableConstants(
                measurement(2.5, 0.0),
                measurement(2.5, 0.0),
                measurement(2.5, 0.0)
            ),
            4,
            0.95
        )
        for mode in (:hist, :pdf, :ecdf, :qq)
            mc_plot = Makie.plot(
                mc_result,
                :R;
                mode,
                data = :both,
                backend = :cairo,
                display_plot = false
            )
            @test mc_plot isa UIPlot
            @test !haskey(mc_plot.controls, :xlog)
            @test !haskey(mc_plot.controls, :ylog)
            test_golden(mc_plot, "mc_$mode")
        end

        line_samples = reshape(collect(1.0:12.0), 1, 1, 3, 4)
        summarize(values) = map(
            index -> SampleSummary(view(values, index.I..., :)),
            CartesianIndices(size(values)[1:3])
        )
        line_statistics = RLCG(
            summarize(line_samples),
            summarize(line_samples .* 1.0e-3),
            summarize(line_samples .* 1.0e-6),
            summarize(line_samples .* 1.0e-4)
        )
        line_distributions = RLCG(
            fill(distribution_model, 1, 1, 3),
            fill(distribution_model, 1, 1, 3),
            fill(distribution_model, 1, 1, 3),
            fill(distribution_model, 1, 1, 3)
        )
        line_mc = LineParametersMC(
            line_statistics,
            RLCG(
                line_samples,
                line_samples .* 1.0e-3,
                line_samples .* 1.0e-6,
                line_samples .* 1.0e-4
            ),
            line_distributions,
            LineParameters(
                measurement_parameters.Z.values[1:1, 1:1, :],
                measurement_parameters.Y.values[1:1, 1:1, :],
                frequency
            ),
            4,
            0.95
        )
        for mode in (:hist, :pdf, :ecdf, :qq)
            line_mc_plot = Makie.plot(
                line_mc,
                :R;
                ijk = (1, 1, 2),
                mode,
                data = :both,
                backend = :cairo,
                display_plot = false
            )
            @test line_mc_plot isa UIPlot
            @test line_mc_plot.page.kwargs.selection == (1, 1, 2)
            @test line_mc_plot.page.kwargs.mode === mode
        end

        library = CablesLibrary()
        load!(library; file_name = joinpath(pkgdir(LineCableModels), "test", "cable_test.json"))
        design = first(values(library.data))
        cable_plot = preview(design; backend = :cairo, display_plot = false)
        @test cable_plot isa UIPlot
        @test !isempty(only(cable_plot.page.views).series)
        @test length(cable_plot.page.kwargs.colorbars) == 3
        @test sort!(collect(keys(cable_plot.controls))) ==
              [:export_svg, :legend, :reset]
        cable_legend = cable_plot.controls[:legend]
        cable_panel = only(cable_plot.panels)
        expected_group_order = unique(
            get(series.attributes, :group, Symbol("series_$index"))
        for (index, series) in enumerate(only(cable_plot.page.views).series)
        )
        @test cable_panel.group_order == expected_group_order
        material_entry = first(last(first(cable_legend.entrygroups[])))
        Makie.toggle_visibility!(material_entry)
        @test any(plot_object -> !plot_object.visible[], only(cable_plot.panels).plots)
        Makie.toggle_visibility!(material_entry)
        @test all(plot_object -> plot_object.visible[], only(cable_plot.panels).plots)
        test_golden(cable_plot, "cable_preview"; tolerance = 0.025)

        position = CablePosition(
            design,
            0.0,
            -0.20,
            Dict(component.id => (index == 1 ? 1 : 0)
            for
            (index, component) in enumerate(design.components))
        )
        system = LineCableSystem("reference-system", 1000.0, position)
        earth = EarthModel(frequency, 100.0, 10.0, 1.0)
        system_plot = preview(
            system;
            earth_model = earth,
            backend = :cairo,
            display_plot = false
        )
        @test system_plot isa UIPlot
        @test only(system_plot.page.views).attributes.aspect === :data
        test_golden(system_plot, "system_preview"; tolerance = 0.025)

        zoomed_system_plot = preview(
            system;
            earth_model = earth,
            zoom_factor = 0.5,
            backend = :cairo,
            display_plot = false
        )
        default_limits = only(system_plot.page.views).attributes.limits
        zoomed_limits = only(zoomed_system_plot.page.views).attributes.limits
        @test zoomed_limits[1][2] - zoomed_limits[1][1] <
              default_limits[1][2] - default_limits[1][1]
        @test zoomed_limits[2][2] - zoomed_limits[2][1] <
              default_limits[2][2] - default_limits[2][1]

        material_plot = show_material_scale(backend = :cairo, display_plot = false)
        @test material_plot isa UIPlot
        @test isempty(material_plot.page.views)
        @test length(material_plot.page.kwargs.colorbars) == 3
        @test collect(keys(material_plot.controls)) == [:export_svg]
        test_golden(material_plot, "material_scale")

        primitive_axis = LineCableModels.PlotBuilder.AxisSpec(
            :x,
            LineCableModels.UnitHandler.QuantityTag{:dimensionless}(),
            LineCableModels.UnitHandler.Units(),
            "index",
            :linear
        )
        primitive_view = LineCableModels.PlotBuilder.ViewSpec(
            primitive_axis,
            LineCableModels.PlotBuilder.AxisSpec(
                :y,
                LineCableModels.UnitHandler.QuantityTag{:dimensionless}(),
                LineCableModels.UnitHandler.Units(),
                "index",
                :linear
            ),
            nothing,
            "Heatmap primitive",
            [
                LineCableModels.PlotBuilder.SeriesSpec(
                :heatmap,
                [1.0, 2.0],
                [1.0, 2.0],
                [1.0 2.0; 3.0 4.0],
                nothing
            ),
            ],
            (; kind = :primitive)
        )
        primitive_page = LineCableModels.PlotBuilder.PageSpec(
            "Heatmap primitive",
            (400, 300),
            :single,
            [primitive_view],
            (;
                display_legend = false,
                controls = LineCableModels.PlotBuilder.control_definitions(
                    reset = false,
                    export_svg = false,
                    xlog = false,
                    ylog = false,
                    legend = false
                )
            )
        )
        primitive_render = LineCableModels.PlotBuilder.RenderSpec(
            LineCableModels.DataModel.MaterialScalePlotSpec,
            [primitive_page]
        )
        primitive_plot = only(ui_components.build(
            primitive_render;
            backend = :cairo,
            display = false
        ))
        @test length(only(primitive_plot.panels).plots) == 1
    end
end
