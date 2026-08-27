@testitem "PlotBuilder / Cairo / golden renders and callbacks" tags=[:visual] setup=[
    PlotBuilderTestSupport, UsePlotBuilderSupport, TestNumerics, TestFixtures] begin
    if get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false")!="true"
        error("set LINECABLEMODELS_TEST_PLOTTING=true to run the visual contract")
    else
        using CairoMakie
        using Measurements: measurement

        const TestPlotBuilder=LineCableModels.PlotBuilder
        const TestUIComponents=Base.get_extension(
            LineCableModels,
            :LineCableModelsMakieExt
        ).UIComponents

        struct TestOwnedPlotDefinition<:TestPlotBuilder.AbstractPlotDefinition end
        struct TestOwnedWidgetDefinition<:TestPlotBuilder.AbstractWidgetDefinition end

        const test_legend_placements=Ref(0)
        const test_colorbar_placements=Ref(0)

        function TestUIComponents.draw!(
                context::TestUIComponents.UIContext,
                ::Type{TestOwnedPlotDefinition},
                page::TestPlotBuilder.PlotPage
        )
            payload=page.payload
            axis=Makie.Axis(context.canvas[1, 1])
            curve=Makie.lines!(axis, payload.x, payload.y)
            TestPlotBuilder.register!(
                context,
                axis;
                groups = (local_curve = (curve,),),
                labels = (local_curve = "local curve",),
                data = ((;
                    xdata = payload.x,
                    ydata = payload.y,
                    group = :local_curve,
                    label = "local curve"
                ),)
            )
            return context
        end

        function TestUIComponents.place_legend!(
                context::TestUIComponents.UIContext,
                ::Type{TestOwnedPlotDefinition},
                page::TestPlotBuilder.PlotPage;
                export_mode::Bool
        )
            test_legend_placements[]+=1
            return TestUIComponents.place_legend!(context, page; export_mode)
        end

        function TestUIComponents.place_colorbars!(
                context::TestUIComponents.UIContext,
                ::Type{TestOwnedPlotDefinition},
                page::TestPlotBuilder.PlotPage
        )
            test_colorbar_placements[]+=1
            return TestUIComponents.place_colorbars!(context, page)
        end

        function TestUIComponents.build_widget!(
                context::TestUIComponents.UIContext,
                toolbar,
                column::Int,
                ::TestOwnedWidgetDefinition
        )
            button=TestUIComponents.toolbar_button!(
                context,
                toolbar,
                column;
                key = :test_action,
                icon = TestUIComponents.MI_REFRESH
            )
            TestUIComponents.bind_widget_callback!(
                context,
                button.clicks,
                _->nothing;
                success = "Test action completed"
            )
            return column+1
        end

        include(joinpath(
            pkgdir(LineCableModels),
            "test",
            "support",
            "golden_fixtures.jl"
        ))
        using .GoldenFixtures: custom_layout_plot

        function pixel_error(current, reference)
            size(current)==size(reference)||return Inf
            channel_error=abs.(Makie.red.(current) .- Makie.red.(reference)) .+
                          abs.(Makie.green.(current) .- Makie.green.(reference)) .+
                          abs.(Makie.blue.(current) .- Makie.blue.(reference))
            return sum(channel_error)/(3length(channel_error))
        end

        function test_golden(handle::UIPlot, name; tolerance = 0.015)
            reference_path=joinpath(
                pkgdir(LineCableModels),
                "test",
                "fixtures",
                "golden",
                "$name.png"
            )
            @test isfile(reference_path)
            reference=CairoMakie.FileIO.load(reference_path)
            current=Makie.colorbuffer(handle.figure)
            alternatives=(
                reference,
                reverse(reference; dims = 1),
                reverse(reference; dims = 2),
                reverse(reverse(reference; dims = 1); dims = 2)
            )
            error_value=minimum(pixel_error(current, candidate)
            for candidate in alternatives)
            @testset "golden: $name" begin
                @test error_value<tolerance
            end
            return nothing
        end

        function block_bounds(block)
            layout=block.layoutobservables
            bounding_box=layout.computedbbox[]
            protrusions=layout.protrusions[]
            return (;
                left = bounding_box.origin[1]-protrusions.left,
                right = bounding_box.origin[1]+bounding_box.widths[1]+protrusions.right,
                bottom = bounding_box.origin[2]-protrusions.bottom,
                top = bounding_box.origin[2]+bounding_box.widths[2]+protrusions.top
            )
        end

        visibility_state(handle) = Bool[plot_object.visible[]
                                        for panel in handle.panels
                                        for plot_object in panel.plots]
        legend_labels(legend) = [entry.label[]
                                 for entry in last(first(legend.entrygroups[]))]

        frequency=[50.0, 100.0, 500.0]
        omega=reshape(2π .* frequency, 1, 1, :)
        resistance_values=reshape([1.0, 0.2, 0.2, 2.0], 2, 2, 1) .*
                          ones(1, 1, length(frequency)) .* 1.0e-4
        inductance_values=fill(2.0e-7, 2, 2, length(frequency))
        conductance_values=fill(3.0e-9, 2, 2, length(frequency))
        capacitance_values=fill(4.0e-10, 2, 2, length(frequency))
        series=SeriesImpedance(complex.(resistance_values, inductance_values .* omega))
        shunt=ShuntAdmittance(complex.(conductance_values, capacitance_values .* omega))
        parameters=LineParameters(series, shunt, frequency)

        rlcg=Makie.plot(
            parameters,
            (R, L, G, C);
            backend = :cairo,
            display_plot = false,
            open_export = false
        )
        cartesian=Makie.plot(
            parameters;
            backend = :cairo,
            display_plot = false
        )
        polar=Makie.plot(
            parameters,
            (abs, angle);
            backend = :cairo,
            display_plot = false
        )
        @test rlcg isa Vector{UIPlot}
        @test cartesian isa Vector{UIPlot}
        @test polar isa Vector{UIPlot}
        @test length(rlcg) == 2
        @test length(cartesian) == 2
        @test length(polar) == 2
        @test all(handle -> length(handle.panels) == 2, cartesian)
        @test all(handle -> length(handle.panels) == 2, polar)
        @test all(handle -> length(handle.panels) == 2, rlcg)
        @test rlcg[1].context !== rlcg[2].context
        @test rlcg[1].context.status !== rlcg[2].context.status
        matrix_pairs=((1, 1), (1, 2), (2, 1), (2, 2))
        @test legend_labels(first(rlcg).context.legend) ==
              ["R[$row,$column], L[$row,$column]" for (row, column) in matrix_pairs]
        @test legend_labels(last(rlcg).context.legend) ==
              ["G[$row,$column], C[$row,$column]" for (row, column) in matrix_pairs]
        @test legend_labels(first(cartesian).context.legend) ==
              ["R[$row,$column], X[$row,$column]" for (row, column) in matrix_pairs]
        @test legend_labels(first(polar).context.legend) ==
              ["|Z|[$row,$column], ∠Z[$row,$column]" for (row, column) in matrix_pairs]
        first_handle=first(rlcg)
        @test fieldnames(typeof(first_handle)) == (:render, :page, :context)
        @test fieldnames(typeof(first_handle.context)) == (
            :backend,
            :interactive,
            :window,
            :figure,
            :shell,
            :canvas,
            :status,
            :panels,
            :widgets,
            :legend,
            :colorbars,
            :responsive_legend,
            :legend_slot_grid,
            :observers,
            :export_state,
            :plot_reference
        )
        @test first_handle.figure === first_handle.context.figure
        @test first_handle.panels === first_handle.context.panels
        @test first_handle.controls === first_handle.context.widgets
        @test first_handle.context.export_state.theme ===
              first_handle.page.export_definition.theme
        @test all(
            name -> !hasproperty(first_handle.page.payload, name),
            (:title, :key, :legend, :colorbars, :widgets, :export_definition)
        )
        series_plots=Makie.plot(
            series,
            frequency;
            backend = :cairo,
            display_plot = false
        )
        shunt_plots=Makie.plot(
            shunt,
            frequency,
            (G, C);
            backend = :cairo,
            display_plot = false
        )
        @test series_plots isa Vector{UIPlot}
        @test shunt_plots isa Vector{UIPlot}
        @test length(series_plots) == 1
        @test length(shunt_plots) == 1
        @test length(only(series_plots).panels) == 2
        @test length(only(shunt_plots).panels) == 2

        comparison_parameters=LineParameters(
            1.01 .* Z(parameters),
            1.01 .* Y(parameters),
            frequency
        )
        comparison_plots=Makie.plot(
            parameters,
            comparison_parameters;
            legend = ("reference", "candidate"),
            backend = :cairo,
            display_plot = false
        )
        @test comparison_plots isa Vector{UIPlot}
        @test length(comparison_plots) == 4
        @test all(plot -> length(plot.panels) == 4, comparison_plots)
        @test length(Makie.plot(
            parameters,
            comparison_parameters;
            legend = ("reference", "candidate"),
            requests = (Z,),
            backend = :cairo,
            display_plot = false
        )) == 2
        @test length(Makie.plot(
            parameters,
            comparison_parameters;
            legend = ("reference", "candidate"),
            requests = (Y,),
            backend = :cairo,
            display_plot = false
        )) == 2
        test_golden(first(rlcg), "line_rlcg")
        test_golden(first(cartesian), "line_zy_cartesian")
        test_golden(first(polar), "line_zy_polar")

        measurement_parameters=LineParameters(
            complex.(measurement.(resistance_values, resistance_values .* 0.05),
                measurement.(inductance_values .* omega, inductance_values .* omega .*
                                                         0.05)),
            complex.(measurement.(conductance_values, conductance_values .* 0.05),
                measurement.(capacitance_values .* omega, capacitance_values .* omega .*
                                                          0.05)),
            frequency
        )
        measurement_plots=Makie.plot(
            measurement_parameters,
            (R, L, G, C);
            backend = :cairo,
            display_plot = false
        )
        measurement_plot=first(measurement_plots)
        measurement_payload=measurement_plot.page.payload
        @test length(first(measurement_plot.panels).plots) >
              prod(length, first(measurement_payload.coordinates)[1:2])
        test_golden(measurement_plot, "line_measurements")

        measurement_conductance=measurement_plots[2]
        measurement_conductance_axis=first(measurement_conductance.panels).axis
        measurement_limits=measurement_conductance_axis.finallimits[]
        measurement_ymin=measurement_limits.origin[2]
        measurement_ymax=measurement_ymin+measurement_limits.widths[2]
        @test measurement_ymin < 2.85e-6
        @test measurement_ymax > 3.15e-6
        measurement_conductance.controls[:ylog].active[]=true
        log_measurement_limits=measurement_conductance_axis.finallimits[]
        log_measurement_ymin=log_measurement_limits.origin[2]
        log_measurement_ymax=log_measurement_ymin+log_measurement_limits.widths[2]
        @test log_measurement_ymin < 2.85e-6
        @test log_measurement_ymax > 3.15e-6

        handle=first(rlcg)
        @test handle.figure.scene.backgroundcolor[] == Makie.to_color(:grey90)
        padding=handle.figure.layout.alignmode[].padding
        @test (padding.left, padding.right, padding.bottom, padding.top) ==
              (20.0, 20.0, 28.0, 28.0)
        @test handle.figure.layout.colsizes[2] == Makie.Auto(true)
        @test handle.figure.layout.addedrowgaps == fill(Makie.Fixed(6), 2)
        @test handle.figure.layout.addedcolgaps == [Makie.Fixed(12)]
        @test handle.controls[:reset].buttoncolor[] == Makie.RGBf(0.94, 0.94, 0.94)
        @test occursin("\\ue5d5", sprint(show, handle.controls[:reset].label[]))
        @test occursin("\\ue161", sprint(show, handle.controls[:export_svg].label[]))
        line_axis=first(handle.panels).axis
        initial_line_limits=line_axis.finallimits[]
        initial_xlabel=line_axis.xlabel[]
        handle.controls[:xlog].active[]=true
        @test line_axis.xscale[] === Makie.log10
        @test line_axis.xlabel[] == "Frequency [Hz]"
        log_line_limits=line_axis.finallimits[]
        @test log_line_limits.origin[1] < minimum(frequency)
        @test log_line_limits.origin[1] + log_line_limits.widths[1] > maximum(frequency)
        handle.controls[:xlog].active[]=false
        @test line_axis.xscale[] === Makie.identity
        @test line_axis.xlabel[] == initial_xlabel
        @test line_axis.finallimits[] == initial_line_limits
        handle.controls[:xlog].active[]=true
        handle.controls[:ylog].active[]=true
        @test line_axis.yscale[] === Makie.log10
        handle.controls[:reset].clicks[]+=1
        @test handle.context.status[] == "Axis limits reset"
        legend=handle.controls[:legend]
        first_entry=first(last(first(legend.entrygroups[])))
        Makie.toggle_visibility!(first_entry)
        @test any(
            plot_object -> !plot_object.visible[],
            Iterators.flatten(panel.plots for panel in handle.panels)
        )
        Makie.xlims!(first(handle.panels).axis, 100.0, 300.0)
        ui_components=Base.get_extension(
            LineCableModels,
            :LineCableModelsMakieExt
        ).UIComponents
        default_export_theme=ui_components._theme(
            export_mode = true,
            export_theme = :default
        )
        publication_export_theme=ui_components._theme(
            export_mode = true,
            export_theme = :publication
        )
        @test !haskey(default_export_theme[:fonts], :regular)
        @test haskey(publication_export_theme[:fonts], :regular)
        @test haskey(publication_export_theme[:fonts], :italic)
        @test publication_export_theme[:Axis][:titlesize][] == 15
        @test publication_export_theme[:Axis][:xticklabelsize][] == 14
        current_recipe=ui_components._current_recipe(handle)
        current_panel=first(only(current_recipe.pages).payload.runtime.panels)
        @test current_panel.xscale === :log10
        @test current_panel.yscale === :log10
        @test !isempty(current_panel.hidden_groups)
        @test current_panel.current_limits !== nothing
        @test collect(current_panel.current_limits[1]) ≈ [100.0, 300.0]
        Makie.toggle_visibility!(first_entry)
        @test all(
            plot_object -> plot_object.visible[],
            Iterators.flatten(panel.plots for panel in handle.panels)
        )

        contrast_conductance=repeat(
            [1.0e-3 1.0e-9; 1.0e-9 2.0e-9],
            1,
            1,
            length(frequency)
        )
        contrast_parameters=LineParameters(
            copy(observe(parameters, Z)),
            complex.(contrast_conductance, capacitance_values .* omega),
            frequency
        )
        contrast_plot=only(Makie.plot(
            contrast_parameters,
            (G,);
            backend = :cairo,
            display_plot = false
        ))
        contrast_axis=only(contrast_plot.panels).axis
        initial_contrast_limits=contrast_axis.finallimits[]
        initial_contrast_max=initial_contrast_limits.origin[2]+
        initial_contrast_limits.widths[2]
        @test initial_contrast_limits.origin[2] < 0
        mktempdir() do directory
            untouched_svg=export_svg(
                contrast_plot;
                path = joinpath(directory, "untouched-linear.svg"),
                open_file = false
            )
            @test filesize(untouched_svg) > 100
            @test occursin("<svg", read(untouched_svg, String))
        end
        contrast_legend=contrast_plot.controls[:legend]
        dominant_entry=first(last(first(contrast_legend.entrygroups[])))
        Makie.toggle_visibility!(dominant_entry)
        Makie.update_state_before_display!(contrast_plot.figure)
        hidden_contrast_limits=contrast_axis.finallimits[]
        hidden_contrast_max=hidden_contrast_limits.origin[2]+
        hidden_contrast_limits.widths[2]
        @test hidden_contrast_max < initial_contrast_max * 1.0e-3
        @test contrast_plot.context.status[] == "Axis limits fitted to visible series"
        current_contrast_recipe=ui_components._current_recipe(contrast_plot)
        current_contrast_panel=only(
            only(current_contrast_recipe.pages).payload.runtime.panels,
        )
        @test first(only(contrast_plot.panels).group_order) in
              current_contrast_panel.hidden_groups
        @test collect(current_contrast_panel.current_limits[2]) ≈ [
            hidden_contrast_limits.origin[2],
            hidden_contrast_max
        ]
        Makie.toggle_visibility!(dominant_entry)
        Makie.update_state_before_display!(contrast_plot.figure)
        restored_contrast_limits=contrast_axis.finallimits[]
        restored_contrast_max=restored_contrast_limits.origin[2]+
        restored_contrast_limits.widths[2]
        @test restored_contrast_max ≈ initial_contrast_max

        susceptance_handle=last(cartesian)
        susceptance_axis=last(susceptance_handle.panels).axis
        @test susceptance_axis.ylabel[] isa Makie.RichText
        @test occursin("−3", sprint(show, susceptance_axis.ylabel[]))
        initial_susceptance_limits=susceptance_axis.finallimits[]
        susceptance_handle.controls[:ylog].active[]=true
        @test susceptance_axis.yscale[] === Makie.log10
        @test susceptance_axis.ytickformat[] === Makie.automatic
        susceptance_quantity=LineCableModels.Units.quantity(B)
        susceptance_unit=LineCableModels.Units.display_unit(
            susceptance_quantity,
            basis(parameters)
        )
        @test susceptance_axis.ylabel[] == LineCableModels.Units.label(
            susceptance_quantity,
            susceptance_unit
        )
        limits=susceptance_axis.finallimits[]
        ymin=limits.origin[2]
        ymax=ymin+limits.widths[2]
        tick_values, tick_labels=Makie.get_ticks(
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
        susceptance_handle.controls[:ylog].active[]=false
        @test susceptance_axis.yscale[] === Makie.identity
        @test susceptance_axis.ylabel[] isa Makie.RichText
        @test susceptance_axis.finallimits[] == initial_susceptance_limits
        @test susceptance_axis.ytickformat[]([1.23456e-3, 9.87654e-3]) ==
              ["1.235", "9.877"]

        threshold_parameters=LineParameters(
            copy(Z(parameters)),
            complex.(conductance_values, 10 .* capacitance_values .* omega),
            frequency
        )
        threshold_plot=only(Makie.plot(
            threshold_parameters,
            (B,);
            backend = :cairo,
            display_plot = false
        ))
        threshold_axis=only(threshold_plot.panels).axis
        @test threshold_axis.ylabel[] isa Makie.RichText
        @test only(threshold_plot.panels).metadata.yaxis.exponent == -2
        @test threshold_axis.ytickformat[]([1.23456e-2, 9.87654e-2]) ==
              ["1.235", "9.877"]

        conductance_handle=rlcg[2]
        conductance_axis=first(conductance_handle.panels).axis
        conductance_handle.controls[:ylog].active[]=true
        conductance_limits=conductance_axis.finallimits[]
        @test conductance_limits.origin[2] == 1.0e-6
        @test conductance_limits.origin[2] + conductance_limits.widths[2] == 1.0e-5
        conductance_ticks, conductance_labels=Makie.get_ticks(
            conductance_axis.yticks[],
            conductance_axis.yscale[],
            conductance_axis.ytickformat[],
            1.0e-6,
            1.0e-5
        )
        @test conductance_ticks == [1.0e-6, 1.0e-5]
        @test all(label -> label isa Makie.RichText, conductance_labels)

        conductance_handle.controls[:ylog].active[]=false
        capacitance_axis=last(rlcg[2].panels).axis
        capacitance_limits=capacitance_axis.finallimits[]
        @test capacitance_limits.origin[2] == 0.38
        @test capacitance_limits.origin[2] + capacitance_limits.widths[2] ≈ 0.42

        zero_conductance_parameters=LineParameters(
            copy(Z(parameters)),
            complex.(zeros(size(conductance_values)), capacitance_values .* omega),
            frequency
        )
        zero_shunt_plot=only(Makie.plot(
            zero_conductance_parameters,
            (G, C);
            xscale = :log10,
            backend = :cairo,
            display_plot = false,
            controls = false
        ))
        zero_conductance_axis=first(zero_shunt_plot.panels).axis
        zero_conductance_limits=zero_conductance_axis.finallimits[]
        @test zero_conductance_limits.origin[2] == -sqrt(eps(Float64))
        @test zero_conductance_limits.origin[2] +
              zero_conductance_limits.widths[2] == sqrt(eps(Float64))
        @test zero_conductance_axis.ylabel[] isa Makie.RichText
        @test occursin("−8", sprint(show, zero_conductance_axis.ylabel[]))
        @test zero_conductance_axis.ytickformat[]([-1.0e-8, 0.0, 1.0e-8]) ==
              ["-1", "0", "1"]

        mktempdir() do directory
            svg_path=joinpath(directory, "line parameters.svg")
            exported=export_svg(handle; path = svg_path, open_file = false)
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
                first_default=export_svg(handle; open_file = false)
                second_default=export_svg(handle; open_file = false)
                @test first_default != second_default
                @test dirname(first_default) == directory
                @test occursin(
                    r"^series_impedance_\d{8}_\d{6}(?:_\d+)?\.svg$",
                    basename(first_default)
                )
                @test occursin("rgb(100%, 100%, 100%)", read(first_default, String))
                @test filesize(second_default) > 100
            end

            publication_path=joinpath(directory, "publication.svg")
            publication=export_svg(
                handle;
                path = publication_path,
                theme = :publication,
                open_file = false
            )
            @test publication == publication_path
            @test filesize(publication) > 100
            @test read(publication) != read(exported)
        end

        fallback_export=cd(pkgdir(LineCableModels)) do
            export_svg(handle; open_file = false)
        end
        @test dirname(fallback_export) ==
              joinpath(tempdir(), "linecablemodels-exports")
        @test filesize(fallback_export) > 100
        rm(fallback_export)

        function combined_monte_carlo_histogram(result)
            target=only(LineCableModels.Grammar.unit_targets(
                (R,),
                basis(result);
                length_prefix = :kilo
            ))
            published=observables(
                result,
                (
                    (samples, R, 1, Colon()),
                    (histograms, R, 1)
                );
                units = (target, target)
            )
            sample, model=published
            density=(;
                values = model.values.density,
                quantity = LineCableModels.Units.Quantity{:probability_density}(),
                unit = inv(target)
            )
            return TestPlotBuilder.plotwindow(
                title = "R histogram",
                backend = :cairo,
                display_plot = false
            ) do context
                axis=TestPlotBuilder.axis!(
                    context,
                    context.canvas[1, 1],
                    sample,
                    density;
                    title = "R histogram",
                    ylabel = "pdf"
                )
                bars=Makie.hist!(
                    axis,
                    sample.values;
                    bins = model.values.edges,
                    normalization = :pdf,
                    label = "samples"
                )
                curve=Makie.stairs!(
                    axis,
                    model.values.edges,
                    [model.values.density; last(model.values.density)];
                    step = :post,
                    color = :red,
                    linewidth = 2,
                    label = "model PDF"
                )
                TestPlotBuilder.register!(
                    context,
                    axis;
                    groups = (samples = (bars,), model = (curve,)),
                    labels = (samples = "samples", model = "model PDF"),
                    data = (
                        (;
                            xdata = sample.values,
                            ydata = nothing,
                            group = :samples,
                            label = "samples"
                        ),
                        (;
                            xdata = model.values.edges,
                            ydata = [model.values.density; last(model.values.density)],
                            group = :model,
                            label = "model PDF"
                        )
                    )
                )
                center=only(unique(model.values.density))
                Makie.ylims!(axis, 0.95center, 1.05center)
            end
        end

        mc_result=TestFixtures.cable_monte_carlo_result()
        histogram=only(histograms(mc_result)).R
        monte_carlo_plots=(
            hist = combined_monte_carlo_histogram(mc_result),
            pdf = Makie.stairs(mc_result, R;
                backend = :cairo, display_plot = false),
            ecdf = Makie.ecdfplot(mc_result, R;
                backend = :cairo, display_plot = false),
            qq = Makie.qqplot(mc_result, R;
                backend = :cairo, display_plot = false)
        )
        for (name, mc_plot) in pairs(monte_carlo_plots)
            @test mc_plot isa UIPlot
            @test !haskey(mc_plot.controls, :xlog)
            @test !haskey(mc_plot.controls, :ylog)
            test_golden(mc_plot, "mc_$name")
        end

        line_samples=reshape(collect(1.0:12.0), 1, 1, 3, 4)
        summarize(values) = map(
            index->SampleSummary(collect(view(values, index.I..., :))),
            CartesianIndices(size(values)[1:3])
        )
        line_statistics=(
            R = summarize(line_samples),
            L = summarize(line_samples .* 1.0e-3),
            C = summarize(line_samples .* 1.0e-6),
            G = summarize(line_samples .* 1.0e-4)
        )
        line_histograms=(
            R = fill(histogram, 1, 1, 3),
            L = fill(histogram, 1, 1, 3),
            C = fill(histogram, 1, 1, 3),
            G = fill(histogram, 1, 1, 3)
        )
        line_representation=LineParameters(
            Z(parameters)[1:1, 1:1, :],
            Y(parameters)[1:1, 1:1, :],
            frequency
        )
        line_sample_values=(
            R = line_samples,
            L = line_samples .* 1.0e-3,
            C = line_samples .* 1.0e-6,
            G = line_samples .* 1.0e-4
        )
        line_formulation=MonteCarlo(
            Formulation(); trials = 4, seed = 2,
            return_samples = true, return_histograms = true
        )
        line_mc=MonteCarloResult(
            line_formulation,
            [line_representation],
            [line_statistics],
            [line_sample_values],
            [line_histograms],
            UInt64(2),
            UInt64[2],
            [4]
        )
        request=@observe R[1, 1, 2]
        line_mc_plots=(
            Makie.hist(line_mc, request;
                backend = :cairo, display_plot = false),
            Makie.stairs(line_mc, request;
                backend = :cairo, display_plot = false),
            Makie.ecdfplot(line_mc, request;
                backend = :cairo, display_plot = false),
            Makie.qqplot(line_mc, request;
                backend = :cairo, display_plot = false)
        )
        for line_mc_plot in line_mc_plots
            @test line_mc_plot isa UIPlot
            @test occursin("R[1,1,2]", line_mc_plot.page.title)
        end

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
        cable_plot=preview(design; backend = :cairo, display_plot = false)
        @test cable_plot isa UIPlot
        @test !isempty(cable_plot.page.payload.polygons)
        @test length(cable_plot.page.colorbars) == 3
        @test sort!(collect(keys(cable_plot.controls))) ==
              [:export_svg, :legend, :reset]
        cable_legend=cable_plot.controls[:legend]
        cable_entries=last(first(cable_legend.entrygroups[]))
        cable_elements=Iterators.flatten(entry.elements for entry in cable_entries)
        for element in cable_elements
            element isa Makie.PolyElement||continue
            @test Makie.to_color(element.attributes[:polystrokecolor][]) ==
                  Makie.to_color(:transparent)
            @test iszero(element.attributes[:polystrokewidth][])
        end
        cable_panel=only(cable_plot.panels)
        expected_group_order=unique(
            polygon.group for polygon in cable_plot.page.payload.polygons
        )
        @test cable_panel.group_order == expected_group_order
        material_entry=first(last(first(cable_legend.entrygroups[])))
        Makie.toggle_visibility!(material_entry)
        @test any(plot_object -> !plot_object.visible[], only(cable_plot.panels).plots)
        Makie.toggle_visibility!(material_entry)
        @test all(plot_object -> plot_object.visible[], only(cable_plot.panels).plots)
        test_golden(cable_plot, "cable_preview"; tolerance = 0.025)

        collection_designs=[
            design,
            equivalent(design; new_id = "equivalent"),
            CableDesign(
                "second detailed design",
                design.components;
                nominal_data = design.nominal_data
            )
        ]
        collection_plot=preview(
            collection_designs;
            backend = :cairo,
            display_plot = false,
            open_export = false
        )
        @test collection_plot isa UIPlot
        @test length(collection_plot.panels) == 3
        @test getproperty.(collection_plot.page.payload.panels, :position) ==
              [(1, 1), (1, 2), (2, 1)]
        @test [panel.axis.title[] for panel in collection_plot.panels] ==
              getproperty.(collection_designs, :cable_id)
        @test sort!(collect(keys(collection_plot.controls))) == [:export_svg, :reset]
        @test collection_plot.context.legend === nothing
        @test length(collection_plot.context.colorbars) == 3
        side_valign=collection_plot.context.shell.side.valign[]
        @test side_valign == convert(typeof(side_valign), :top)
        collection_replay=ui_components._current_recipe(collection_plot)
        @test all(
            panel -> panel.payload.runtime !== nothing,
            only(collection_replay.pages).payload.panels
        )
        replayed_collection=only(ui_components.build(
            collection_replay;
            backend = :cairo,
            display = false,
            controls = false,
            export_mode = true
        ))
        @test length(replayed_collection.panels) == 3

        compact_cable_plot=preview(
            design;
            size = (900, 350),
            backend = :cairo,
            display_plot = false,
            open_export = false
        )
        compact_legend=compact_cable_plot.controls[:legend]
        compact_panel=only(compact_cable_plot.panels)
        complete_labels=[compact_panel.group_labels[group]
                         for group in compact_panel.group_order
                         if haskey(compact_panel.group_labels, group)]
        compact_entries=last(first(compact_legend.entrygroups[]))
        compact_labels=legend_labels(compact_legend)
        @test last(compact_labels) == "(...)"
        @test compact_labels[1:(end - 1)] ==
              complete_labels[1:(length(compact_labels) - 1)]
        @test length(compact_labels) < length(complete_labels)
        @test length(filter(block -> block isa Makie.Colorbar,
            compact_cable_plot.figure.content)) == 3

        first_compact_entry=first(compact_entries)
        Makie.toggle_visibility!(first_compact_entry)
        hidden_state=visibility_state(compact_cable_plot)
        @test any(!, hidden_state)
        Makie.resize!(compact_cable_plot.figure, 900, 700)
        Makie.update_state_before_display!(compact_cable_plot.figure)
        @test legend_labels(compact_legend) == complete_labels
        @test visibility_state(compact_cable_plot) == hidden_state
        Makie.resize!(compact_cable_plot.figure, 900, 350)
        Makie.update_state_before_display!(compact_cable_plot.figure)
        compact_entries=last(first(compact_legend.entrygroups[]))
        @test last(compact_entries).label[] == "(...)"
        @test visibility_state(compact_cable_plot) == hidden_state
        Makie.toggle_visibility!(last(compact_entries))
        @test visibility_state(compact_cable_plot) == hidden_state

        legend_bounds=block_bounds(compact_legend)
        canvas_bounds=block_bounds(compact_panel.axis)
        scale_blocks=filter(
            block->block isa Makie.Colorbar,
            compact_cable_plot.figure.content
        )
        scale_bounds=block_bounds.(scale_blocks)
        material_labels=filter(compact_cable_plot.figure.content) do block
            block isa Makie.Label||return false
            block_bounds(block).left>=legend_bounds.left-1
        end
        material_label_bounds=block_bounds.(material_labels)
        @test length(material_labels) == 3
        @test canvas_bounds.right < legend_bounds.left
        @test all(bounds -> bounds.right < scale_bounds[1].left,
            material_label_bounds)
        @test all(bounds -> bounds.left >= canvas_bounds.right, material_label_bounds)
        @test all(bounds -> bounds.right <= legend_bounds.right, scale_bounds)
        @test legend_bounds.bottom > maximum(getfield.(scale_bounds, :top))
        @test all(
            block -> iszero(block.layoutobservables.protrusions[].left) &&
                     iszero(block.layoutobservables.protrusions[].right),
            scale_blocks)

        current_compact_recipe=ui_components._current_recipe(compact_cable_plot)
        @test !isempty(only(current_compact_recipe.pages).payload.runtime.hidden_groups)
        exported_compact=only(ui_components.build(
            current_compact_recipe;
            backend = :cairo,
            display = false,
            controls = false,
            export_mode = true
        ))
        exported_legend=only(filter(
            block->block isa Makie.Legend,
            exported_compact.figure.content
        ))
        @test legend_labels(exported_legend) == complete_labels
        fitted_export_size=ui_components._fit_export_content!(
            exported_compact.figure,
            exported_compact.page
        )
        @test fitted_export_size[2] > compact_cable_plot.page.size[2]
        exported_legend_bounds=block_bounds(exported_legend)
        exported_scale_bounds=block_bounds.(filter(
            block->block isa Makie.Colorbar,
            exported_compact.figure.content
        ))
        @test exported_legend_bounds.bottom >
              maximum(getfield.(exported_scale_bounds, :top))
        @test any(!plot_object.visible[] for panel in exported_compact.panels
        for plot_object in panel.plots)

        mktempdir() do directory
            responsive_svg=export_svg(
                compact_cable_plot;
                path = joinpath(directory, "responsive-preview.svg"),
                open_file = false
            )
            svg=read(responsive_svg, String)
            height_match=match(r"<svg[^>]*height=\"([0-9]+)\"", svg)
            @test height_match !== nothing
            @test Base.parse(Int, only(height_match.captures)) >
                  compact_cable_plot.page.size[2]
            @test filesize(responsive_svg) > 100
        end
        @test last(last(first(compact_legend.entrygroups[]))).label[] == "(...)"
        @test visibility_state(compact_cable_plot) == hidden_state
        Makie.toggle_visibility!(first(last(first(compact_legend.entrygroups[]))))
        @test all(visibility_state(compact_cable_plot))
        test_golden(compact_cable_plot, "cable_preview_compact"; tolerance = 0.025)

        position=CablePosition(
            design,
            0.0,
            -0.20,
            Dict(component.id=>(index==1 ? 1 : 0)
            for
            (index, component) in enumerate(design.components))
        )
        system=LineCableSystem("reference-system", 1000.0, position)
        earth=EarthModel(100.0, 10.0, 1.0)
        system_plot=preview(
            system;
            earth_model = earth,
            backend = :cairo,
            display_plot = false
        )
        @test system_plot isa UIPlot
        @test only(system_plot.panels).metadata.aspect === :data
        test_golden(system_plot, "system_preview"; tolerance = 0.025)

        zoomed_system_plot=preview(
            system;
            earth_model = earth,
            zoom_factor = 0.5,
            backend = :cairo,
            display_plot = false
        )
        default_limits=system_plot.page.payload.limits
        zoomed_limits=zoomed_system_plot.page.payload.limits
        @test zoomed_limits[1][2] - zoomed_limits[1][1] <
              default_limits[1][2] - default_limits[1][1]
        @test zoomed_limits[2][2] - zoomed_limits[2][1] <
              default_limits[2][2] - default_limits[2][1]

        material_plot=LineCableModels.DataModel.show_material_scale(
            backend = :cairo,
            display_plot = false
        )
        @test material_plot isa UIPlot
        @test length(material_plot.page.colorbars) == 3
        @test collect(keys(material_plot.controls)) == [:export_svg]
        test_golden(material_plot, "material_scale")

        local_payload=(;
            x = [1.0, 2.0, 3.0],
            y = [2.0, 3.0, 5.0]
        )
        local_recipe=TestPlotBuilder.PlotRecipe(
            TestOwnedPlotDefinition,
            [TestPlotBuilder.PlotPage(
                "Test-owned definition",
                (400, 300),
                (; kind = :test_owned),
                local_payload;
                legend = TestPlotBuilder.LegendDefinition(),
                widgets = (TestOwnedWidgetDefinition(),),
                export_definition = TestPlotBuilder.ExportDefinition(
                    name = "test-owned-definition",
                    open_file = false
                )
            )]
        )
        local_plot=only(TestUIComponents.build(
            local_recipe;
            backend = :cairo,
            display = false,
            controls = true
        ))
        @test local_plot.page.payload === local_payload
        @test length(only(local_plot.panels).plots) == 1
        @test legend_labels(local_plot.context.legend) == ["local curve"]
        @test test_legend_placements[] == 1
        @test test_colorbar_placements[] == 1
        @test haskey(local_plot.controls, :test_action)
        local_plot.controls[:test_action].clicks[]=1
        @test local_plot.context.status[] == "Test action completed"

        primitive_plot=LineCableModels.PlotBuilder.plotwindow(
            title = "Heatmap primitive",
            size = (400, 300),
            backend = :cairo,
            display_plot = false,
            controls = false,
            legend = false
        ) do ui
            axis=Makie.Axis(ui.canvas[1, 1]; xgridvisible = false)
            image=Makie.heatmap!(
                axis,
                [1.0, 2.0],
                [1.0, 2.0],
                [1.0 2.0; 3.0 4.0]
            )
            LineCableModels.PlotBuilder.register!(
                ui,
                axis;
                groups = (heatmap = image,),
                data = ((;
                    xdata = [1.0, 2.0],
                    ydata = [1.0, 2.0],
                    group = :heatmap,
                    label = nothing
                ),)
            )
        end
        @test length(only(primitive_plot.panels).plots) == 1
        @test !only(primitive_plot.panels).axis.xgridvisible[]

        custom_plot=custom_layout_plot()
        @test sprint(show, MIME"text/plain"(), custom_plot) ==
              "UIPlot(title=\"Nested PlotBuilder layout\", panels=3, backend=:cairo)"
        @test custom_plot.page.key == (; kind = :native_canvas)
        @test hasproperty(custom_plot.page.payload, :callback)
        @test length(custom_plot.panels) == 3
        @test all(panel -> panel.axis.xlabel[] == "x", custom_plot.panels)
        @test all(panel -> panel.axis.ylabel[] == "y", custom_plot.panels)
        top_width=custom_plot.panels[1].axis.scene.viewport[].widths[1]
        bottom_widths=[panel.axis.scene.viewport[].widths[1]
                       for panel in custom_plot.panels[2:3]]
        @test top_width > maximum(bottom_widths)
        legend=custom_plot.controls[:legend]
        colorbar=only(filter(block->block isa Makie.Colorbar, custom_plot.figure.content))
        legend_box=legend.layoutobservables.computedbbox[]
        colorbar_box=colorbar.layoutobservables.computedbbox[]
        colorbar_label=only(filter(custom_plot.figure.content) do block
            block isa Makie.Label||return false
            block_bounds(block).left>=legend_box.origin[1]-1
        end)
        colorbar_label_box=colorbar_label.layoutobservables.computedbbox[]
        @test legend.halign[] === :left
        @test colorbar_label_box.origin[1] ≈ legend_box.origin[1] atol = 1
        @test colorbar_label_box.origin[1]+colorbar_label_box.widths[1] <
              colorbar_box.origin[1]
        @test legend_box.widths[1] ≈ ui_components.LEGEND_DOCK_WIDTH atol = 1
        @test colorbar_box.widths[1] ≈ ui_components.COLORBAR_WIDTH atol = 1

        frequency_observation=observables(
            parameters,
            (frequencies,)
        )|>only
        native_observation=(
            values = [1.0, 2.0],
            quantity = frequency_observation.quantity,
            unit = frequency_observation.unit
        )
        scientific_window=LineCableModels.PlotBuilder.plotwindow(;
            title = "Scientific native axis",
            size = (400, 300),
            backend = :cairo,
            display_plot = false,
            controls = false,
            legend = false,
            open_export = false
        ) do ui
            axis=LineCableModels.PlotBuilder.axis!(
                ui,
                ui.canvas[1, 1],
                native_observation,
                native_observation
            )
            lines!(axis, native_observation.values, native_observation.values)
        end
        scientific_axis=only(scientific_window.panels).axis
        expected_scientific_label=LineCableModels.Units.label(
            native_observation.quantity,
            native_observation.unit
        )
        @test scientific_axis.xlabel[] == expected_scientific_label
        @test scientific_axis.ylabel[] == expected_scientific_label

        small_observation=(;
            values = [1.0e-8, 2.0e-8],
            quantity = native_observation.quantity,
            unit = native_observation.unit
        )
        exponent_window=LineCableModels.PlotBuilder.plotwindow(;
            title = "Scientific exponent",
            size = (400, 300),
            backend = :cairo,
            display_plot = false,
            controls = false,
            legend = false,
            open_export = false
        ) do ui
            axis=LineCableModels.PlotBuilder.axis!(
                ui,
                ui.canvas[1, 1],
                native_observation,
                small_observation
            )
            lines!(axis, native_observation.values, small_observation.values)
        end
        exponent_panel=only(exponent_window.panels)
        @test exponent_panel.metadata.yaxis.exponent == -8
        @test exponent_panel.axis.ylabel[] isa Makie.RichText
        @test exponent_panel.axis.ytickformat[]([1.0e-8, 2.0e-8]) == ["1", "2"]

        test_golden(custom_plot, "custom_layout"; tolerance = 0.02)
        custom_export_render=custom_layout_plot(
            controls = false,
            export_mode = true
        )
        @test custom_export_render.figure.layout.rowsizes[1] == Makie.Fixed(0)
        @test length(custom_export_render.figure.layout.rowsizes) == 2 ||
              custom_export_render.figure.layout.rowsizes[3] == Makie.Fixed(0)
        mktempdir() do directory
            custom_svg=export_svg(
                custom_plot;
                path = joinpath(directory, "custom-layout.svg"),
                open_file = false
            )
            @test filesize(custom_svg) > 100
            @test occursin("<svg", read(custom_svg, String))
        end
    end
end
