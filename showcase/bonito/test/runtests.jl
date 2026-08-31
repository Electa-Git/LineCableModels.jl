using Test
using LineCableModels
using Bonito
using WGLMakie
using CairoMakie
using GraphMakie
using Graphs
using Markdown
using NetworkLayout
import Measurements
import PowerImpedance

CairoMakie.activate!()
include(joinpath(@__DIR__, "..", "app.jl"))

const LiveManual = find_presentation("live-manual", PRESENTATION_DESCRIPTORS)
const ICHQP2026 = find_presentation("ICHQP2026", PRESENTATION_DESCRIPTORS)
const FundamentalsPage = getfield(
    ICHQP2026.namespace,
    :ICHQP2026FundamentalsDeck
)
const CableConstantsPage = getfield(
    ICHQP2026.namespace,
    :ICHQP2026CableConstantsDeck
)
const LineParametersPage = getfield(
    ICHQP2026.namespace,
    :ICHQP2026LineParametersDeck
)
const ApplicationCasePage = getfield(
    ICHQP2026.namespace,
    :ICHQP2026ApplicationCaseDeck
)
const CablePage = getfield(LiveManual.namespace, :CableDesignPage)
const TransitionPage = getfield(LiveManual.namespace, :OHLUGCTransitionPage)
const ReplPage = getfield(LiveManual.namespace, :LiveJuliaPage)
const PowerImpedanceDiagramExt = getfield(
    LiveManual.namespace,
    :PowerImpedanceDiagramExt
)
const TemplatePresentation = load_presentation_descriptor(
    joinpath(@__DIR__, "..", "presentations", "_template");
    slug = "starter-template"
)

@testset "Bonito cable showcase" begin
    @testset "page authoring" begin
        layout = webgrid(
            [:controls :figure; :controls :prose];
            columns = "2fr 5fr",
            rows = "1fr 1fr",
            controls = webpart("controls"; kind = :controls),
            figure = webpart("figure"; kind = :plot),
            prose = PageAuthoring.prose(md"Inline mathematics: ``x^2``.")
        )
        page = deck_page(
            "First view";
            id = "first-view",
            build = (session, state) -> (; body = state)
        )
        descriptor = deck_descriptor(
            id = "fixture-deck",
            group = "Fixtures",
            title = "Fixture deck",
            order = 1,
            setup = identity,
            pages = (page,),
            fixture = true
        )
        resource_state = preflight_state()
        resource = preflight_resource(
            id = "fixture-resource",
            title = "Fixture resource",
            state = resource_state,
            activate = () -> set_preflight!(
                resource_state,
                :hot,
                1.0,
                "Fixture resource hot"
            )
        )
        preparation = preparation_status()
        masters = (
            layout_single("main"),
            layout_two_columns("left", "right"),
            layout_three_columns("left", "center", "right"),
            layout_top_split("top", "left", "right"),
            layout_split_bottom("left", "right", "bottom"),
            layout_sidebar_main("sidebar", "main"),
            layout_controls_plot("controls", "plot"; footer = "footer"),
            layout_controls_plot(
                "controls",
                "plot";
                footer = "footer",
                footer_placement = :under_plot
            ),
            action_row("Reset plot view", "Reset selectors"),
            layout_two_rows("top", "bottom"; ratio = :short_tall),
            layout_quad("tl", "tr", "bl", "br")
        )
        selector = Observable(9)
        reset_selectors!(selector => 3)
        silent_figure = Figure(size = (64, 64))
        silent_wgl = wgl_figure(silent_figure)
        silent_spinner = silent_wgl.config.spinner

        @test layout isa Bonito.Hyperscript.Node
        @test all(node -> node isa Bonito.Hyperscript.Node, masters)
        @test selector[] == 3
        @test silent_wgl isa WGLMakie.WithConfig
        @test silent_wgl.config.resize_to === :parent
        @test silent_spinner isa Bonito.Hyperscript.Node
        @test getfield(silent_spinner, :attrs)["hidden"] === true
        @test getfield(silent_spinner, :attrs)["class"] == "wglmakie-spinner"
        @test_throws ArgumentError wgl_figure(silent_figure; spinner = nothing)
        @test_throws ArgumentError layout_two_columns("left", "right"; ratio = :bad)
        @test_throws ArgumentError layout_two_rows("top", "bottom"; ratio = :bad)
        @test_throws ArgumentError layout_controls_plot(
            "controls",
            "plot";
            footer = "footer",
            footer_placement = :bad
        )
        @test preparation isa Bonito.Hyperscript.Node
        @test preflight_summary((resource,)).phase === :cold
        resource.activate()
        @test preflight_ready(resource_state)
        @test preflight_summary((resource,)).phase === :hot
        @test preflight_summary(()).phase === :hot
        set_preflight!(resource_state, :cold, 0.0, "Cold again")
        PlaygroundHome.activate!((; resources = (resource,)))
        @test preflight_ready(resource_state)
        @test descriptor.render
        @test descriptor.fixture
        @test descriptor.class == ""
        @test only(descriptor.pages).id == "first-view"
        @test_throws ArgumentError webgrid(
            [:area :area; :area :other];
            area = "area",
            other = "other"
        )
        @test_throws ArgumentError webgrid([:left :right]; left = "left")
        @test_throws ArgumentError webgrid(
            reshape([:left], 1, 1);
            left = "left",
            unexpected = "unexpected"
        )
        @test_throws ArgumentError webgrid(reshape(["bad area"], 1, 1); bad = "bad")
        @test_throws ArgumentError preflight_resource(
            id = "Bad resource",
            title = "Bad",
            activate = () -> nothing
        )

        authored_app = App(_ -> layout)
        parent = Session(NoConnection(); asset_server = NoServer())
        io = IOBuffer()
        rendered = Bonito.show_html(io, authored_app; parent)
        html = String(take!(io))
        try
            for marker in (
                "lc-canvas",
                "lc-grid-cell-controls",
                "lc-grid-cell-figure",
                "lc-bonito-markdown",
                "data-math-renderer=\"bonito\"",
                "katex-container",
                "--lc-authored-areas"
            )
                @test occursin(marker, html)
            end
        finally
            close(rendered)
            close(parent)
        end

        @test TemplatePresentation.slug == "starter-template"
        @test length(TemplatePresentation.decks) == 1
        template_deck = only(TemplatePresentation.decks)
        @test template_deck.id == "starter-deck"
        @test length(rendered_deck_pages(template_deck)) == 9
        @test length(TemplatePresentation.resources) == 1
        @test only(TemplatePresentation.resources).id == "starter-resource"
        template_source = read(only(getproperty.(TemplatePresentation.all_decks, :source)), String)
        for layout_name in (
            "layout_single",
            "layout_two_columns",
            "layout_three_columns",
            "layout_top_split",
            "layout_split_bottom",
            "layout_sidebar_main",
            "layout_controls_plot",
            "layout_two_rows",
            "layout_quad"
        )
            @test occursin(layout_name, template_source)
        end
        @test occursin("preflight_resource(", template_source)
        @test occursin("resources = (PREFLIGHT_RESOURCE,)", template_source)
    end

    @testset "domain construction" begin
        cable = CablePage.design(12.5, 8.0)
        values = CablePage.geometry(cable)

        @test cable isa CableDesign
        @test cable.root isa Stack
        @test length(cable.root.items) == 2
        @test cable.root.items[1] isa Group
        @test cable.root.items[2] isa Region
        @test length(cable.geometry.regions) == 2
        @test cable.terminal_order == [:core]
        @test values.core_r_in_m == 0.0
        @test values.core_r_ex_m ≈ 12.5e-3
        @test values.insulation_r_in_m ≈ values.core_r_ex_m
        @test values.insulation_r_ex_m ≈ 20.5e-3
        @test values.core_area_m2 > 0
        @test values.core_area_m2 ≈ π * values.core_r_ex_m^2

        wider_core = CablePage.design(15.0, 8.0)
        thicker_insulation = CablePage.design(12.5, 10.0)
        @test wider_core !== cable
        @test thicker_insulation !== cable
        @test CablePage.geometry(wider_core).core_r_ex_m ≈ 15.0e-3
        @test CablePage.geometry(thicker_insulation).insulation_r_ex_m ≈ 22.5e-3

        for invalid in (0.0, -1.0, Inf, NaN)
            @test_throws DomainError CablePage.design(invalid, 8.0)
            @test_throws DomainError CablePage.design(12.5, invalid)
        end
    end

    @testset "plot construction" begin
        state = Observable(CablePage.design(12.5, 8.0))
        plot = CablePage.cable_figure(state)
        figure = plot.figure
        insulation_plot = plot.insulation_plot
        core_plot = plot.core_plot

        @test figure isa Figure
        @test length(plot.axis.scene.plots) == 2
        @test plot.insulation_circle[].r ≈ 20.5f0
        @test plot.core_circle[].r ≈ 12.5f0
        @test plot.figure.scene.backgroundcolor[] ==
              Makie.to_color(CablePage.PLOT_BACKGROUND)
        @test plot.axis.scene.backgroundcolor[] ==
              Makie.to_color(CablePage.PLOT_BACKGROUND)
        @test plot.insulation_plot.color[] == CablePage.XLPE_COLOR
        @test plot.core_plot.color[] == CablePage.COPPER_COLOR
        @test plot.insulation_plot.strokecolor[] == :black
        @test plot.core_plot.strokecolor[] == :black

        state[] = CablePage.design(20.0, 3.0)
        @test plot.figure === figure
        @test plot.insulation_plot === insulation_plot
        @test plot.core_plot === core_plot
        @test plot.insulation_circle[].r ≈ 23.0f0
        @test plot.core_circle[].r ≈ 20.0f0

        mktempdir() do directory
            output = joinpath(directory, "cable.svg")
            Makie.save(output, figure)
            @test isfile(output)
            @test filesize(output) > 0
        end
    end

    @testset "isolated stateful Julia worker" begin
        worker = ReplPage.start_worker()
        try
            first_events = Any[]
            ReplPage.evaluate!(worker, "x = 21; println(\"worker output\"); 2x") do event
                push!(first_events, event)
            end
            @test ReplPage.worker_running(worker)
            @test any(
                event -> event.kind === :stdout &&
                         String(event.data) == "worker output\n",
                first_events
            )
            @test any(
                event -> event.kind === :result && event.data == "42",
                first_events
            )
            @test last(first_events).kind === :done

            second_events = Any[]
            ReplPage.evaluate!(event -> push!(second_events, event), worker, "x + 1")
            @test getproperty.(second_events, :kind) == [:result, :done]
            @test first(second_events).data == "22"

            expression_events = Any[]
            ReplPage.evaluate!(event -> push!(expression_events, event), worker, ":(1 + 1)")
            ReplPage.evaluate!(
                event -> push!(expression_events, event),
                worker,
                "ans isa Expr"
            )
            @test any(
                event -> event.kind === :result && event.data == "true",
                expression_events
            )

            error_events = Any[]
            ReplPage.evaluate!(event -> push!(error_events, event), worker, "error(\"boom\")")
            @test any(event -> event.kind === :error, error_events)
            @test occursin(
                "boom",
                only(event.data for event in error_events if event.kind === :error)
            )

            suppressed_events = Any[]
            ReplPage.evaluate!(event -> push!(suppressed_events, event), worker, "x + 2;")
            @test getproperty.(suppressed_events, :kind) == [:done]
            @test occursin("julia>", ReplPage.prompt("1 + 1", 1))
        finally
            ReplPage.stop_worker!(worker)
        end
        @test !ReplPage.worker_running(worker)

        session = Session(NoConnection(); asset_server = NoServer())
        page = ReplPage.build(session, nothing)
        page_worker = page.controller.worker
        try
            @test page.editor isa CodeEditor
            @test page.terminal isa TerminalOutput
            @test !isnothing(page_worker)
            page.editor.onchange[] = "page_state = 40; println(\"from page\"); page_state + 2"
            page.run_button.value[] = true
            @test timedwait(() -> !page.controller.busy, 10.0) === :ok
            @test occursin("from page", page.terminal.content[])
            @test occursin("42", page.terminal.content[])
            @test occursin("ready", page.status[])
        finally
            close(session)
        end
        @test page.controller.closed
        @test !ReplPage.worker_running(page_worker)
    end

    @testset "parametric OHL/UGC application case" begin
        network = TransitionPage.transition_network(0.25)
        endpoint_ohl = TransitionPage.transition_network(0.0)
        endpoint_ugc = TransitionPage.transition_network(1.0)

        @test network isa PowerImpedance.NetworkBuilder.NetworkState
        @test propertynames(network.elements) ==
              (:g1, :c1, :c2, :ugc, :ohl, :g4, :tl1, :tl78)
        @test length(network.topology.connections) == 22
        @test network.elements.ugc.element_model.length ≈ 25.0e3
        @test network.elements.ohl.element_model.length ≈ 75.0e3
        @test endpoint_ohl.elements.ugc.element_model.length ≈ 100.0
        @test endpoint_ohl.elements.ohl.element_model.length ≈ 99.9e3
        @test endpoint_ugc.elements.ugc.element_model.length ≈ 99.9e3
        @test endpoint_ugc.elements.ohl.element_model.length ≈ 100.0
        @test TransitionPage.transition_network(0.25) !== network
        @test_throws DomainError TransitionPage.transition_network(NaN)

        response = PowerImpedance.FrequencyResponseResult(
            PowerImpedance.NodalImpedance(),
            :nodal_impedance,
            reshape(ComplexF64[1 + 0im, 10 + 0im, 100 + 0im], 1, 1, 3),
            2π .* [100.0, 1000.0, 5000.0],
            [:B5],
            nothing,
            (;)
        )
        curve = TransitionPage.harmonic_curve(response)
        @test curve.frequencies_hz ≈ [100.0, 1000.0, 5000.0]
        @test curve.magnitude_db ≈ [0.0, 20.0, 40.0]

        impedance_plot = TransitionPage.impedance_figure(response)
        impedance_line = impedance_plot.line
        diagram = PowerImpedanceDiagramExt.diagram(
            network;
            layout = NetworkLayout.Stress(; seed = 7),
            interactive = false
        )
        TransitionPage.use_node_labels!(diagram)
        corridor_edges = TransitionPage.corridor_edge_indices(diagram)
        figure = impedance_plot.figure
        graph_plot = diagram.plots.graph
        line = impedance_line

        @test length(corridor_edges) == 4
        @test length(TransitionPage.corridor_edge_indices(diagram, :ohl)) == 2
        @test length(TransitionPage.corridor_edge_indices(diagram, :ugc)) == 2
        @test !isnothing(diagram.plots.legend)

        next_curve = (;
            frequencies_hz = curve.frequencies_hz,
            magnitude_db = [5.0, 25.0, 45.0]
        )
        fixed_network_model = Ref(:fixed_linearization)
        reference_powerflow = Ref(:reference_powerflow)
        runtime = TransitionPage.TransitionRuntime(
            network,
            fixed_network_model,
            reference_powerflow,
            1,
            impedance_plot,
            impedance_line,
            diagram,
            Dict(50 => curve, 60 => next_curve),
            0.50,
            ReentrantLock()
        )
        TransitionPage.update!(
            runtime,
            0.60;
            frequency_range = (100.0, 5.0e3, 3)
        )

        @test runtime.network === network
        @test runtime.network_model === fixed_network_model
        @test runtime.reference_powerflow === reference_powerflow
        @test runtime.impedance_plot.figure === figure
        @test runtime.impedance_line === line
        @test runtime.diagram.plots.graph === graph_plot
        @test runtime.share == 0.60
        @test network.elements.ohl.element_model.length ≈ 40.0e3
        @test network.elements.ugc.element_model.length ≈ 60.0e3
        @test last.(line[1][]) == next_curve.magnitude_db
        @test all(
            ==(TransitionPage.OHL_COLOR),
            diagram.plots.graph.edge_color[][
                TransitionPage.corridor_edge_indices(diagram, :ohl)
            ]
        )
        @test all(
            ==(TransitionPage.UGC_COLOR),
            diagram.plots.graph.edge_color[][
                TransitionPage.corridor_edge_indices(diagram, :ugc)
            ]
        )
        @test all(==(6.0f0), diagram.plots.graph.edge_width[][corridor_edges])
        b4_key = TransitionPage.diagram_bus_key(diagram, :B4)
        bx_key = TransitionPage.diagram_bus_key(diagram, :BX)
        b5_key = TransitionPage.diagram_bus_key(diagram, :B5)
        @test diagram.positions[bx_key] ≈
              diagram.positions[b4_key] +
              0.60 * (diagram.positions[b5_key] - diagram.positions[b4_key])
        @test_throws DomainError TransitionPage.set_corridor_lengths!(network, NaN)

        mktempdir() do directory
            output = joinpath(directory, "network.svg")
            Makie.save(output, diagram.figure)
            @test isfile(output)
            @test filesize(output) > 0
        end
    end

    @testset "ICHQP deterministic OHL/UGC boundaries" begin
        bounds = ApplicationCasePage.transition_bounds(0.50, 5.0)
        @test bounds.lower ≈ 0.475
        @test bounds.nominal ≈ 0.50
        @test bounds.upper ≈ 0.525
        @test ApplicationCasePage.transition_bounds(1.0, 10.0).upper == 1.0
        @test_throws DomainError ApplicationCasePage.transition_bounds(0.5, NaN)
        @test_throws DomainError ApplicationCasePage.transition_bounds(0.5, -1.0)

        frequencies = [100.0, 1000.0, 5000.0]
        lower_curve = (; frequencies_hz = frequencies, magnitude_db = [38.0, 45.0, 52.0])
        nominal_curve = (; frequencies_hz = frequencies, magnitude_db = [40.0, 47.0, 54.0])
        upper_curve = (; frequencies_hz = frequencies, magnitude_db = [42.0, 49.0, 56.0])
        bundle = (;
            bounds,
            curves = (;
                lower = lower_curve,
                nominal = nominal_curve,
                upper = upper_curve
            )
        )
        plot = ApplicationCasePage.impedance_bounds_figure(bundle, 5.0)
        identities = (
            plot.figure,
            plot.axes.lower,
            plot.lines.lower,
            plot.lines.nominal,
            plot.lines.upper
        )
        shifted = (;
            bounds = ApplicationCasePage.transition_bounds(0.60, 5.0),
            curves = (;
                lower = merge(lower_curve, (; magnitude_db = [39.0, 46.0, 53.0])),
                nominal = merge(nominal_curve, (; magnitude_db = [41.0, 48.0, 55.0])),
                upper = merge(upper_curve, (; magnitude_db = [43.0, 50.0, 57.0]))
            )
        )
        ApplicationCasePage.update_impedance_bounds!(plot, shifted, 5.0)
        @test identities == (
            plot.figure,
            plot.axes.lower,
            plot.lines.lower,
            plot.lines.nominal,
            plot.lines.upper
        )
        @test last.(plot.lines.nominal[1][]) == [41.0, 48.0, 55.0]

        source = read(
            joinpath(
                @__DIR__,
                "..",
                "presentations",
                "ICHQP2026",
                "content",
                "050_applicationcase.jl"
            ),
            String
        )
        @test !occursin("Measurements", source)
        @test occursin("async_latest(input, 1)", source)
        @test occursin("Case.harmonic_response", source)
        @test occursin("Case.update_corridor_diagram!", source)
        @test occursin(
            "network_model = Case.linearized_network(network, prepared.powerflow)",
            source
        )
    end

    @testset "fundamentals numerical kernels" begin
        skin_model = FundamentalsPage.SkinEffectModel(
            grid_points = 32,
            radial_points = 128
        )
        skin_cache = FundamentalsPage.build_skin_cache(skin_model)
        skin = FundamentalsPage.solve_skin_state(
            1.0e3,
            5.8e7,
            skin_model,
            skin_cache
        )
        @test size(skin.field_db) == (32, 32)
        @test length(skin.profile_db) == 128
        @test isfinite(skin.skin_depth)
        @test isfinite(skin.ac_over_dc)
        @test isfinite(real(skin.impedance))
        @test isfinite(imag(skin.impedance))
        @test_throws ArgumentError FundamentalsPage.solve_skin_state(
            0.0,
            5.8e7,
            skin_model,
            skin_cache
        )

        earth_model = FundamentalsPage.EarthReturnModel(
            local_ny = 32,
            local_nz = 32,
            earth_ny = 32,
            earth_nz = 32,
            spectrum_points = 128
        )
        earth_cache = FundamentalsPage.build_earth_cache(earth_model)
        earth = FundamentalsPage.solve_earth_state(
            50.0,
            500.0,
            100.0,
            10.0,
            0.0,
            earth_model,
            earth_cache
        )
        @test size(earth.field_db) == (32, 32)
        @test size(earth.current_db) == (32, 32)
        @test isfinite(earth.skin_depth)
        @test isfinite(earth.conduction_ratio)
        @test isfinite(real(earth.earth_self))
        @test isfinite(imag(earth.earth_self))
        @test isfinite(earth.maximum_current)
        @test_throws ArgumentError FundamentalsPage.solve_earth_state(
            50.0,
            500.0,
            -1.0,
            10.0,
            0.0,
            earth_model,
            earth_cache
        )
    end

    @testset "package-backed cable constants dashboard" begin
        design = CableConstantsPage.build_design(4.7e-3, 8.0e-3, 2)
        @test design isa CableDesign
        @test eltype(design) === Float64
        @test design.terminal_order == [:core]
        @test count(
            region -> region.terminal === :core,
            design.geometry.regions
        ) == CableConstantsPage.strand_count(2) == 19

        @test CableConstantsPage.strand_count(1) == 7
        @test CableConstantsPage.strand_count(3) == 37
        @test CableConstantsPage.strand_count(10) == 331
        @test last(CableConstantsPage.STRAND_LAYER_RANGE) == 10
        @test last(CableConstantsPage.UNCERTAINTY_RANGE_PCT) == 10.0
        three_layer_design = CableConstantsPage.build_design(4.7e-3, 8.0e-3, 3)
        @test count(
            region -> region.terminal === :core,
            three_layer_design.geometry.regions
        ) == 37
        ten_layer_design = CableConstantsPage.build_design(4.7e-3, 8.0e-3, 10)
        @test count(
            region -> region.terminal === :core,
            ten_layer_design.geometry.regions
        ) == 331

        projection = CableConstantsPage.preview_projection(design)
        @test length(projection.strand_centers) == 19
        @test length(projection.strand_polygons) == 19
        @test length(projection.layers) == 1
        @test projection.strand_diameter_mm ≈ 4.7 rtol = 2.0e-4
        @test projection.outer_diameter_mm ≈ 39.5

        boundaries = CableConstantsPage.boundary_projection(design)
        @test length(boundaries) == 21
        wire_envelope = CableConstantsPage.envelope_projection(
            design,
            CableConstantsPage.DEFAULT_INPUT,
            :wire_diameter
        )
        insulation_envelope = CableConstantsPage.envelope_projection(
            design,
            CableConstantsPage.DEFAULT_INPUT,
            :insulation_thickness
        )
        @test length(wire_envelope.affected) == 21
        @test length(insulation_envelope.affected) == 1
        @test wire_envelope.vertex_count > insulation_envelope.vertex_count > 0
        @test wire_envelope.face_count == wire_envelope.vertex_count
        zero_envelope = CableConstantsPage.envelope_projection(
            design,
            merge(
                CableConstantsPage.DEFAULT_INPUT,
                (; wire_uncertainty = 0.0)
            ),
            :wire_diameter
        )
        @test isempty(zero_envelope.affected)
        @test zero_envelope.vertex_count == 0

        measured_wire = Measurements.measurement(4.7e-3, 4.7e-5)
        measured_insulation = Measurements.measurement(8.0e-3, 8.0e-5)
        uncertain_design = CableConstantsPage.build_design(
            measured_wire,
            measured_insulation,
            2
        )
        @test eltype(uncertain_design) === Measurements.Measurement{Float64}

        diagnostics = CableConstantsPage.conductor_diagnostics(uncertain_design)
        expected_area = 19 * π * (measured_wire / 2)^2
        @test diagnostics.conductor_area ≈ expected_area
        @test diagnostics.dc_resistance ≈
              CableConstantsPage.COPPER.rho / expected_area
        @test Measurements.uncertainty(diagnostics.conductor_area) > 0
        @test Measurements.uncertainty(diagnostics.dc_resistance) > 0

        result = CableConstantsPage.compute_constants(
            CableConstantsPage.DEFAULT_INPUT
        )
        @test result.design isa CableDesign
        @test result.constants isa CableConstants
        @test result.input.frequency == 1.0e3
        @test result.input.strand_layers == 2
        @test all(
            number -> number isa Measurements.Measurement{Float64},
            values(result.observations)
        )
        @test all(
            number -> Measurements.uncertainty(number) > 0,
            values(result.observations)
        )
        @test result.diagnostics.conductor_area ≈ diagnostics.conductor_area
        @test occursin(
            "mm²",
            CableConstantsPage.display_physical_measurement(
                result.diagnostics.conductor_area,
                1.0e6,
                "mm²"
            )
        )
        @test occursin(
            "Ω/km",
            CableConstantsPage.display_physical_measurement(
                result.diagnostics.dc_resistance,
                1.0e3,
                "Ω/km"
            )
        )
        frequency_only = merge(
            CableConstantsPage.DEFAULT_INPUT,
            (; frequency = 2.0e3)
        )
        @test CableConstantsPage.numerical_design_key(frequency_only) ==
              CableConstantsPage.numerical_design_key(
            CableConstantsPage.DEFAULT_INPUT
        )
        another_layer = merge(
            CableConstantsPage.DEFAULT_INPUT,
            (; strand_layers = 3)
        )
        @test CableConstantsPage.numerical_design_key(another_layer) !=
              CableConstantsPage.numerical_design_key(
            CableConstantsPage.DEFAULT_INPUT
        )
        @test occursin(
            "Ω/km",
            CableConstantsPage.display_measurement(result.observations.R, :R)
        )
        @test occursin(
            "mH/km",
            CableConstantsPage.display_measurement(result.observations.L, :L)
        )
        @test occursin(
            "μF/km",
            CableConstantsPage.display_measurement(result.observations.C, :C)
        )

        preview = CableConstantsPage.preview_figure(design)
        figure = preview.figure
        axis = preview.axis
        layer_plot = preview.layer_plot
        strand_plot = preview.strand_plot
        clone_layer_plot = preview.clone_layer_plot
        clone_strand_plot = preview.clone_strand_plot
        next_design = CableConstantsPage.build_design(6.0e-3, 10.0e-3, 3)
        next_projection = CableConstantsPage.update_preview!(preview, next_design)
        @test preview.figure === figure
        @test preview.axis === axis
        @test preview.layer_plot === layer_plot
        @test preview.strand_plot === strand_plot
        @test preview.clone_layer_plot === clone_layer_plot
        @test preview.clone_strand_plot === clone_strand_plot
        @test length(axis.scene.plots) == 4
        @test length(next_projection.strand_centers) == 37
        @test length(next_projection.strand_polygons) == 37
        @test strand_plot[1][] == next_projection.strand_polygons
        translated = CableConstantsPage.translated_preview_projection(
            next_projection,
            Point2f(CableConstantsPage.CLONE_LATERAL_OFFSET_MM, 0)
        )
        @test clone_layer_plot[1][] == translated.layers
        @test clone_strand_plot[1][] == translated.strand_polygons
        @test all(
            isapprox(
                translated.strand_centers[index][1] -
                next_projection.strand_centers[index][1],
                1000.0
            ) for index in eachindex(next_projection.strand_centers)
        )
        @test next_projection.outer_diameter_mm ≈ 62.0
        @test preview.radius_limit_mm[] ≈
              CableConstantsPage.PREVIEW_MARGIN * 31.0

        uncertainty_preview = CableConstantsPage.uncertainty_figure(
            design,
            CableConstantsPage.DEFAULT_INPUT
        )
        uncertainty_plots = (
            uncertainty_preview.layer_plot,
            uncertainty_preview.strand_plot,
            uncertainty_preview.wire_envelope_plot,
            uncertainty_preview.insulation_envelope_plot,
            uncertainty_preview.outline_plot
        )
        @test length(uncertainty_preview.axis.scene.plots) == 5
        CableConstantsPage.update_uncertainty_nominal!(
            uncertainty_preview,
            next_design,
            3
        )
        CableConstantsPage.update_uncertainty_source!(
            uncertainty_preview,
            merge(CableConstantsPage.DEFAULT_INPUT, (;
                strand_layers = 3,
                wire_diameter = 6.0e-3,
                insulation_thickness = 10.0e-3
            )),
            :wire_diameter
        )
        @test uncertainty_plots == (
            uncertainty_preview.layer_plot,
            uncertainty_preview.strand_plot,
            uncertainty_preview.wire_envelope_plot,
            uncertainty_preview.insulation_envelope_plot,
            uncertainty_preview.outline_plot
        )
        @test_throws DomainError CableConstantsPage.build_design(0.0, 8.0e-3, 2)
        @test_throws DomainError CableConstantsPage.build_design(4.7e-3, Inf, 2)
        @test_throws DomainError CableConstantsPage.build_design(4.7e-3, 8.0e-3, -1)
    end

    @testset "package-backed line-parameter dashboards" begin
        frequencies = [0.5, 50.0, 2500.0]
        model = LineParametersPage.build_problem(
            LineParametersPage.DEFAULT_INPUT;
            frequencies
        )
        @test model.cable_1 isa CableDesign
        @test model.cable_2 isa CableDesign
        @test model.cable_1.terminal_order == [:core]
        @test model.cable_2.terminal_order == [:core]
        @test length(model.cable_1.geometry.regions) == 2
        @test length(model.cable_2.geometry.regions) == 2
        @test model.system.terminal_order == [
            (cable = 1, terminal = :core),
            (cable = 2, terminal = :core)
        ]
        @test model.system.connection_order == [1, 2]
        @test Measurements.value(model.system.positions[1].y) == -1.0
        @test Measurements.value(model.system.positions[2].y) == -1.0
        @test Measurements.value(
            model.system.positions[2].x - model.system.positions[1].x
        ) ≈ LineParametersPage.DEFAULT_INPUT.lateral_distance
        @test Measurements.uncertainty(model.system.line_length) > 0
        @test Measurements.uncertainty(model.earth.layers[2].rho) > 0

        result = LineParametersPage.compute_line_parameters(
            LineParametersPage.DEFAULT_INPUT;
            frequencies
        )
        @test result.parameters isa LineParameters
        @test basis(result.parameters.Z) === :total
        @test basis(result.parameters.Y) === :total
        @test size(result.parameters.Z) == (2, 2, 3)
        @test size(result.parameters.Y) == (2, 2, 3)
        @test result.frequencies == frequencies
        for element in (:self, :mutual), quantity in (:R, :L, :C, :G)
            curve = getproperty(getproperty(result.curves, element), quantity)
            @test length(curve.nominal) == length(frequencies)
            @test all(isfinite, curve.nominal)
            @test all(curve.lower .<= curve.nominal)
            @test all(curve.nominal .<= curve.upper)
        end
        @test any(
            getproperty(result.curves.self, quantity).lower !=
            getproperty(result.curves.self, quantity).upper
        for quantity in (:R, :L, :C, :G))

        self_plot = LineParametersPage.line_parameter_figure(result, :self)
        mutual_plot = LineParametersPage.line_parameter_figure(result, :mutual)
        @test length(self_plot.axes) == 4
        @test length(mutual_plot.axes) == 4
        identities = (
            self_plot.figure,
            self_plot.axes[1],
            self_plot.curves.R.line_plot,
            self_plot.curves.R.band_plot
        )
        LineParametersPage.update_line_parameter_figure!(self_plot, result)
        @test identities == (
            self_plot.figure,
            self_plot.axes[1],
            self_plot.curves.R.line_plot,
            self_plot.curves.R.band_plot
        )
        @test_throws DomainError LineParametersPage.build_cable(
            "invalid",
            0.0,
            8.0e-3
        )
        @test_throws DomainError LineParametersPage.measured(1.0, 101.0)
    end

    @testset "presentation and deck discovery" begin
        @test getproperty.(PRESENTATION_DESCRIPTORS, :slug) == [
            "ICHQP2026",
            "live-manual"
        ]
        @test ICHQP2026.title == "ICHQP2026"
        @test getproperty.(ICHQP2026.all_decks, :id) ==
              [
            "opening",
            "fundamentals",
            "cable-constants",
            "line-parameters",
            "application-case",
            "conclusion"
        ]
        @test getproperty.(ICHQP2026.decks, :id) ==
              [
            "opening",
            "fundamentals",
            "cable-constants",
            "line-parameters",
            "application-case",
            "conclusion"
        ]
        opening_deck = find_deck("opening", ICHQP2026.decks)
        fundamentals_deck = find_deck("fundamentals", ICHQP2026.decks)
        cable_constants_deck = find_deck("cable-constants", ICHQP2026.decks)
        line_parameters_deck = find_deck("line-parameters", ICHQP2026.decks)
        application_case_deck = find_deck("application-case", ICHQP2026.decks)
        conclusion_deck = find_deck("conclusion", ICHQP2026.decks)
        @test opening_deck.render
        @test opening_deck.expand_navigation
        @test fundamentals_deck.id == "fundamentals"
        @test getproperty.(opening_deck.pages, :id) ==
              ("title", "introduction")
        @test first(opening_deck.pages).class == "lc-title-page"
        @test occursin("lc-fill-page", last(opening_deck.pages).class)
        @test getproperty.(fundamentals_deck.pages, :id) ==
              (
            "sources-of-uncertainties",
            "internal-impedances",
            "earth-return",
            "uncertainty-quantification"
        )
        @test all(
            page -> occursin("lc-application-slide", page.class),
            fundamentals_deck.pages[2:3]
        )
        @test all(
            page -> occursin("lc-fluid-type", page.class),
            (fundamentals_deck.pages[1], fundamentals_deck.pages[4])
        )
        @test isempty(opening_deck.resources)
        @test length(fundamentals_deck.resources) == 1
        fundamentals_resource = only(fundamentals_deck.resources)
        @test fundamentals_resource.id == "fundamentals-dashboards"
        @test fundamentals_resource.state === FundamentalsPage.PREFLIGHT_STATE
        @test getproperty.(cable_constants_deck.pages, :id) ==
              (
            "geometry-uncertainty-envelopes",
            "cable-constants-under-uncertainty"
        )
        @test length(cable_constants_deck.resources) == 1
        cable_constants_resource = only(cable_constants_deck.resources)
        @test cable_constants_resource.id == "cable-constants-dashboard"
        @test cable_constants_resource.state === CableConstantsPage.PREFLIGHT_STATE
        @test getproperty.(line_parameters_deck.pages, :id) ==
              ("self-line-parameters", "mutual-line-parameters")
        @test length(line_parameters_deck.resources) == 1
        line_parameters_resource = only(line_parameters_deck.resources)
        @test line_parameters_resource.id == "line-parameters-dashboard"
        @test line_parameters_resource.state === LineParametersPage.PREFLIGHT_STATE
        @test getproperty.(application_case_deck.pages, :id) ==
              ("parametric-ohl-ugc-transition",)
        @test length(application_case_deck.resources) == 1
        application_case_resource = only(application_case_deck.resources)
        @test application_case_resource.id == "ichqp2026-application-case"
        @test application_case_resource.state === ApplicationCasePage.PREFLIGHT_STATE
        @test getproperty.(conclusion_deck.pages, :id) ==
              ("concluding-remarks",)
        @test isempty(conclusion_deck.resources)
        @test ICHQP2026.resources == [
            fundamentals_resource,
            cable_constants_resource,
            line_parameters_resource,
            application_case_resource
        ]
        cable_constants_source = read(cable_constants_deck.source, String)
        @test occursin("LCM.DataModel.preview_shapes", cable_constants_source)
        @test occursin("LCM.CableConstants(", cable_constants_source)
        @test occursin("LCM.observe(constants, LCM.R)", cable_constants_source)
        @test occursin("count = iszero(layer) ? 1 : 6layer", cable_constants_source)
        @test occursin("COPPER.rho / conductor_area", cable_constants_source)
        @test occursin("boundary_strip_mesh", cable_constants_source)
        @test occursin("perturbed_designs", cable_constants_source)
        @test occursin("Makie.update!", cable_constants_source)
        @test occursin("projection.strand_polygons", cable_constants_source)
        @test !occursin("scatter!", cable_constants_source)
        @test occursin("Threads.@spawn compute_constants", cable_constants_source)
        @test !occursin("bessel", lowercase(cable_constants_source))
        @test !occursin("montecarlo", lowercase(cable_constants_source))
        line_parameters_source = read(line_parameters_deck.source, String)
        @test occursin("LCM.LineParametersProblem(", line_parameters_source)
        @test occursin("LCM.compute(", line_parameters_source)
        @test occursin("LCM.observe(parameters", line_parameters_source)
        @test occursin("output_basis = :total", line_parameters_source)
        @test occursin("Threads.@spawn compute_line_parameters", line_parameters_source)
        @test !occursin("pollaczek", lowercase(line_parameters_source))
        @test !occursin("bessel", lowercase(line_parameters_source))
        @test preflight_summary(ICHQP2026.resources).phase in
              (:cold, :preparing, :hot)
        @test all(
            isfile(joinpath(ICHQP2026.directory, "assets", asset))
        for asset in (
            "ETCH_LOGO_RGB_NEG.svg",
            "ENERGYVILLE-LOGO.svg",
            "kul_logo.svg",
            "skeffect.png",
            "earthreturn.png",
            "cable_dark_mode.svg"
        ))
        @test LiveManual.title == "Live manual"
        @test LiveManual.content == abspath(joinpath(
            @__DIR__, "..", "presentations", "live-manual", "content"
        ))
        @test parentmodule(CablePage) === LiveManual.namespace
        @test length(DECK_DESCRIPTORS) == 3
        @test getproperty.(DECK_DESCRIPTORS, :id) == [
            "core-and-insulation",
            "parametric-ohl-ugc-transition",
            "live-julia-workspace"
        ]
        @test all(deck -> deck.render, DECK_DESCRIPTORS)
        @test allunique(deck.id for deck in DECK_DESCRIPTORS)
        @test allunique(deck.source for deck in DECK_DESCRIPTORS)
        transition_deck = find_deck("parametric-ohl-ugc-transition", DECK_DESCRIPTORS)
        @test length(transition_deck.resources) == 1
        @test only(transition_deck.resources).id == "powerflow-linearization"
        @test only(transition_deck.resources).state === TransitionPage.PREFLIGHT_STATE
        @test LiveManual.resources == collect(transition_deck.resources)
        @test preflight_summary(LiveManual.resources).phase in (:cold, :preparing, :hot)
        @test getproperty.(rendered_deck_pages(transition_deck), :id) ==
              ("corridor", "b5-impedance", "linearization")
        @test all(
            startswith(basename(deck.source), lpad(string(deck.order), 3, '0'))
        for deck in DECK_DESCRIPTORS)
        @test find_deck("core-and-insulation", DECK_DESCRIPTORS).export_figure isa
              Function

        hidden_cable = map(DECK_DESCRIPTORS) do deck
            deck.id == "core-and-insulation" ? merge(deck, (; render = false)) : deck
        end
        selected = select_rendered_decks(hidden_cable)
        @test getproperty.(selected, :id) ==
              ["parametric-ohl-ugc-transition", "live-julia-workspace"]
        @test_throws ArgumentError select_rendered_decks(
            map(deck -> merge(deck, (; render = false)), DECK_DESCRIPTORS)
        )

        hidden_page_deck = merge(transition_deck,
            (;
                pages = map(transition_deck.pages) do page
                page.id == "b5-impedance" ? merge(page, (; render = false)) : page
            end
            ))
        @test getproperty.(rendered_deck_pages(hidden_page_deck), :id) ==
              ("corridor", "linearization")

        mktempdir() do directory
            write(
                joinpath(directory, "020_visible.jl"),
                """module VisibleDiscoveryFixture
                setup(value) = value
                build(session, state) = (; body = state)
                const DECK = (id = \"visible-fixture\", group = \"Fixtures\", title = \"Visible\", order = 20, render = true, class = \"fixture\", setup = setup, pages = ((id = \"visible-page\", title = \"Visible page\", render = true, class = \"\", build = build),))
                end
                VisibleDiscoveryFixture.DECK
                """
            )
            write(
                joinpath(directory, "010_hidden.jl"),
                """module HiddenDiscoveryFixture
                setup(value) = value
                build(session, state) = (; body = state)
                const DECK = (id = \"hidden-fixture\", group = \"Fixtures\", title = \"Hidden\", order = 10, render = false, class = \"fixture\", setup = setup, pages = ((id = \"hidden-page\", title = \"Hidden page\", render = true, class = \"\", build = build),))
                end
                HiddenDiscoveryFixture.DECK
                """
            )
            discovered = load_deck_descriptors(directory)
            @test getproperty.(discovered, :id) ==
                  ["hidden-fixture", "visible-fixture"]
            @test getproperty.(select_rendered_decks(discovered), :id) ==
                  ["visible-fixture"]
        end
    end

    @testset "application DOM" begin
        assets = app_assets()
        transition_setup_count = Ref(0)
        test_preflight_state = preflight_state(
            phase = :hot,
            progress = 1.0,
            label = "Fixture preflight hot"
        )
        test_preflight_resource = preflight_resource(
            id = "fixture-preflight",
            title = "Fixture preflight",
            state = test_preflight_state,
            activate = () -> nothing
        )
        transition_setup = _ -> begin
            transition_setup_count[] += 1
            return nothing
        end
        corridor_stub = (_,
            _) -> (;
            body = (
            DOM.h1("Corridor composition"),
            DOM.div(
                DOM.div(; id = "ugc-share-control"),
                DOM.div(; id = "transition-network-diagram");
                class = "lc-interactive-grid lc-transition-grid"
            )
        )
        )
        impedance_stub = (_,
            _) -> (;
            body = (
            DOM.h1("Driving-point impedance at B5"),
            DOM.div(; id = "b5-impedance-plot")
        )
        )
        linearization_stub = (_,
            _) -> (;
            body = (
            DOM.h1("Operating point and linearization"),
            DOM.button(; id = "recache-powerflow-control")
        )
        )
        repl_stub = (_,
            _) -> (;
            body = (
            DOM.h1("Live Julia workspace"),
            DOM.div(
                DOM.button(; id = "live-julia-run"),
                DOM.button(; id = "live-julia-interrupt"),
                DOM.button(; id = "live-julia-restart"),
                DOM.button(; id = "live-julia-clear"),
                DOM.div(; id = "live-julia-editor"),
                DOM.div(; id = "live-julia-console"),
                DOM.div(; id = "live-julia-terminal");
                class = "lc-repl-grid"
            )
        )
        )
        test_decks = map(RENDERED_DECKS) do deck
            if deck.id == "parametric-ohl-ugc-transition"
                pages = map(deck.pages) do page
                    build = page.id == "corridor" ? corridor_stub :
                            page.id == "b5-impedance" ? impedance_stub :
                            linearization_stub
                    merge(page, (; build))
                end
                merge(deck, (;
                    setup = transition_setup,
                    pages,
                    resources = (test_preflight_resource,)
                ))
            elseif deck.id == "live-julia-workspace"
                pages = map(page -> merge(page, (; build = repl_stub)), deck.pages)
                merge(deck, (; pages))
            else
                deck
            end
        end
        test_presentation = merge(
            LiveManual,
            (;
                all_decks = test_decks,
                decks = test_decks,
                resources = presentation_resources(test_decks)
            )
        )
        test_presentations = [test_presentation]
        app = playground_app(assets; presentations = test_presentations)
        routes = playground_routes(assets; presentations = test_presentations)
        @test app isa App
        @test isnothing(app.loading_page)
        @test Set(keys(routes.routes)) == Set((
            "/",
            "/presentations/live-manual/core-and-insulation",
            "/presentations/live-manual/parametric-ohl-ugc-transition",
            "/presentations/live-manual/live-julia-workspace"
        ))
        @test app.indicator.size == 8
        @test app.indicator.top == "20px"
        @test app.indicator.right == "16px"
        @test assets.app_css isa String
        @test occursin(".lc-documenter", assets.app_css)
        @test length(RENDERED_DECKS) == 3
        @test allunique(deck.id for deck in RENDERED_DECKS)
        @test allunique(deck.title for deck in RENDERED_DECKS)
        @test last(RENDERED_DECKS).title == "Live Julia workspace"

        theme = read(joinpath(@__DIR__, "..", "assets", "theme.css"), String)
        @test occursin("--lc-bg: #1f2424;", theme)
        @test occursin("--lc-sidebar-bg: #282f2f;", theme)
        @test occursin("--lc-sidebar-width: 18rem;", theme)
        @test occursin("grid-template-columns: var(--lc-sidebar-width)", theme)
        @test occursin("height: 100dvh;", theme)
        @test occursin("height: 100% !important;", theme)
        @test occursin("position: fixed;", theme)
        @test occursin("--lc-type-unit: clamp", theme)
        @test occursin("var(--lc-viewport-width, 100%)", theme)
        @test occursin("var(--lc-viewport-height, 100%)", theme)
        @test occursin(".lc-preparation-progress", theme)
        @test occursin(".lc-home-activity", theme)
        @test occursin(
            "--lc-effective-rows: minmax(0, 3fr) minmax(0, 2fr);",
            theme
        )
        @test occursin("overscroll-behavior: contain;", theme)
        @test occursin("overflow-x: hidden;", theme)
        @test occursin("margin: 0 !important;", theme)
        @test occursin("padding-right: 1.5rem;", theme)
        @test occursin("@media print", theme)
        @test occursin(".lc-documenter:fullscreen", theme)
        @test occursin(".lc-documenter.is-sidebar-collapsed", theme)
        @test occursin(".lc-navbar-actions", theme)
        @test occursin(".lc-repl-slide", theme)
        @test !occursin("reveal", lowercase(theme))
        @test !occursin("linear-gradient", theme)
        @test !occursin("backdrop-filter", theme)

        source = read(joinpath(@__DIR__, "..", "app.jl"), String)
        @test occursin("PageAuthoring.jl", source)
        @test occursin("load_presentation_descriptors()", source)
        @test occursin("joinpath(directory, \"content\")", source)
        @test occursin("deck -> deck.render", source)
        @test occursin("page -> page.render", source)
        @test occursin("function playground_routes", source)
        @test occursin("presentation_route(presentation, deck)", source)
        @test occursin("state = deck.setup(session)", source)
        @test occursin("window.addEventListener(\"hashchange\"", source)
        @test occursin("activateDeckPage", source)
        @test occursin("window.location.assign", source)
        @test occursin("root.requestFullscreen", source)
        @test occursin("window.sessionStorage", source)
        @test occursin("linecable:playground:presentation-mode", source)
        @test occursin("setFocusMode(readPresentationMode(), false)", source)
        @test !occursin("fallbackFocus", source)
        @test occursin("window.print()", source)
        @test occursin("event.pointerType === \"touch\"", source)
        @test occursin("window.visualViewport", source)
        @test occursin("window.devicePixelRatio", source)
        @test occursin("document.documentElement.clientWidth", source)
        @test occursin("resolution:", source)
        @test occursin("window.setTimeout(applyViewport, 240)", source)
        @test occursin("function deck_content_app", source)
        @test occursin("function loading_deck_content", source)
        @test occursin("selected_content = Observable{Any}", source)
        @test occursin("selected_content[] = loading_deck_content(target)", source)
        @test occursin("selectedDeck.notify(targetDeckId)", source)
        @test occursin("activate_resources!(resources)", source)
        @test occursin("function presentation_resources", source)
        @test occursin("function preflight_observable", source)
        @test !occursin("loading_page =", source)
        @test occursin("new MutationObserver", source)
        @test occursin("window.history.pushState", source)
        @test occursin("window.addEventListener(\"popstate\"", source)
        @test !occursin("ES6Module", source)
        @test !occursin("Reveal", source)
        @test !occursin("dispatchEvent(new Event(\"resize\"))", source)
        @test !occursin("wglmakie_screen", source)

        cable_source = read(
            joinpath(
                @__DIR__, "..", "presentations", "live-manual", "content",
                "010_core_and_insulation.jl"
            ),
            String
        )
        @test occursin("module CableDesignPage", cable_source)
        @test occursin("using ..PageAuthoring", cable_source)
        @test occursin("const DECK = deck_descriptor(", cable_source)
        @test occursin("render = true", cable_source)
        @test !occursin("DOM.", cable_source)
        @test !occursin("LineCableModels.PlotBuilder", cable_source)
        @test !occursin("preview(", cable_source)

        case_source = read(
            joinpath(
                @__DIR__, "..", "presentations", "live-manual", "content",
                "020_ohl_ugc_transition.jl"
            ),
            String
        )
        @test occursin("module OHLUGCTransitionPage", case_source)
        @test occursin("using ..PageAuthoring", case_source)
        @test occursin("const DECK = deck_descriptor(", case_source)
        @test occursin("setup = setup", case_source)
        @test occursin("id = \"corridor\"", case_source)
        @test occursin("id = \"b5-impedance\"", case_source)
        @test occursin("id = \"linearization\"", case_source)
        @test occursin("render = true", case_source)
        @test occursin("``Z_{B_5,B_5}(f)``", case_source)
        @test !occursin(raw"$Z_{B_5,B_5}(f)$", case_source)
        @test !occursin("DOM.", case_source)
        @test occursin("async_latest(requested_share, 1)", case_source)
        @test occursin("LinearizationProblem(network, powerflow)", case_source)
        @test occursin("cached_reference_powerflow(; force = true)", case_source)
        @test occursin("set_corridor_lengths!(runtime.network", case_source)
        @test occursin("harmonic_response(runtime.network_model", case_source)
        @test occursin("Bonito.Button(", case_source)
        @test occursin("recache-powerflow-control", case_source)
        @test occursin("using ..PowerImpedanceDiagramExt", case_source)
        @test occursin("PowerImpedanceDiagramExt.diagram", case_source)
        @test occursin("Channel{NamedTuple}(16)", case_source)
        @test occursin(
            "Threads.@spawn compute_prepared_case(updates)",
            case_source
        )
        @test !occursin("addprocs(", case_source)
        @test !occursin("RemoteChannel", case_source)
        @test occursin(":powerflow,", case_source)
        @test occursin("Solving the reference power flow", case_source)
        @test occursin("preflight_resource(", case_source)
        @test occursin("id = \"powerflow-linearization\"", case_source)
        @test occursin("state = PREFLIGHT_STATE", case_source)
        @test occursin("network = deepcopy(prepared.network)", case_source)
        @test occursin(
            "network_model = linearized_network(network, prepared.powerflow)",
            case_source
        )
        @test !occursin(
            "deepcopy((prepared.network, prepared.network_model))",
            case_source
        )
        @test !occursin("PowerImpedance.plot(", case_source)
        @test !occursin("Observable{NetworkState}", case_source)

        repl_source = read(
            joinpath(
                @__DIR__, "..", "presentations", "live-manual", "content",
                "030_live_julia.jl"
            ),
            String
        )
        worker_source = read(
            joinpath(
                @__DIR__, "..", "presentations", "live-manual", "support",
                "repl_worker.jl"
            ),
            String)
        @test occursin("module LiveJuliaPage", repl_source)
        @test occursin("using ..PageAuthoring", repl_source)
        @test occursin("const DECK = deck_descriptor(", repl_source)
        @test occursin("CodeEditor(", repl_source)
        @test occursin("TerminalOutput(", repl_source)
        @test occursin("session.on_close", repl_source)
        @test occursin("render = true", repl_source)
        @test !occursin("Core.eval", repl_source)
        @test occursin("Core.eval(workspace", worker_source)
        @test occursin("Meta.parseall", worker_source)
        @test occursin("MAX_OUTPUT_BYTES", worker_source)

        project = read(joinpath(@__DIR__, "..", "Project.toml"), String)
        manifest = read(joinpath(@__DIR__, "..", "Manifest.toml"), String)
        line_cable_models_revision =
            "b2b074712f291fac89fd2f4936724c565afa3953"
        public_revision = "457ba27831a3841c97e14b9c832840390df946f4"
        diagram_revision = "71183fe1251cd178bc6c8704594d197d4d988414"
        @test occursin(
            "PowerImpedance = {url = \"https://github.com/Electa-Git/PowerImpedance.jl.git\"",
            project
        )
        @test occursin("Distributed =", project)
        @test occursin("Serialization =", project)
        @test occursin("rev = \"$public_revision\"", project)
        @test occursin("repo-rev = \"$public_revision\"", manifest)
        @test occursin("rev = \"$line_cable_models_revision\"", project)
        @test occursin(
            "repo-rev = \"$line_cable_models_revision\"",
            manifest
        )
        @test occursin(
            "repo-url = \"https://github.com/Electa-Git/LineCableModels.jl.git\"",
            manifest
        )
        @test occursin(
            "repo-url = \"https://github.com/Electa-Git/PowerImpedance.jl.git\"",
            manifest
        )
        @test !occursin("gitlab.kuleuven.be", project)

        diagram_source = read(
            joinpath(
                @__DIR__, "..", "presentations", "live-manual", "support",
                "PowerImpedanceDiagramExt.jl"
            ),
            String
        )
        @test occursin("module PowerImpedanceDiagramExt", diagram_source)
        @test occursin(diagram_revision, diagram_source)
        @test occursin("include(\"PowerImpedanceDiagramExt/projection.jl\")", diagram_source)
        @test occursin("include(\"PowerImpedanceDiagramExt/render.jl\")", diagram_source)
        @test !occursin("plotbuilder.jl", lowercase(diagram_source))
        @test !occursin("PowerImpedance.PlotBuilder", diagram_source)

        serve_source = read(joinpath(@__DIR__, "..", "serve.jl"), String)
        launch_source = read(joinpath(@__DIR__, "..", "launch.sh"), String)
        @test occursin("Bonito.route!(server, playground_routes())", serve_source)
        @test !occursin("start_deck_preparations!()", serve_source)
        @test occursin("Sys.which(\"xdg-open\")", serve_source)
        @test occursin("--threads=auto", launch_source)
        @test stat(joinpath(@__DIR__, "..", "launch.sh")).mode & 0o111 != 0

        home_source = read(joinpath(@__DIR__, "..", "home.jl"), String)
        @test occursin("module PlaygroundHome", home_source)
        @test occursin("Bonito.Dropdown(", home_source)
        @test occursin("TerminalOutput", home_source)
        @test occursin("playground-presentation-activate", home_source)
        @test occursin("playground-presentation-open", home_source)
        @test occursin("layout_split_bottom(", home_source)
        @test occursin("LineCableModels playground", home_source)

        function route_html(route_app)
            parent = Session(NoConnection(); asset_server = NoServer())
            io = IOBuffer()
            rendered = Bonito.show_html(io, route_app; parent)
            try
                return String(take!(io))
            finally
                close(rendered)
                close(parent)
            end
        end

        route_htmls = Dict(
            route => route_html(route_app)
        for (route, route_app) in routes.routes
        )
        for html in values(route_htmls)
            @test occursin("--lc-bg: #1f2424", html)
            @test !occursin("-theme.css", html)
            for marker in (
                "lc-live-docs",
                "lc-sidebar",
                "page-search",
                "page-navigation",
                "lc-navbar",
                "playground-app",
                "playground-pages",
                "lc-manual",
                "lc-previous-page",
                "lc-next-page",
                "lc-page-position",
                "lc-fullscreen-toggle",
                "lc-print",
                "lc-brand-logo"
            )
                @test occursin(marker, html)
            end
        end

        home_html = route_htmls["/"]
        @test count("class=\"lc-page ", home_html) == 1
        @test occursin("home-page-home", home_html)
        @test occursin("LineCableModels playground", home_html)
        @test occursin("playground-presentation-select", home_html)
        @test occursin("playground-presentation-activate", home_html)
        @test occursin("playground-presentation-open", home_html)
        @test occursin("playground-preflight-console", home_html)
        @test occursin("playground-preflight-terminal", home_html)
        @test occursin("Preflight activity", home_html)
        @test occursin("class=\"lc-preparation-progress\"", home_html)
        @test occursin("lc-bonito-markdown", home_html)
        @test occursin("data-math-renderer=\"bonito\"", home_html)
        @test !occursin("core-radius-control", home_html)
        @test !occursin("deck-parametric-ohl-ugc-transition-page-corridor", home_html)

        cable_html = route_htmls["/presentations/live-manual/core-and-insulation"]
        @test count("class=\"lc-page ", cable_html) == 1
        @test occursin("deck-core-and-insulation-page-core-and-insulation", cable_html)
        @test occursin("lc-canvas", cable_html)
        @test occursin("core-radius-control", cable_html)
        @test occursin("insulation-thickness-control", cable_html)
        @test occursin("cable-plot", cable_html)
        @test occursin("href=\"./core-and-insulation#core-and-insulation\"", cable_html)
        @test occursin("href=\"./parametric-ohl-ugc-transition#corridor\"", cable_html)
        @test occursin("href=\"./live-julia-workspace#live-julia-workspace\"", cable_html)
        @test occursin("href=\"../../\"", cable_html)
        @test !occursin("home-page-home", cable_html)
        @test !occursin("ugc-share-control", cable_html)

        transition_html = route_htmls[
        "/presentations/live-manual/parametric-ohl-ugc-transition"
]
        @test transition_setup_count[] == 1
        @test count("class=\"lc-page ", transition_html) == 3
        @test occursin("deck-parametric-ohl-ugc-transition-page-corridor", transition_html)
        @test occursin("deck-parametric-ohl-ugc-transition-page-b5-impedance", transition_html)
        @test occursin("deck-parametric-ohl-ugc-transition-page-linearization", transition_html)
        @test occursin("data-page-id=\"b5-impedance\"", transition_html)
        @test occursin("hidden=\"true\"", transition_html)
        @test occursin("Application case - parametric OHL/UGC transition", transition_html)
        @test occursin("ugc-share-control", transition_html)
        @test occursin("recache-powerflow-control", transition_html)
        @test occursin("transition-network-diagram", transition_html)
        @test occursin("b5-impedance-plot", transition_html)
        @test !occursin("home-page-home", transition_html)
        @test !occursin("core-radius-control", transition_html)

        repl_html = route_htmls["/presentations/live-manual/live-julia-workspace"]
        @test count("class=\"lc-page ", repl_html) == 1
        @test occursin("deck-live-julia-workspace-page-live-julia-workspace", repl_html)
        @test occursin("Live Julia workspace", repl_html)
        @test occursin("live-julia-run", repl_html)
        @test occursin("live-julia-interrupt", repl_html)
        @test occursin("live-julia-restart", repl_html)
        @test occursin("live-julia-clear", repl_html)
        @test occursin("live-julia-editor", repl_html)
        @test occursin("live-julia-console", repl_html)
        @test occursin("live-julia-terminal", repl_html)
        @test !occursin("core-radius-control", repl_html)

        preparing_state = preflight_state(label = "Fixture activation required")
        preparing_resource = preflight_resource(
            id = "fixture-preflight",
            title = "Fixture preflight",
            state = preparing_state,
            activate = () -> nothing
        )
        preparing_decks = map(test_decks) do deck
            deck.id == "parametric-ohl-ugc-transition" ?
            merge(deck, (; resources = (preparing_resource,))) : deck
        end
        preparing_presentation = merge(
            test_presentation,
            (;
                all_decks = preparing_decks,
                decks = preparing_decks,
                resources = presentation_resources(preparing_decks)
            )
        )
        preparing_app = playground_app(
            assets;
            presentations = [preparing_presentation],
            presentation = preparing_presentation,
            deck = find_deck("parametric-ohl-ugc-transition", preparing_decks)
        )
        preparing_html = route_html(preparing_app)
        @test occursin("lc-documenter is-preparing", preparing_html)
        @test occursin("Preparing presentation", preparing_html)
        @test occursin("lc-preparation-progress", preparing_html)
        @test !occursin("transition-network-diagram", preparing_html)

        hidden_test_decks = select_rendered_decks(map(test_decks) do deck
            deck.id == "core-and-insulation" ?
            merge(deck, (; render = false)) : deck
        end)
        hidden_presentation = merge(
            test_presentation,
            (;
                all_decks = hidden_test_decks,
                decks = hidden_test_decks,
                resources = presentation_resources(hidden_test_decks)
            )
        )
        hidden_routes = playground_routes(
            assets;
            presentations = [hidden_presentation]
        )
        @test Set(keys(hidden_routes.routes)) ==
              Set((
            "/",
            "/presentations/live-manual/parametric-ohl-ugc-transition",
            "/presentations/live-manual/live-julia-workspace"
        ))
        hidden_html = route_html(hidden_routes.routes["/"])
        @test !occursin("deck-core-and-insulation-page-core-and-insulation", hidden_html)
        @test !occursin("core-radius-control", hidden_html)
        @test occursin("home-page-home", hidden_html)
        @test !occursin("deck-parametric-ohl-ugc-transition-page-corridor", hidden_html)
    end
end
