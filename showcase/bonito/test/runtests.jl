using Test
using LineCableModels
using Bonito
using WGLMakie
using CairoMakie
using GraphMakie
using Graphs
using Markdown
using NetworkLayout
import PowerImpedance

CairoMakie.activate!()
include(joinpath(@__DIR__, "..", "app.jl"))

const CablePage = CableDesignPage
const TransitionPage = OHLUGCTransitionPage
const ReplPage = LiveJuliaPage

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
        preparation = preparation_status()

        @test layout isa Bonito.Hyperscript.Node
        @test preparation isa Bonito.Hyperscript.Node
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
        @test TransitionPage.transition_color(0.0) == TransitionPage.OHL_COLOR
        @test TransitionPage.transition_color(1.0) == TransitionPage.UGC_COLOR

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
            ==(TransitionPage.transition_color(0.60)),
            diagram.plots.graph.edge_color[][corridor_edges]
        )
        @test all(==(6.0f0), diagram.plots.graph.edge_width[][corridor_edges])
        @test_throws DomainError TransitionPage.set_corridor_lengths!(network, NaN)

        mktempdir() do directory
            output = joinpath(directory, "network.svg")
            Makie.save(output, diagram.figure)
            @test isfile(output)
            @test filesize(output) > 0
        end
    end

    @testset "deck discovery" begin
        @test length(DECK_DESCRIPTORS) == 4
        @test getproperty.(DECK_DESCRIPTORS, :id) == [
            "overview",
            "core-and-insulation",
            "parametric-ohl-ugc-transition",
            "live-julia-workspace"
        ]
        @test all(deck -> deck.render, DECK_DESCRIPTORS)
        @test allunique(deck.id for deck in DECK_DESCRIPTORS)
        @test allunique(deck.source for deck in DECK_DESCRIPTORS)
        transition_deck = find_deck("parametric-ohl-ugc-transition", DECK_DESCRIPTORS)
        @test hasproperty(transition_deck, :prepare)
        @test hasproperty(transition_deck, :ready)
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
              ["overview", "parametric-ohl-ugc-transition", "live-julia-workspace"]
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
        transition_setup = _ -> begin
            transition_setup_count[] += 1
            return nothing
        end
        corridor_stub = (_, _) -> (;
            body = (
            DOM.h1("Corridor composition"),
            DOM.div(
                DOM.div(; id = "ugc-share-control"),
                DOM.div(; id = "transition-network-diagram");
                class = "lc-interactive-grid lc-transition-grid"
            )
        )
        )
        impedance_stub = (_, _) -> (;
            body = (
            DOM.h1("Driving-point impedance at B5"),
            DOM.div(; id = "b5-impedance-plot")
        )
        )
        linearization_stub = (_, _) -> (;
            body = (
            DOM.h1("Operating point and linearization"),
            DOM.button(; id = "recache-powerflow-control")
        )
        )
        repl_stub = (_, _) -> (;
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
                    ready = () -> true
                ))
            elseif deck.id == "live-julia-workspace"
                pages = map(page -> merge(page, (; build = repl_stub)), deck.pages)
                merge(deck, (; pages))
            else
                deck
            end
        end
        app = cable_app(assets; decks = test_decks)
        routes = cable_routes(assets; decks = test_decks)
        @test app isa App
        @test isnothing(app.loading_page)
        @test Set(keys(routes.routes)) == Set((
            "/",
            "/core-and-insulation",
            "/parametric-ohl-ugc-transition",
            "/live-julia-workspace"
        ))
        @test app.indicator.size == 8
        @test app.indicator.top == "20px"
        @test app.indicator.right == "16px"
        @test assets.app_css isa String
        @test occursin(".lc-documenter", assets.app_css)
        @test length(RENDERED_DECKS) == 4
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
        @test occursin("load_deck_descriptors()", source)
        @test occursin("deck -> deck.render", source)
        @test occursin("page -> page.render", source)
        @test occursin("function cable_routes", source)
        @test occursin("deck_route(deck, index)", source)
        @test occursin("state = deck.setup(session)", source)
        @test occursin("window.addEventListener(\"hashchange\"", source)
        @test occursin("activateDeckPage", source)
        @test occursin("window.location.assign", source)
        @test occursin("root.requestFullscreen", source)
        @test occursin("window.sessionStorage", source)
        @test occursin("linecable:presentation-mode", source)
        @test occursin("setFocusMode(readPresentationMode(), false)", source)
        @test !occursin("fallbackFocus", source)
        @test occursin("window.print()", source)
        @test occursin("event.pointerType === \"touch\"", source)
        @test occursin("window.visualViewport", source)
        @test occursin("window.devicePixelRatio", source)
        @test occursin("document.documentElement.clientWidth", source)
        @test occursin("resolution:", source)
        @test occursin("window.setTimeout(applyViewport, 240)", source)
        @test occursin("start_deck_preparations!", source)
        @test !occursin("loading_page =", source)
        @test !occursin("new MutationObserver", source)
        @test !occursin("window.history.pushState", source)
        @test !occursin("ES6Module", source)
        @test !occursin("Reveal", source)
        @test !occursin("dispatchEvent(new Event(\"resize\"))", source)
        @test !occursin("wglmakie_screen", source)

        cable_source = read(
            joinpath(@__DIR__, "..", "pages", "020_core_and_insulation.jl"),
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
            joinpath(@__DIR__, "..", "pages", "030_ohl_ugc_transition.jl"),
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
        @test occursin("addprocs(1", case_source)
        @test occursin("RemoteChannel", case_source)
        @test occursin("remotecall(compute_prepared_case", case_source)
        @test occursin(":powerflow,", case_source)
        @test occursin("Solving the reference power flow", case_source)
        @test occursin("deepcopy((prepared.network, prepared.network_model))", case_source)
        @test !occursin("PowerImpedance.plot(", case_source)
        @test !occursin("Observable{NetworkState}", case_source)

        repl_source = read(
            joinpath(@__DIR__, "..", "pages", "040_live_julia.jl"),
            String
        )
        worker_source = read(joinpath(@__DIR__, "..", "repl_worker.jl"), String)
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
        @test occursin(
            "repo-url = \"https://github.com/Electa-Git/PowerImpedance.jl.git\"",
            manifest
        )
        @test !occursin("gitlab.kuleuven.be", project)

        diagram_source = read(
            joinpath(@__DIR__, "..", "PowerImpedanceDiagramExt.jl"),
            String
        )
        @test occursin("module PowerImpedanceDiagramExt", diagram_source)
        @test occursin(diagram_revision, diagram_source)
        @test occursin("include(\"PowerImpedanceDiagramExt/projection.jl\")", diagram_source)
        @test occursin("include(\"PowerImpedanceDiagramExt/render.jl\")", diagram_source)
        @test !occursin("plotbuilder.jl", lowercase(diagram_source))
        @test !occursin("PowerImpedance.PlotBuilder", diagram_source)

        serve_source = read(joinpath(@__DIR__, "..", "serve.jl"), String)
        @test occursin("Bonito.route!(server, cable_routes())", serve_source)
        @test occursin("start_deck_preparations!()", serve_source)
        @test occursin("Sys.which(\"xdg-open\")", serve_source)
        @test stat(joinpath(@__DIR__, "..", "launch.sh")).mode & 0o111 != 0

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
                "linecable-app",
                "linecable-pages",
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
            @test occursin("href=\"./core-and-insulation#core-and-insulation\"", html)
            @test occursin("href=\"./parametric-ohl-ugc-transition#corridor\"", html)
            @test occursin("href=\"./live-julia-workspace#live-julia-workspace\"", html)
        end

        overview_html = route_htmls["/"]
        @test count("class=\"lc-page ", overview_html) == 1
        @test occursin("deck-overview-page-overview", overview_html)
        @test occursin("lc-bonito-markdown", overview_html)
        @test occursin("data-math-renderer=\"bonito\"", overview_html)
        @test occursin("lc-preparation-progress", overview_html)
        @test !occursin("core-radius-control", overview_html)
        @test !occursin("deck-parametric-ohl-ugc-transition-page-corridor", overview_html)

        cable_html = route_htmls["/core-and-insulation"]
        @test count("class=\"lc-page ", cable_html) == 1
        @test occursin("deck-core-and-insulation-page-core-and-insulation", cable_html)
        @test occursin("lc-canvas", cable_html)
        @test occursin("core-radius-control", cable_html)
        @test occursin("insulation-thickness-control", cable_html)
        @test occursin("cable-plot", cable_html)
        @test !occursin("deck-overview-page-overview", cable_html)
        @test !occursin("ugc-share-control", cable_html)

        transition_html = route_htmls["/parametric-ohl-ugc-transition"]
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
        @test !occursin("deck-overview-page-overview", transition_html)
        @test !occursin("core-radius-control", transition_html)

        repl_html = route_htmls["/live-julia-workspace"]
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

        preparing_decks = map(test_decks) do deck
            deck.id == "parametric-ohl-ugc-transition" ?
            merge(deck, (; ready = () -> false)) : deck
        end
        preparing_app = cable_app(
            assets;
            decks = preparing_decks,
            deck = find_deck("parametric-ohl-ugc-transition", preparing_decks)
        )
        preparing_html = route_html(preparing_app)
        @test occursin("lc-documenter is-preparing", preparing_html)
        @test occursin("Preparing application case", preparing_html)
        @test occursin("lc-preparation-progress", preparing_html)
        @test !occursin("transition-network-diagram", preparing_html)

        hidden_test_decks = select_rendered_decks(map(test_decks) do deck
            deck.id == "core-and-insulation" ?
            merge(deck, (; render = false)) : deck
        end)
        hidden_routes = cable_routes(assets; decks = hidden_test_decks)
        @test Set(keys(hidden_routes.routes)) ==
              Set(("/", "/parametric-ohl-ugc-transition", "/live-julia-workspace"))
        hidden_html = route_html(hidden_routes.routes["/"])
        @test !occursin("deck-core-and-insulation-page-core-and-insulation", hidden_html)
        @test !occursin("core-radius-control", hidden_html)
        @test occursin("deck-overview-page-overview", hidden_html)
        @test !occursin("deck-parametric-ohl-ugc-transition-page-corridor", hidden_html)
    end
end
