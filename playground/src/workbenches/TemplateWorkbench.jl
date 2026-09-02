module TemplateWorkbench

using Bonito
using ..WorkbenchUI

export Application, app

const TEMPLATE_STYLES_PATH = joinpath(@__DIR__, "template_workbench.css")
include_dependency(TEMPLATE_STYLES_PATH)
const TEMPLATE_STYLES = read(TEMPLATE_STYLES_PATH, String)
const PLAYGROUND_LOGO = Bonito.Asset(
    normpath(joinpath(@__DIR__, "..", "..", "assets", "logo.svg"))
)

struct Application <: AbstractWorkbench end

struct State
    active_view::Observable{Symbol}
    active_dock::Observable{Symbol}
    layer::Observable{Symbol}
    depth::Observable{Float64}
    spacing::Observable{Float64}
    scenario::Observable{String}
    earth_model::Observable{Symbol}
    loading::Observable{Bool}
    status::Observable{String}
    messages::Observable{Vector{NamedTuple{(:channel, :text),Tuple{Symbol,String}}}}
end

struct SelectView <: AbstractWorkbenchAction
    id::Symbol
end

struct SetLayer <: AbstractWorkbenchAction
    id::Symbol
end

struct FitView <: AbstractWorkbenchAction end
struct ResetView <: AbstractWorkbenchAction end
struct NewSession <: AbstractWorkbenchAction end
struct AppendOutput <: AbstractWorkbenchAction end
struct ClearOutput <: AbstractWorkbenchAction end
struct ToggleLoading <: AbstractWorkbenchAction end

struct OverviewView <: AbstractWorkbenchView end

struct GeometryViewport <: AbstractWorkbenchView
    state::State
end

struct GeometryScene <: AbstractWorkbenchView
    state::State
end

struct GeometryInputs <: AbstractWorkbenchView
    state::State
end

struct SweepViewport <: AbstractWorkbenchView
    state::State
end

struct EmptyResultsView <: AbstractWorkbenchView end

struct PropertiesView <: AbstractWorkbenchView
    state::State
end

struct OutputView <: AbstractWorkbenchView
    state::State
end

struct JobActivityView <: AbstractWorkbenchView
    state::State
end

struct DiagnosticsView <: AbstractWorkbenchView end
struct BrokerActivityView <: AbstractWorkbenchView end

function append_message!(state::State, channel::Symbol, text::AbstractString)
    state.messages[] = [state.messages[]; (; channel, text=string(text))]
    return nothing
end

function WorkbenchUI.initialize(::Application, session)
    messages = [
        (; channel=:info, text="Template workbench session initialized"),
        (; channel=:output, text="All viewports mounted · no numerical work dispatched"),
        (; channel=:broker, text="Broker capability not installed for this template"),
    ]
    return State(
        Observable(:geometry),
        Observable(:output),
        Observable(:nominal),
        Observable(1.0),
        Observable(6.0),
        Observable("Baseline corridor"),
        Observable(:default),
        Observable(false),
        Observable("ready · local template"),
        Observable(messages)
    )
end

function WorkbenchUI.compose(::Application, state::State)
    identity = Identity(
        "LineCableModels.jl",
        "Engineering workbench";
        mark=DOM.img(src=PLAYGROUND_LOGO, alt="")
    )

    menu = MenuBar(
        Menu(
            "File",
            Command(NewSession(); icon=:document, label="New mock session"),
            Command(nothing; icon=:download, label="Export", disabled=true)
        ),
        Menu(
            "View",
            Command(SelectView(:overview); icon=:workbench, label="Overview"),
            Command(SelectView(:geometry); icon=:geometry, label="Geometry"),
            Command(SelectView(:frequency); icon=:chart, label="Frequency sweep")
        ),
        Menu(
            "Help",
            Command(AppendOutput(); icon=:document, label="Record template event")
        )
    )

    toolbar = Toolbar(
        Command(FitView(); icon=:fit, tooltip="Fit geometry to viewport"),
        Command(ResetView(); icon=:reset, label="Reset view"),
        CommandChoice(
            "Layer",
            ChoiceOption(:nominal, "Nominal", SetLayer(:nominal)),
            ChoiceOption(:uncertainty, "Uncertainty", SetLayer(:uncertainty)),
            ChoiceOption(:interfaces, "Interfaces", SetLayer(:interfaces));
            selected=state.layer,
            icon=:layers,
            tooltip="Displayed geometry layer"
        ),
        Command(ToggleLoading(); icon=:play, label="Mock activity"),
        Command(nothing; icon=:download, label="Export", disabled=true);
        label="Geometry commands",
        size=:medium
    )

    navigation = Sidebar(
        NavGroup(
            "Workspace",
            NavItem(
                :overview,
                "System overview",
                SelectView(:overview);
                icon=:workbench
            ),
            NavItem(
                :geometry,
                "Cable geometry",
                SelectView(:geometry);
                icon=:geometry
            ),
            NavItem(
                :frequency,
                "Frequency sweep",
                SelectView(:frequency);
                icon=:chart
            )
        ),
        NavGroup(
            "Analysis",
            NavItem(
                :results,
                "Result archive",
                SelectView(:results);
                icon=:archive
            ),
            NavItem(
                :remote,
                "Remote execution",
                nothing;
                icon=:cloud,
                tooltip="Remote execution unavailable in the template",
                disabled=true
            )
        );
        active=state.active_view,
        footer=DOM.div(
            DOM.span(; class="lc-wb-demo-status-dot"),
            DOM.span("local template · no services"; class="lc-wb-demo-footer-copy");
            class="lc-wb-demo-footer"
        )
    )

    workspace = ViewStack(
        View(:overview, "System overview", OverviewView()),
        View(:geometry, "Cable geometry", GeometryViewport(state)),
        View(:frequency, "Frequency sweep", SweepViewport(state)),
        View(:results, "Result archive", EmptyResultsView());
        active=state.active_view
    )

    inspector = Inspector(
        "Cable section",
        PropertiesView(state);
        subtitle="Inspector"
    )

    output = Dock(
        DockTab(:output, "Output", OutputView(state)),
        DockTab(:jobs, "Jobs", JobActivityView(state)),
        DockTab(:diagnostics, "Diagnostics", DiagnosticsView()),
        DockTab(:broker, "Broker", BrokerActivityView());
        active=state.active_dock
    )

    status = StatusBar(
        DOM.span(state.status),
        right=DOM.span("template_workbench · session-local")
    )

    return Workbench(;
        namespace=:template_workbench,
        identity=identity,
        menu=menu,
        toolbar=toolbar,
        navigation=navigation,
        workspace=workspace,
        inspector=inspector,
        output=output,
        status=status
    )
end

function WorkbenchUI.handle!(::Application, state::State, action::SelectView)
    state.active_view[] = action.id
    state.status[] = "view · $(replace(string(action.id), '_' => ' '))"
    return nothing
end

function WorkbenchUI.handle!(::Application, state::State, action::SetLayer)
    state.layer[] = action.id
    state.status[] = "geometry layer · $(action.id)"
    append_message!(state, :output, "Display layer changed to $(action.id)")
    return nothing
end

function WorkbenchUI.handle!(::Application, state::State, ::FitView)
    state.status[] = "geometry viewport · fitted"
    append_message!(state, :output, "Viewport fitted to declared scene bounds")
    return nothing
end

function WorkbenchUI.handle!(::Application, state::State, ::ResetView)
    state.layer[] = :nominal
    state.depth[] = 1.0
    state.spacing[] = 6.0
    state.scenario[] = "Baseline corridor"
    state.earth_model[] = :default
    state.loading[] = false
    state.status[] = "geometry viewport · reset"
    append_message!(state, :output, "View and local selectors returned to defaults")
    return nothing
end

function WorkbenchUI.handle!(::Application, state::State, ::NewSession)
    state.active_view[] = :geometry
    state.active_dock[] = :output
    state.layer[] = :nominal
    state.depth[] = 1.0
    state.spacing[] = 6.0
    state.scenario[] = "Baseline corridor"
    state.earth_model[] = :default
    state.loading[] = false
    state.messages[] = [
        (; channel=:info, text="Mock session reinitialized"),
    ]
    state.status[] = "ready · local template"
    return nothing
end

function WorkbenchUI.handle!(::Application, state::State, ::AppendOutput)
    append_message!(state, :info, "Explicit template event recorded")
    state.active_dock[] = :output
    state.status[] = "output · event appended"
    return nothing
end

function WorkbenchUI.handle!(::Application, state::State, ::ClearOutput)
    state.messages[] = NamedTuple{(:channel, :text),Tuple{Symbol,String}}[]
    state.status[] = "output · cleared"
    return nothing
end

function WorkbenchUI.handle!(::Application, state::State, ::ToggleLoading)
    state.loading[] = !state.loading[]
    state.active_dock[] = :jobs
    state.status[] = state.loading[] ? "mock job · running" : "mock job · idle"
    return nothing
end

function Bonito.jsrender(session::Session, ::OverviewView)
    return Bonito.jsrender(
        session,
        DOM.div(
            DOM.style(TEMPLATE_STYLES),
            DOM.header(
                DOM.span("WORKBENCH TEMPLATE"; class="lc-wb-demo-kicker"),
                DOM.h1("Browser-hosted engineering workbench"),
                DOM.p(
                    "One persistent shell coordinates views, properties, commands, " *
                    "operational output, and status without executing scientific code."
                )
            ),
            DOM.div(
                DOM.article(
                    DOM.span("01"),
                    DOM.h2("Persistent views"),
                    DOM.p("Workspace and dock contents remain mounted while selection changes.")
                ),
                DOM.article(
                    DOM.span("02"),
                    DOM.h2("Typed commands"),
                    DOM.p("Local action values reach Julia through ordinary multiple dispatch.")
                ),
                DOM.article(
                    DOM.span("03"),
                    DOM.h2("Explicit boundary"),
                    DOM.p("A real application may dispatch jobs; this template has no broker client.")
                );
                class="lc-wb-demo-principles"
            );
            class="lc-wb-demo-overview"
        )
    )
end

function Bonito.jsrender(session::Session, view::GeometryViewport)
    depth_label = map(session, view.state.depth) do depth
        "$(round(depth; digits=2)) m burial depth"
    end
    layout = SplitPane(
        GeometryScene(view.state),
        GeometryInputs(view.state);
        orientation=:horizontal,
        ratio=0.72,
        min_first="24rem",
        min_second="15rem"
    )
    return Bonito.jsrender(
        session,
        DOM.div(
            DOM.header(
                DOM.div(
                    DOM.span("SCENE VIEWPORT"; class="lc-wb-demo-kicker"),
                    DOM.h1("Cable geometry")
                ),
                DOM.span(depth_label; class="lc-wb-demo-machine")
            ),
            layout;
            class="lc-wb-demo-viewport"
        )
    )
end

function Bonito.jsrender(session::Session, view::GeometryScene)
    scene_class = map(session, view.state.layer) do layer
        "lc-wb-demo-scene layer-$(layer)"
    end
    scene_attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => "Mock cross-section of two underground cables",
    )
    centered_text = Dict{Symbol,Any}(Symbol("text-anchor") => "middle")
    spacing_label = SVG.text(
        "$(round(view.state.spacing[]; digits=2)) m";
        centered_text...,
        x="400",
        y="116",
        class="label"
    )
    synchronize = js"""
    (() => {
        const label = $(spacing_label);
        $(view.state.spacing).on(value => {
            label.textContent = `${Number(value).toFixed(2)} m`;
        });
    })();
    """
    return Bonito.jsrender(
        session,
        DOM.figure(
            DOM.figcaption(
                DOM.strong("Two-conductor line arrangement"),
                DOM.span("local mock geometry · no simulation")
            ),
            SVG.svg(
                SVG.line(; x1="70", y1="88", x2="730", y2="88", class="ground"),
                SVG.line(; x1="250", y1="88", x2="250", y2="218", class="dimension"),
                SVG.line(; x1="550", y1="88", x2="550", y2="218", class="dimension"),
                SVG.line(; x1="250", y1="126", x2="550", y2="126", class="dimension"),
                spacing_label,
                SVG.g(
                    SVG.circle(; cx="250", cy="218", r="65", class="envelope"),
                    SVG.circle(; cx="250", cy="218", r="54", class="jacket"),
                    SVG.circle(; cx="250", cy="218", r="34", class="insulation"),
                    SVG.circle(; cx="250", cy="218", r="18", class="conductor"),
                    SVG.circle(; cx="250", cy="218", r="34", class="interface"),
                    SVG.circle(; cx="250", cy="218", r="18", class="interface")
                ),
                SVG.g(
                    SVG.circle(; cx="550", cy="218", r="65", class="envelope"),
                    SVG.circle(; cx="550", cy="218", r="54", class="jacket"),
                    SVG.circle(; cx="550", cy="218", r="34", class="insulation"),
                    SVG.circle(; cx="550", cy="218", r="18", class="conductor"),
                    SVG.circle(; cx="550", cy="218", r="34", class="interface"),
                    SVG.circle(; cx="550", cy="218", r="18", class="interface")
                ),
                SVG.text("CABLE 01"; centered_text..., x="250", y="308", class="node"),
                SVG.text("CABLE 02"; centered_text..., x="550", y="308", class="node");
                scene_attributes...,
                viewBox="0 0 800 340",
                role="img"
            ),
            DOM.script(synchronize);
            class=scene_class
        )
    )
end

function Bonito.jsrender(session::Session, view::GeometryInputs)
    scenario_changed = Observable(view.state.scenario[])
    spacing_changed = Observable(string(view.state.spacing[]))
    earth_changed = Observable(string(view.state.earth_model[]))

    on(session, scenario_changed) do value
        view.state.scenario[] = string(value)
        return nothing
    end
    on(session, spacing_changed) do value
        parsed = tryparse(Float64, value)
        isnothing(parsed) || (view.state.spacing[] = clamp(parsed, 1.0, 10.0))
        return nothing
    end
    on(session, earth_changed) do value
        view.state.earth_model[] = Symbol(value)
        return nothing
    end

    scenario_attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => "Scenario name",
    )
    spacing_attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => "Lateral spacing in metres",
    )
    earth_attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => "Earth-return model",
    )
    scenario = DOM.input(
        ;
        scenario_attributes...,
        type="text",
        value=view.state.scenario[],
        oninput=js"event => $(scenario_changed).notify(event.currentTarget.value)"
    )
    spacing = DOM.input(
        ;
        spacing_attributes...,
        type="number",
        min="1",
        max="10",
        step="0.25",
        value=string(view.state.spacing[]),
        oninput=js"event => $(spacing_changed).notify(event.currentTarget.value)"
    )
    earth = DOM.select(
        DOM.option("Default earth"; value="default", selected=true),
        DOM.option("Carson"; value="carson"),
        DOM.option("Wedepohl"; value="wedepohl");
        earth_attributes...,
        onchange=js"event => $(earth_changed).notify(event.currentTarget.value)"
    )
    synchronize = js"""
    (() => {
        const scenario = $(scenario);
        const spacing = $(spacing);
        const earth = $(earth);
        $(view.state.scenario).on(value => { scenario.value = String(value); });
        $(view.state.spacing).on(value => { spacing.value = String(value); });
        $(view.state.earth_model).on(value => { earth.value = String(value); });
    })();
    """
    spacing_readout = map(session, view.state.spacing) do value
        "$(round(value; digits=2)) m"
    end
    return Bonito.jsrender(
        session,
        DOM.section(
            DOM.header(
                DOM.span("DATA ENTRY"; class="lc-wb-demo-kicker"),
                DOM.h2("Arrangement")
            ),
            DOM.div(
                DOM.label("Scenario", scenario),
                DOM.label(
                    DOM.span("Lateral spacing"),
                    DOM.div(spacing, DOM.span("m"); class="lc-wb-demo-input-with-unit")
                ),
                DOM.label("Earth model", earth);
                class="lc-wb-demo-fields"
            ),
            DOM.dl(
                DOM.div(DOM.dt("Spacing"), DOM.dd(spacing_readout)),
                DOM.div(DOM.dt("Execution"), DOM.dd("not dispatched"; class="is-accent"));
                class="lc-wb-demo-readout"
            ),
            DOM.p(
                "This local block owns presentation state only. SplitPane can be " *
                "nested to arrange additional views or input regions."
            ),
            DOM.script(synchronize);
            class="lc-wb-demo-input-panel"
        )
    )
end

function Bonito.jsrender(session::Session, view::SweepViewport)
    state_class = map(session, view.state.loading) do loading
        loading ? "lc-wb-demo-sweep is-loading" : "lc-wb-demo-sweep"
    end
    state_copy = map(session, view.state.loading) do loading
        loading ? "Mock job running" : "No job dispatched"
    end
    return Bonito.jsrender(
        session,
        DOM.div(
            DOM.header(
                DOM.span("PERSISTENT VIEW"; class="lc-wb-demo-kicker"),
                DOM.h1("Frequency sweep"),
                DOM.p("This inert scene demonstrates loading and idle presentation states.")
            ),
            DOM.div(
                DOM.div(; class="lc-wb-demo-spinner"),
                DOM.span(state_copy),
                DOM.div(
                    (DOM.i(; style="--sample: $(index)") for index in 1:24)...;
                    class="lc-wb-demo-trace"
                );
                class=state_class
            );
            class="lc-wb-demo-generic-view"
        )
    )
end

function Bonito.jsrender(session::Session, ::EmptyResultsView)
    return Bonito.jsrender(
        session,
        DOM.div(
            icon(:archive; class="lc-wb-demo-empty-icon"),
            DOM.h1("No archived results"),
            DOM.p("Empty states are finite and explicit. Nothing is silently loading.");
            class="lc-wb-demo-empty"
        )
    )
end

function Bonito.jsrender(session::Session, view::PropertiesView)
    input_attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => "Mock burial depth",
    )
    input = DOM.input(
        ;
        input_attributes...,
        type="range",
        min="0.5",
        max="3.0",
        step="0.1",
        value=string(view.state.depth[]),
        class="lc-wb-demo-range"
    )
    changed = Observable(string(view.state.depth[]))
    on(session, changed) do value
        parsed = tryparse(Float64, value)
        isnothing(parsed) || (view.state.depth[] = parsed)
        return nothing
    end
    sync = js"""
    (() => {
        const input = $(input);
        const depth = $(view.state.depth);
        depth.on(value => { input.value = String(value); });
    })();
    """
    depth_value = map(session, view.state.depth) do depth
        "$(round(depth; digits=1)) m"
    end
    return Bonito.jsrender(
        session,
        DOM.div(
            DOM.div(
                DOM.label("Burial depth"),
                DOM.output(depth_value),
                input,
                DOM.script(sync);
                class="lc-wb-demo-property"
            ),
            DOM.div(
                DOM.span("Construction"),
                DOM.strong("Single-core cable")
            ),
            DOM.div(
                DOM.span("Spacing"),
                DOM.strong(map(session, view.state.spacing) do spacing
                    "$(round(spacing; digits=2)) m"
                end)
            ),
            DOM.div(
                DOM.span("Execution"),
                DOM.strong("not dispatched"; class="is-accent")
            );
            class="lc-wb-demo-properties"
        )
    )
end

function Bonito.jsrender(session::Session, view::OutputView)
    lines = map(session, view.state.messages) do messages
        isempty(messages) && return DOM.div(
            "No output messages";
            class="lc-wb-demo-log-empty"
        )
        return DOM.div(
            (
                DOM.div(
                    DOM.span(uppercase(string(message.channel))),
                    DOM.code(message.text);
                    class="channel-$(message.channel)"
                ) for message in messages
            )...;
            class="lc-wb-demo-log-lines"
        )
    end
    clear = Observable(false)
    on(session, clear) do _
        WorkbenchUI.handle!(Application(), view.state, ClearOutput())
        return nothing
    end
    return Bonito.jsrender(
        session,
        DOM.div(
            DOM.button(
                "Clear";
                type="button",
                onclick=js"event => $(clear).notify(true)"
            ),
            lines;
            class="lc-wb-demo-output"
        )
    )
end

function Bonito.jsrender(session::Session, view::JobActivityView)
    class = map(session, view.state.loading) do loading
        loading ? "lc-wb-demo-job is-running" : "lc-wb-demo-job"
    end
    label = map(session, view.state.loading) do loading
        loading ? "RUNNING" : "IDLE"
    end
    return Bonito.jsrender(
        session,
        DOM.div(
            DOM.div(
                DOM.span(label),
                DOM.strong("Mock frequency sweep"),
                DOM.small("Template-only lifecycle · no broker request")
            ),
            DOM.progress(value="42", max="100");
            class
        )
    )
end

function Bonito.jsrender(session::Session, ::DiagnosticsView)
    return Bonito.jsrender(
        session,
        DOM.div(
            icon(:warning),
            DOM.span("No diagnostics"),
            DOM.small("The template contains no domain validation errors.");
            class="lc-wb-demo-dock-empty"
        )
    )
end

function Bonito.jsrender(session::Session, ::BrokerActivityView)
    return Bonito.jsrender(
        session,
        DOM.div(
            DOM.span("OFFLINE"),
            DOM.strong("Broker capability not installed"),
            DOM.small("The workbench shell remains fully usable without NATS.");
            class="lc-wb-demo-broker"
        )
    )
end

app() = workbench_app(
    Application();
    title="Workbench template · LineCableModels playground"
)

end
