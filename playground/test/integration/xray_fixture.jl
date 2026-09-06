using LineCableModelsPlayground, Bonito
const LCM = LineCableModelsPlayground
const WB = LCM.WorkbenchUI
const X = LCM.ComponentXRay
const FIXTURE_CSS = read(joinpath(@__DIR__, "xray_fixture.css"), String)

struct PreviewCard
    id::String
    counter::Observable{Int}
end
X.css_editors(::PreviewCard) = Dict(
    "gap" => X.CssEditor(:length; units=["px", "rem"], minimum=0, maximum=60, step=2))
X.inspection(card::PreviewCard) = X.ComponentInspection(card;
    source=X.source_reference(@__MODULE__, @__FILE__, @__LINE__),
    parameters=[X.PropertyInspection(:name, card.id), X.PropertyInspection(:enabled, true)],
    bindings=[X.BindingInspection(:counter, card.counter)],
    actions=[X.ActionInspection(:click, :increment, identity)],
    css_scopes=[".xp-card"])
function Bonito.jsrender(session::Bonito.Session, card::PreviewCard)
    count = card.counter
    node = DOM.section(DOM.h2(card.id), DOM.output(count),
        DOM.button("Increment"; type="button", onclick=js"event => $(count).notify($(count).value + 1)"),
        DOM.input(; type="text", value="Domain input", class="lc-control-input");
        id=card.id, class="xp-card", style="--application-value: 42;")
    rendered = Bonito.jsrender(session, X.instrument(session, node, card))
    Bonito.evaljs(session, js"window.__xrayFixtureReady = true")
    return rendered
end
struct PreviewGroup{T}
    id::String
    content::T
end
X.inspection(group::PreviewGroup) = X.ComponentInspection(group;
    parameters=[X.PropertyInspection(:name, group.id)], css_scopes=[".xp-group"])
function Bonito.jsrender(session::Bonito.Session, group::PreviewGroup)
    node = DOM.div(group.content; id=group.id, class="xp-group")
    return Bonito.jsrender(session, X.instrument(session, node, group))
end
struct PreviewWorkbench <: WB.AbstractWorkbench
    tree::Bool
end
PreviewWorkbench() = PreviewWorkbench(false)
preview_contents(app::PreviewWorkbench, state) = app.tree ? (
    PreviewGroup("group-a", DOM.div(PreviewCard("card-a", state.counter),
        PreviewGroup("group-inner", PreviewCard("card-b", state.counter)))),
    PreviewCard("card-c", state.counter),
) : (PreviewCard("card-a", state.counter), PreviewCard("card-b", state.counter))
WB.initialize(::PreviewWorkbench, session) = (active=Observable(:preview), counter=Observable(0), dock=Observable(:output))
WB.handle!(::PreviewWorkbench, state, action) = nothing
WB.compose(app::PreviewWorkbench, state) = WB.Workbench(
    namespace=:xray_fixture, identity=WB.Identity("X-ray preview fixture", "No backend"),
    navigation=WB.Sidebar(WB.NavGroup("Test", WB.NavItem(:preview, "Preview")); active=state.active),
    workspace=WB.ViewStack(WB.View(:preview, "Preview", DOM.div(
        DOM.style(FIXTURE_CSS; var"data-lcm-css-source"="test/integration/xray_fixture.css"),
        DOM.style(FIXTURE_CSS; var"data-lcm-css-source"="test/integration/xray_fixture.css"),
        preview_contents(app, state)...;
        class="xp-fixture")); active=state.active),
    output=WB.Dock(WB.DockTab(:output, "Output", DOM.p("CSS preview fixture · no execution")); active=state.dock))

server = Bonito.Server("127.0.0.1", parse(Int, ARGS[1]))
LCM.register_static_site_routes!(server)
Bonito.route!(server, "/fixture/xray" => WB.workbench_app(PreviewWorkbench(); xray=X.XRayPolicy(true)))
Bonito.route!(server, "/fixture/tree" => WB.workbench_app(PreviewWorkbench(true); xray=X.XRayPolicy(true)))
Bonito.route!(server, "/fixture/readonly" => WB.workbench_app(PreviewWorkbench();
    xray=X.XRayPolicy(; permitted=true, enabled=true, css_preview=false)))
Bonito.route!(server, "/workbenches/template" => LCM.TemplateWorkbench.app(; xray=true))
try
    wait(Condition())
finally
    close(server)
end
