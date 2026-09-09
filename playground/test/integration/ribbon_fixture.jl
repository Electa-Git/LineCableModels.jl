# Browser-only fixtures: real components, no gallery stylesheet or backend.
using LineCableModelsPlayground, Bonito
const LCM = LineCableModelsPlayground
const WB = LCM.WorkbenchUI

function controls()
    return (
        ToolbarButton(:normal; icon=:cable, label="Cable"),
        ToolbarButton(:bus; icon=:busbar, label="Bus", active=true),
        ToolbarButton(:busy; icon=:chart, label="Compiling", busy=true),
        ToolbarButton(:disabled; icon=:play, label="Remote", disabled=true),
        ToolbarToggle(:snap; icon=:grid, label="Snap", checked=true),
        ToolbarDropdown(:layer, [:nominal => "Nominal", :results => "Results"];
            icon=:layers, label="Layer"),
        ToolbarNumber(:length; label="Length", value=12, minimum=0, maximum=100),
    )
end

function specimen()
    events = Observable(0)
    binding = ToolbarBinding(_ -> (events[] += 1); namespace=:ribbon_theme)
    tabs = map((:small, :medium, :large)) do size
        RibbonTab(size, string(size), RibbonGroup("Controls", controls()...; size))
    end
    ribbon = Ribbon(tabs...; binding,
        quick_access=(ToolbarButton(:quick; icon=:save, tooltip="Save"),))
    return DOM.div(
        ribbon,
        LCM.toolbar(collect(controls()); binding),
        DOM.output(events; id="fixture-events");
        id="fixture-components"
    )
end

struct RibbonWorkbench <: WB.AbstractWorkbench end
WB.initialize(::RibbonWorkbench, session) = Observable(:ribbon)
WB.compose(::RibbonWorkbench, active) = WB.Workbench(
    namespace=:ribbon_theme,
    identity=WB.Identity("Ribbon theme audit", "Test fixture"),
    navigation=WB.Sidebar(WB.NavGroup("Fixture", WB.NavItem(:ribbon, "Ribbon")); active),
    workspace=WB.ViewStack(WB.View(:ribbon, "Ribbon", specimen()); active),
)
WB.handle!(::RibbonWorkbench, state, action) = nothing

routes = Bonito.HTTPServer.Routes()
LCM.register_static_site_routes!(routes)
LCM.register_widget_routes!(routes)
Bonito.route!(routes, "/workbenches/template" => LCM.TemplateWorkbench.app(; xray=true))
Bonito.route!(routes, "/fixture/standalone" => App() do session
    DOM.div(DOM.style(LCM.BRAND_THEME), DOM.style(LCM.CONTROL_CONTRACT),
        LCM.widget_theme_script(), specimen())
end)
Bonito.route!(routes, "/fixture/workbench" => WB.workbench_app(RibbonWorkbench()))
Bonito.route!(routes, "/fixture/embedded" => App() do session
    DOM.div(DOM.style(LCM.BRAND_THEME), LCM.widget_theme_script(),
        DOM.iframe(; src="/widgets/ribbon", width="1200", height="850", id="ribbon-frame"))
end)
server = Bonito.Server("127.0.0.1", parse(Int, ARGS[1]); routes)
try
    wait(Condition())
finally
    close(server)
end
