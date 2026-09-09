# This check must precede package loading: setpriv's parent-death signal and
# this launch-race guard together bind the process to its actual supervisor.
ccall(:getppid, Cint, ()) == parse(Int, ENV["LCM_UI_PARENT_PID"]) || exit(70)
using LineCableModelsPlayground
const LCM = LineCableModelsPlayground
length(ARGS) == 2 || error("Expected approved UI driver and diagnostic flag")
driver = Symbol(ARGS[1])
ARGS[2] in ("true", "false") || error("Invalid diagnostic flag")
xray = ARGS[2] == "true"
routes = if driver == :template
    ("/workbenches/template" => LCM.TemplateWorkbench.app(; xray),)
elseif driver == :showcase
    LCM.Showcase.routes(LCM.RuntimeClient(LCM.UUID(ENV["LCM_RUN_ID"])))
elseif driver == :cable_study
    ("/workbenches/cable-study" => LCM.CableStudy.app(LCM.RuntimeClient(LCM.UUID(ENV["LCM_RUN_ID"])); xray),)
elseif driver == :starter
    ("/widgets/control-panel" => LCM.control_panel_widget(),)
elseif driver == :specimen
    ("/presentations/probe" => LCM.presentation_probe_widget(),
     "/widgets/control-panel" => LCM.control_panel_widget())
elseif driver == :gallery
    # Consume the exact gallery factories, not copies of their markup or CSS.
    Tuple(route => factory() for (route, factory) in LCM.WIDGET_ROUTES)
else
    error("UI driver is not installed")
end
LCM.serve_owned_ui(routes)
