"""Concrete scientific workbench; reuses the showcase's views and runtime controls."""
module CableStudy

using Bonito
using ..WorkbenchUI
using ..ScientificViews
using ..ScientificViews.StudyCases
using ..ComponentXRay
using ..LineCableModelsPlayground: RuntimeClient, JuliaTerminal, WorkerDiagnostics

export Application, app

"""Bind the scientific workbench to an explicitly provided owned application run."""
struct Application <: AbstractWorkbench
    "Same-origin runtime context; does not allocate resources."
    client::RuntimeClient
end

struct SelectView <: AbstractWorkbenchAction
    id::Symbol
end

function WorkbenchUI.initialize(application::Application, session)
    client = application.client
    return (active=Observable(:runtime), dock=Observable(:diagnostics), views=(
        runtime=StudyRuntime(client), geometry=CableGeometry(),
        parameters=ScientificView(session, LineParameters(), client),
        corridor=ScientificView(session, CorridorImpedance(), client),
        terminal=JuliaTerminal(client, :terminal; rows=24)))
end

function WorkbenchUI.compose(application::Application, state)
    ids = (:runtime, :geometry, :parameters, :corridor, :terminal)
    labels = ("Workers and preparation", "Cable construction", "Line parameters", "OHL / UGC case", "Julia terminal")
    icons = (:cloud, :geometry, :chart, :chart, :terminal)
    navigation = Sidebar(NavGroup("Cable study",
        (NavItem(id, label, SelectView(id); icon) for (id,label,icon) in zip(ids,labels,icons))...);
        active=state.active, footer=DOM.a("Playground home"; href="/", target="_top"))
    workspace = ViewStack((View(id, label, getproperty(state.views, id)) for (id,label) in zip(ids,labels))...;
        active=state.active)
    return Workbench(; namespace=:cable_study,
        identity=Identity("LineCableModels", "CableStudy · scientific application"),
        navigation, workspace,
        toolbar=Toolbar(Command(SelectView(:runtime); icon=:cloud, label="Workers"),
            Command(SelectView(:parameters); icon=:chart, label="Line parameters"),
            Command(SelectView(:corridor); icon=:chart, label="OHL / UGC"),
            Command(SelectView(:terminal); icon=:terminal, label="Julia terminal"); label="Study views"),
        output=Dock(DockTab(:diagnostics, "Diagnostics", WorkerDiagnostics(application.client)); active=state.dock),
        status=StatusBar(DOM.span("UI ready · calculations require explicit assignment, preparation and Run")))
end

function WorkbenchUI.handle!(::Application, state, action::SelectView)
    hasproperty(state.views, action.id) || throw(ArgumentError("unknown CableStudy view"))
    state.active[] = action.id
    return nothing
end

"""Create an isolated-session workbench using the supplied run; no implicit worker allocation."""
app(client::RuntimeClient; xray::Bool=false) = workbench_app(Application(client);
    title="CableStudy · LineCableModels", xray=ComponentXRay.XRayPolicy(xray))

end
