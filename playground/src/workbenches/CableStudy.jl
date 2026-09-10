"""Concrete scientific workbench; reuses the showcase's views and runtime controls."""
module CableStudy

using Bonito
using ..WorkbenchUI
using ..Toolkit: NavigationButton
using ..ScientificViews
using ..ScientificViews.StudyCases
using ..ComponentXRay
using ..LineCableModelsPlayground: RuntimeClient, JuliaTerminal, WorkerDiagnostics
using ..LineCableModelsPlayground: PLAYGROUND_ROOT

export Application, app

"""
    Application(client; return_button=NavigationButton("Home"; href="/", icon=icon(:home)))

Bind the scientific workbench to an owned run. `return_button` configures the
sidebar's destination independently of scientific state or worker assignment.
"""
struct Application <: AbstractWorkbench
    "Same-origin runtime context; does not allocate resources."
    client::RuntimeClient
    "Reusable navigation control supplied by the workbench owner."
    return_button::NavigationButton
end

default_return_button() = NavigationButton("Home"; href="/", icon=icon(:home))
Application(client::RuntimeClient; return_button::NavigationButton=default_return_button()) =
    Application(client, return_button)

struct SelectView <: AbstractWorkbenchAction
    id::Symbol
end

function WorkbenchUI.initialize(application::Application, session)
    client = application.client
    return (active=Observable(:runtime), dock=Observable(:diagnostics), views=(
        runtime=StudyRuntime(client; diagnostics=false), geometry=CableGeometry(),
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
        active=state.active, footer=application.return_button)
    workspace = ViewStack((View(id, label, getproperty(state.views, id)) for (id,label) in zip(ids,labels))...;
        active=state.active)
    return Workbench(; namespace=:cable_study,
        identity=Identity("LineCableModels.jl", "Cable study";
            mark=DOM.img(src=Bonito.Asset(joinpath(PLAYGROUND_ROOT, "assets", "logo.svg")), alt="")),
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

"""
    app(client; xray=false, return_button=NavigationButton("Home"; href="/", icon=icon(:home)))

Create an isolated-session workbench without implicit worker allocation. Supply
a `NavigationButton` to change the sidebar return label, destination and icon.
"""
app(client::RuntimeClient; xray::Bool=false,
    return_button::NavigationButton=default_return_button()) = workbench_app(Application(client; return_button);
    title="CableStudy · LineCableModels", xray=ComponentXRay.XRayPolicy(xray))

end
