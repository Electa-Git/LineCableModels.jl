"""Bonito frame factories for the registered Quarto scientific deck."""
module Showcase

using Bonito
using ..ScientificViews
using ..ScientificViews.StudyCases
using ..LineCableModelsPlayground: RuntimeClient, JuliaTerminal, widget_shell

export routes

# Factories run in the owning Bonito session. No global Observables or worker
# allocations are shared between viewers, frames, or application runs.
frame(factory, title) = App(; title) do session
    widget_shell("SCIENTIFIC SHOWCASE", title, factory(session); header=false)
end

"""Return approved live routes sharing one run context; rendering never prepares work."""
routes(client::RuntimeClient) = (
    "/science/runtime" => frame(_ -> StudyRuntime(client), "Workers and preparation"),
    "/science/geometry" => frame(_ -> CableGeometry(), "Cable construction"),
    "/science/line-parameters" => frame(session -> ScientificView(session, LineParameters(), client), "Line parameters"),
    "/science/corridor" => frame(session -> ScientificView(session, CorridorImpedance(), client), "OHL / UGC case"),
    "/science/terminal" => frame(_ -> JuliaTerminal(client, :terminal; rows=18), "Private Julia terminal"),
)

end
