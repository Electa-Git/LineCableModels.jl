"""Shared explicit selection/preparation and diagnostics for both scientific consumers."""
struct StudyRuntime
    "Owned run context; no worker is allocated by construction or rendering."
    client::RuntimeClient
end

function Bonito.jsrender(session::Session, panel::StudyRuntime)
    cases = (LineParameters(), CorridorImpedance())
    node = DOM.section(DOM.h2("Prepare the application"),
        DOM.p("Assign the scientific roles independently, then prepare each before presenting. Worker connection, preparation and execution are distinct states."; class="lc-study-note"),
        (ViewportFrame(case_title(case), DOM.div(
            DOM.p(preparation_note(case); class="lc-study-note"),
            WorkerSelector(panel.client, role(case); profiles=(profile(case),)),
            PreparationStatus(panel.client, role(case); parameters=preparation_inputs(case)))) for case in cases)...,
        ViewportFrame("Private Julia terminal · optional", WorkerSelector(panel.client, :terminal; profiles=("julia-terminal",))),
        WorkerDiagnostics(panel.client); class="lc-study-runtime")
    return Bonito.jsrender(session, DOM.div(styles(), ComponentXRay.instrument(session, node, panel); style="display: contents;"))
end

function ComponentXRay.inspection(panel::StudyRuntime)
    return ComponentXRay.ComponentInspection(panel; name="StudyRuntime",
        source=ComponentXRay.source_reference(@__MODULE__, @__FILE__, @__LINE__),
        css_scopes=[".lc-study-runtime", ".lc-study-note"],
        notes=["Composes the same role controls and diagnostics as the scientific views; allocates nothing on entry."])
end
