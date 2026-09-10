"""Shared explicit selection/preparation and diagnostics for both scientific consumers."""
struct StudyRuntime
    "Owned run context; no worker is allocated by construction or rendering."
    client::RuntimeClient
end

function Bonito.jsrender(session::Session, panel::StudyRuntime)
    cases = (LineParameters(), CorridorImpedance())
    node = WorkspacePage("Workers and preparation", DOM.div(
        (ViewportFrame(case_title(case), DOM.div(
            DOM.p(preparation_note(case); class="lc-study-note"),
            WorkerSelector(panel.client, role(case); profiles=(profile(case),)),
            PreparationStatus(panel.client, role(case); parameters=preparation_inputs(case));
            class="lc-content-stack lc-panel-content"); sizing=:content) for case in cases)...,
        ViewportFrame("Private Julia terminal · optional", DOM.div(
            WorkerSelector(panel.client, :terminal; profiles=("julia-terminal",));
            class="lc-panel-content"); sizing=:content),
        Disclosure("Worker diagnostics", WorkerDiagnostics(panel.client));
        class="lc-content-stack"); eyebrow="RUNTIME",
        description="Assign each role, then prepare its executor. Connection, preparation and execution are separate states.")
    root = DOM.section(node; class="lc-study-runtime")
    return Bonito.jsrender(session, DOM.div(styles(), ComponentXRay.instrument(session, root, panel); style="display: contents;"))
end

function ComponentXRay.inspection(panel::StudyRuntime)
    return ComponentXRay.ComponentInspection(panel; name="StudyRuntime",
        source=ComponentXRay.source_reference(@__MODULE__, @__FILE__, @__LINE__),
        css_scopes=[".lc-study-note"],
        notes=["Composes the same role controls and diagnostics as the scientific views; allocates nothing on entry."])
end
