const RUNTIME_ASSETS = Dict(
    "brand.css" => (joinpath(@__DIR__, "..", "..", "assets", "brand.css"), "text/css"),
    "control-contract.css" => (joinpath(@__DIR__, "..", "..", "assets", "control-contract.css"), "text/css"),
    "published-text.css" => (joinpath(@__DIR__, "..", "..", "assets", "published-text.css"), "text/css"),
    "forms.css" => (joinpath(@__DIR__, "..", "..", "assets", "forms.css"), "text/css"),
    "data-views.css" => (joinpath(@__DIR__, "..", "..", "assets", "data-views.css"), "text/css"),
    "workspace.css" => (joinpath(@__DIR__, "..", "..", "assets", "workspace.css"), "text/css"),
    "run.css" => (joinpath(@__DIR__, "..", "ui", "run.css"), "text/css"),
    "run.js" => (joinpath(@__DIR__, "..", "ui", "run.js"), "text/javascript"),
    "control.js" => (joinpath(@__DIR__, "..", "ui", "control.js"), "text/javascript"),
    "runtime-client.js" => (joinpath(@__DIR__, "..", "..", "assets", "runtime-client.js"), "text/javascript"),
    "runtime-controls.js" => (joinpath(@__DIR__, "..", "..", "assets", "runtime-controls.js"), "text/javascript"),
    "runtime-controls.css" => (joinpath(@__DIR__, "..", "..", "assets", "runtime-controls.css"), "text/css"),
    "runtime-terminal-client.js" => (joinpath(@__DIR__, "..", "..", "assets", "runtime-terminal-client.js"), "text/javascript"),
    "runtime-terminal.js" => (joinpath(@__DIR__, "..", "..", "assets", "runtime-terminal.js"), "text/javascript"),
    "runtime-terminal.css" => (joinpath(@__DIR__, "..", "..", "assets", "runtime-terminal.css"), "text/css"),
    "runtime-terminal.bundle.js" => (joinpath(@__DIR__, "..", "..", "assets", "vendor", "runtime-terminal.bundle.js"), "text/javascript"),
    "runtime-terminal.bundle.css" => (joinpath(@__DIR__, "..", "..", "assets", "vendor", "runtime-terminal.bundle.css"), "text/css"),
)
const RUNTIME_THEME_INIT = joinpath(@__DIR__, "..", "..", "assets", "theme-init.html")

html_text(value) = replace(string(value), '&'=>"&amp;", '<'=>"&lt;", '>'=>"&gt;", '"'=>"&quot;", '\''=>"&#39;")

function run_surface(stream, supervisor::UIHostSupervisor, run::RunRecord; status::Int=200, automatic::Bool=true)
    app = get(supervisor.registry.definitions, run.application, nothing)
    title = isnothing(app) ? "Application unavailable" : app.title
    entrypoint = isnothing(app) ? "/" : app.entrypoint
    kind = isnothing(app) ? "workbench" : String(app.kind)
    entry_surface = isnothing(app) ? "published" : String(app.entry_surface)
    busy = run.state in (:reserved, :starting)
    reason = isempty(run.reason) && busy ? "Preparing the isolated UI host…" : run.reason
    body = """<!doctype html>
    <html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
    <title>$(html_text(title)) · Runtime</title>
    $(read(RUNTIME_THEME_INIT, String))
    <link rel="stylesheet" href="/runtime/assets/brand.css">
    <link rel="stylesheet" href="/runtime/assets/control-contract.css">
    <link rel="stylesheet" href="/runtime/assets/published-text.css">
    <link rel="stylesheet" href="/runtime/assets/forms.css">
    <link rel="stylesheet" href="/runtime/assets/data-views.css">
    <link rel="stylesheet" href="/runtime/assets/workspace.css">
    <link rel="stylesheet" href="/runtime/assets/run.css"></head>
    <body class="lc-runtime-page"><main class="lc-runtime-surface" data-run="$(run.id)"
      data-application="$(html_text(run.application))" data-entrypoint="$(html_text(entrypoint))"
      data-kind="$kind" data-entry-surface="$entry_surface" data-automatic="$automatic" data-deadline="$(supervisor.limits.startup_seconds + 15)">
      <p class="lc-runtime-eyebrow">APPLICATION RUN</p><h1>$(html_text(title))</h1>
      <p id="runtime-status" class="lc-activity-status" data-busy="$busy" role="status" aria-live="polite">$(html_text(run.state))</p>
      <p id="runtime-reason">$(html_text(reason))</p>
      <p id="runtime-elapsed" class="lc-runtime-hint" aria-live="off" $(busy ? "" : "hidden")>Waiting for application readiness. First startup may take longer while Julia loads and compiles.</p>
      <p class="lc-runtime-hint">The public site remains available. Restart creates a clean run;
        unsaved UI state and terminal memory cannot be recovered after process loss.</p>
      <div class="lc-runtime-actions"><a id="runtime-open" hidden>Open application</a>
        <button class="lc-button lc-button-secondary" id="runtime-restart" type="button" hidden>Start a clean run</button>
        <button class="lc-button lc-button-secondary" id="runtime-stop" type="button">Stop run</button>
        <a href="/" target="_top">Playground home</a></div>
      <p><a href="/runtime/control?run=$(run.id)">Worker control and diagnostics</a></p>
      <p class="lc-runtime-hint"><code>$(run.id)</code></p>
    </main><script type="module" src="/runtime/assets/run.js"></script></body></html>"""
    return gateway_response(stream, status, body; content_type="text/html; charset=utf-8")
end

"""
    control_surface(stream, supervisor, principal)

Render the protected worker control page using the same assets as the Bonito
components. An optional `run` query selects an authorized run; its role controls
come from the registered application definition. Rendering allocates no resources.
"""
function control_surface(stream, supervisor::UIHostSupervisor, principal::Principal)
    query = URIs.queryparampairs(URIs.URI(stream.message.target))
    length(query) <= 1 && all(p -> first(p) == "run", query) ||
        throw(AccessDenied(400, "Expected at most one run selection"))
    run = isempty(query) ? nothing :
        get_run(supervisor.store, principal, requested_uuid(last(only(query))))
    definition = run === nothing ? nothing : get(supervisor.registry.definitions, run.application, nothing)
    roles = definition === nothing ? () : Tuple((role=r.role, profiles=r.profiles) for r in definition.requirements)
    configuration = JSON3.write((kind="panel", run_id=run === nothing ? nothing : string(run.id), roles))
    runs = list_runs(supervisor.store, principal)
    links = join(("<li><a href=\"/runtime/control?run=$(item.id)\">$(html_text(item.application)) · $(item.id)</a> · $(html_text(item.state))</li>"
        for item in runs), "\n")
    selected = run === nothing ? "Inventory only · choose an owned run for role assignments." :
        "$(run.application) · $(run.id) · $(run.state)"
    body = """<!doctype html>
    <html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Worker control · LineCableModels</title>
    $(read(RUNTIME_THEME_INIT, String))
    <link rel="stylesheet" href="/runtime/assets/brand.css">
    <link rel="stylesheet" href="/runtime/assets/control-contract.css">
    <link rel="stylesheet" href="/runtime/assets/published-text.css">
    <link rel="stylesheet" href="/runtime/assets/forms.css">
    <link rel="stylesheet" href="/runtime/assets/data-views.css">
    <link rel="stylesheet" href="/runtime/assets/workspace.css">
    <link rel="stylesheet" href="/runtime/assets/run.css">
    <link rel="stylesheet" href="/runtime/assets/runtime-controls.css"></head>
    <body class="lc-runtime-page"><main class="lc-runtime-surface">
      <p class="lc-runtime-eyebrow">RUNTIME CONTROL</p><h1>Workers and owned runs</h1>
      <p>$(html_text(selected))</p>
      <label class="lc-field">Theme<select class="lc-control-select lc-form-control" data-lcm-theme-selector>
        <option value="system">System</option><option value="dark">Dark</option><option value="light">Light</option>
      </select></label>
      <details><summary>Choose an owned run</summary><ul>$links</ul>
        <a href="/runtime/control">Inventory only</a>
        <p class="lc-runtime-hint">Reload this page to refresh the run list. Worker status updates live.</p></details>
      <div class="lc-runtime-controls" data-lcm-runtime-controls="$(html_text(configuration))">Loading worker controls…</div>
      <nav class="lc-runtime-actions"><a href="/">Playground home</a><a href="/dev/">Developer gallery</a></nav>
    </main><script defer src="/runtime/assets/runtime-client.js"></script>
      <script defer src="/runtime/assets/runtime-controls.js"></script>
      <script defer src="/runtime/assets/control.js"></script></body></html>"""
    return gateway_response(stream, 200, body; content_type="text/html; charset=utf-8")
end
