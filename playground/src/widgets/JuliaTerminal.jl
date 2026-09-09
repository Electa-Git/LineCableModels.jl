const TERMINAL_ASSET_NAMES = ("runtime-client.js", "vendor/runtime-terminal.bundle.js",
    "runtime-terminal-client.js", "runtime-terminal.js")
const TERMINAL_SCRIPT_ASSETS = map(TERMINAL_ASSET_NAMES) do name
    path = joinpath(PLAYGROUND_ROOT, "assets", name)
    include_dependency(path)
    Bonito.Asset(path)
end
const TERMINAL_STYLE_NAMES = ("runtime-controls.css", "vendor/runtime-terminal.bundle.css", "runtime-terminal.css")
const TERMINAL_STYLES = map(TERMINAL_STYLE_NAMES) do name
    path = joinpath(PLAYGROUND_ROOT, "assets", name)
    include_dependency(path)
    read(path, String)
end

"""
    JuliaTerminal(client, role; title="Julia REPL", rows=18)

Render a private, disposable Julia terminal for one registered application role.
Compose [`WorkerSelector`](@ref) separately to select its container-backed worker.
Construction and rendering allocate no terminal; connection, interruption, stop
and restart are explicit browser actions. The gateway authorizes each action.

`title` is a display label of at most 120 characters. `rows` sets the initial
viewport height in terminal rows, between 6 and 60. Width follows its parent.
The shared theme controls the terminal palette in all presentation surfaces.

Input and output travel over the private same-origin socket, never through
Bonito Observables, X-ray metadata or scientific jobs. A lost acknowledgement
stops input without replay. Reconnect permits output inspection; explicitly
resume input or restart after reviewing an unconfirmed action. Reconnection may retain memory only within the
worker's grace period; reload, stop and restart do not restore Julia variables.
"""
struct JuliaTerminal <: AbstractRuntimeControl
    "Same-origin owned application-run context."
    client::RuntimeClient
    "Registered terminal role, not an executable command."
    role::String
    "Public display label."
    title::String
    "Initial viewport height in terminal rows."
    rows::Int
    function JuliaTerminal(client::RuntimeClient, role; title::AbstractString="Julia REPL", rows::Integer=18)
        name = runtime_control_token(string(role))
        1 <= length(title) <= 120 && !any(iscntrl, title) || throw(ArgumentError("expected a short terminal title"))
        !(rows isa Bool) && 6 <= rows <= 60 || throw(ArgumentError("terminal rows must be between 6 and 60"))
        new(client, name, String(title), Int(rows))
    end
end

runtime_control_options(c::JuliaTerminal) = (kind="terminal", role=c.role, title=c.title, rows=c.rows)
runtime_control_actions(::JuliaTerminal) = (:connect, :disconnect, :interrupt, :stop, :restart, :clear, :resumeInput)

function ComponentXRay.inspection(component::JuliaTerminal)
    ComponentXRay.ComponentInspection(component;
        name="JuliaTerminal", source=ComponentXRay.source_reference(@__MODULE__, @__FILE__, @__LINE__),
        parameters=[ComponentXRay.PropertyInspection(name, value) for (name, value) in pairs(runtime_control_options(component))],
        actions=[ComponentXRay.ActionInspection(action, action == :clear ? "xterm.Terminal.clear (local display)" :
            action in (:stop, :restart) ? "RuntimeTerminal.control(:$action)" : "RuntimeTerminal.$action", nothing)
            for action in runtime_control_actions(component)],
        css_scopes=[".lc-runtime-controls", ".lc-runtime-terminal", ".lc-terminal-heading", ".lc-terminal-phase",
            ".lc-terminal-status", ".lc-terminal-viewport", ".lc-terminal-screen", ".lc-terminal-confirm"],
        notes=["No terminal input, output, writer identity, credentials or history is exposed to X-ray.",
            "Gateway authorization and worker resource limits remain authoritative."])
end

function Bonito.jsrender(session::Session, component::JuliaTerminal)
    scripts = map(asset -> DOM.script(src=asset), TERMINAL_SCRIPT_ASSETS)
    configuration = JSON3.write(merge(runtime_control_options(component),
        (run_id=component.client.run_id === nothing ? nothing : string(component.client.run_id),)))
    target = DOM.div("Loading private terminal controls…"; class="lc-runtime-controls lc-runtime-terminal",
        var"data-lcm-runtime-terminal"=configuration)
    Bonito.onload(session, target, js"""
        element => {
            const scripts = $(scripts);
            const start = () => {
                if (!element.isConnected || !globalThis.LineCableModelsRuntimeClient ||
                    !globalThis.LineCableModelsTerminalVendor || !globalThis.LineCableModelsTerminalClient ||
                    !globalThis.LineCableModelsTerminal) return;
                try { globalThis.LineCableModelsTerminal.mount(element); }
                catch { element.textContent = "Private terminal controls could not be initialized."; }
            };
            for (const script of scripts) {
                script.addEventListener("load", start, {once:true});
                script.addEventListener("error", () => { element.textContent = "Private terminal library is unavailable."; }, {once:true});
            }
            start();
        }
    """)
    styles = map(zip(TERMINAL_STYLE_NAMES, TERMINAL_STYLES)) do (name, content)
        DOM.style(content; var"data-lcm-css-source"="assets/" * name)
    end
    Bonito.jsrender(session, DOM.div(styles..., scripts..., ComponentXRay.instrument(session, target, component);
        style="display: contents;"))
end
