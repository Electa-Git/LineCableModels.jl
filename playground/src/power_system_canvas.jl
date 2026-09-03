const POWER_SYSTEM_CANVAS_STYLES_PATH = normpath(
    joinpath(@__DIR__, "..", "assets", "power-system-canvas.css")
)
const POWER_SYSTEM_CANVAS_LIBRARY_STYLES_PATH = normpath(
    joinpath(@__DIR__, "..", "assets", "vendor", "power-system-canvas.bundle.css")
)
const POWER_SYSTEM_CANVAS_SCRIPT_PATH = normpath(
    joinpath(@__DIR__, "..", "assets", "vendor", "power-system-canvas.bundle.js")
)
include_dependency(POWER_SYSTEM_CANVAS_STYLES_PATH)
include_dependency(POWER_SYSTEM_CANVAS_LIBRARY_STYLES_PATH)
const POWER_SYSTEM_CANVAS_STYLES = read(POWER_SYSTEM_CANVAS_STYLES_PATH, String)
const POWER_SYSTEM_CANVAS_LIBRARY_STYLES = read(
    POWER_SYSTEM_CANVAS_LIBRARY_STYLES_PATH,
    String
)
const POWER_SYSTEM_CANVAS_SCRIPT = Bonito.Asset(POWER_SYSTEM_CANVAS_SCRIPT_PATH)

const POWER_SYSTEM_DIAGRAM_MAX_BYTES = 2 * 1024 * 1024
const POWER_SYSTEM_TOPOLOGY_MAX_ITEMS = 10_000

"An editable, browser-local power-system single-line diagram surface."
struct PowerSystemCanvas
    diagram::String
    read_only::Bool
    topology::Observable{String}
    bridge_payload::Observable{String}
end

function require_json_array(document, key::Symbol, context::AbstractString)
    hasproperty(document, key) || throw(ArgumentError("$context is missing `$key`"))
    value = getproperty(document, key)
    value isa JSON3.Array || throw(ArgumentError("$context `$key` must be an array"))
    length(value) <= POWER_SYSTEM_TOPOLOGY_MAX_ITEMS || throw(ArgumentError(
        "$context `$key` exceeds the item limit"
    ))
    return value
end

function validate_power_system_diagram(payload::AbstractString)
    ncodeunits(payload) <= POWER_SYSTEM_DIAGRAM_MAX_BYTES || throw(ArgumentError(
        "power-system diagram exceeds the $(POWER_SYSTEM_DIAGRAM_MAX_BYTES)-byte limit"
    ))
    document = JSON3.read(payload)
    hasproperty(document, :version) || throw(ArgumentError(
        "power-system diagram is missing `version`"
    ))
    string(document.version) == "1" || throw(ArgumentError(
        "power-system diagram version must be `1`"
    ))
    require_json_array(document, :elements, "power-system diagram")
    for key in (:buses, :junctions, :wires, :annotations, :customKinds)
        hasproperty(document, key) || continue
        require_json_array(document, key, "power-system diagram")
    end
    return String(payload)
end

function validate_power_system_topology(payload::AbstractString)
    ncodeunits(payload) <= POWER_SYSTEM_DIAGRAM_MAX_BYTES || throw(ArgumentError(
        "compiled topology exceeds the $(POWER_SYSTEM_DIAGRAM_MAX_BYTES)-byte limit"
    ))
    document = JSON3.read(payload)
    hasproperty(document, :schema) || throw(ArgumentError(
        "compiled topology is missing `schema`"
    ))
    string(document.schema) == "lcm.power-system-topology/1" || throw(ArgumentError(
        "unsupported compiled-topology schema"
    ))
    hasproperty(document, :diagram) || throw(ArgumentError(
        "compiled topology is missing `diagram`"
    ))
    validate_power_system_diagram(JSON3.write(document.diagram))
    hasproperty(document, :topology) || throw(ArgumentError(
        "compiled topology is missing `topology`"
    ))
    topology = document.topology
    for key in (:nodes, :terminal_to_node, :diagnostics)
        require_json_array(topology, key, "compiled topology")
    end
    return String(payload)
end

function PowerSystemCanvas(
        diagram;
        read_only::Bool=false
    )
    encoded = diagram isa AbstractString ? String(diagram) : JSON3.write(diagram)
    validate_power_system_diagram(encoded)
    topology = Observable("")
    bridge_payload = Observable("")
    on(bridge_payload) do payload
        isempty(payload) && return nothing
        accepted = try
            validate_power_system_topology(payload)
        catch error
            @warn "Rejected invalid browser power-system topology" exception=(
                error,
                catch_backtrace()
            )
            return nothing
        end
        topology[] = accepted
        return nothing
    end
    return PowerSystemCanvas(encoded, read_only, topology, bridge_payload)
end

"Return the latest validated browser topology, or `nothing` before mounting."
function snapshot(component::PowerSystemCanvas)
    isempty(component.topology[]) && return nothing
    return JSON3.read(component.topology[], Dict{String,Any})
end

function Bonito.jsrender(session::Session, component::PowerSystemCanvas)
    onboarding = DOM.script("""
        try { window.localStorage.setItem('ole-onboarding-dismissed', '1'); }
        catch (_) {}
    """)
    library = DOM.script(src=POWER_SYSTEM_CANVAS_SCRIPT)
    reset = DOM.button(
        "Reset specimen";
        type="button",
        var"data-power-role"="reset",
        class="lc-power-system-reset"
    )
    summary = DOM.output(
        "Compiling topology…";
        var"data-power-role"="summary",
        class="lc-power-system-summary"
    )
    editor_mount = DOM.div(
        var"data-power-role"="editor",
        class="lc-power-system-editor-mount"
    )
    root = DOM.section(
        DOM.style(POWER_SYSTEM_CANVAS_LIBRARY_STYLES),
        DOM.style(POWER_SYSTEM_CANVAS_STYLES),
        onboarding,
        library,
        DOM.header(
            DOM.div(
                DOM.span("IEC single-line workspace"; class="lc-power-system-label"),
                summary;
                class="lc-power-system-copy"
            ),
            reset;
            class="lc-power-system-command-bar"
        ),
        editor_mount,
        DOM.footer(
            "Diagram edits remain browser-local. Julia receives a validated topology snapshot; no calculation is dispatched.";
            class="lc-power-system-boundary"
        );
        class="lc-power-system-canvas"
    )
    configuration = JSON3.write((
        diagram=JSON3.read(component.diagram),
        readOnly=component.read_only,
    ))
    configuration_hex = bytes2hex(codeunits(configuration))

    Bonito.onload(session, root, js"""
        (element) => {
            const library = $(library);
            const editor = $(editor_mount);
            const reset = $(reset);
            const summary = $(summary);
            const topology = $(component.bridge_payload);
            const configurationBytes = Uint8Array.from(
                $configuration_hex.match(/../g),
                byte => Number.parseInt(byte, 16),
            );
            const configuration = JSON.parse(
                new TextDecoder().decode(configurationBytes),
            );
            let mounted = false;

            const start = () => {
                if (mounted || !globalThis.LCMPowerSystemCanvas) return;
                mounted = true;
                try {
                    element._lcmPowerSystemCanvas = globalThis.LCMPowerSystemCanvas.mount(
                        editor,
                        configuration,
                        topology,
                        summary,
                    );
                    reset.addEventListener("click", () => {
                        element._lcmPowerSystemCanvas?.reset();
                    });
                } catch (error) {
                    summary.dataset.kind = "error";
                    summary.textContent = error instanceof Error
                        ? error.message
                        : "Unable to initialize the single-line editor";
                }
            };

            if (globalThis.LCMPowerSystemCanvas) start();
            else {
                library.addEventListener("load", start, {once: true});
                library.addEventListener("error", () => {
                    summary.dataset.kind = "error";
                    summary.textContent = "Unable to load the single-line editor";
                }, {once: true});
            }

            window.addEventListener("pagehide", () => {
                element._lcmPowerSystemCanvas?.destroy();
            }, {once: true});
        }
    """)
    return Bonito.jsrender(session, root)
end
