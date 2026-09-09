const RUNTIME_CONTROLS_PATH = joinpath(PLAYGROUND_ROOT, "assets", "runtime-controls.css")
const RUNTIME_CLIENT_SCRIPT_PATH = joinpath(PLAYGROUND_ROOT, "assets", "runtime-client.js")
const RUNTIME_CONTROLS_SCRIPT_PATH = joinpath(PLAYGROUND_ROOT, "assets", "runtime-controls.js")
include_dependency(RUNTIME_CONTROLS_PATH)
include_dependency(RUNTIME_CLIENT_SCRIPT_PATH)
include_dependency(RUNTIME_CONTROLS_SCRIPT_PATH)
const RUNTIME_CONTROLS_CSS = read(RUNTIME_CONTROLS_PATH, String)
const RUNTIME_CLIENT_SCRIPT = Bonito.Asset(RUNTIME_CLIENT_SCRIPT_PATH)
const RUNTIME_CONTROLS_SCRIPT = Bonito.Asset(RUNTIME_CONTROLS_SCRIPT_PATH)

"""
    RuntimeClient(run_id=nothing)

Identify an owned application run for same-origin browser controls. Construction
does not connect to a broker, poll HTTP, allocate resources, or retain credentials.
`nothing` permits inventory inspection but disables assignment. The gateway
authorizes every request; this value does not confer access to a run.
"""
struct RuntimeClient
    "Owned application-run identity, or inventory-only context."
    run_id::Union{Nothing,UUID}
end
RuntimeClient() = RuntimeClient(nothing)

abstract type AbstractRuntimeControl end

"""
    WorkerSelector(client, role; profiles)

Select an approved profile and automatic, pinned, or dedicated placement for one
declared application role. Assignment and release require explicit actions.
`profiles` contains stable profile IDs; it does not install or prepare them.
The coordinator independently validates the run's registered requirements.
"""
struct WorkerSelector{P<:Tuple} <: AbstractRuntimeControl
    "Same-origin owned-run context."
    client::RuntimeClient
    "Registered application role."
    role::String
    "Profile IDs offered for this role."
    profiles::P
    function WorkerSelector(client::RuntimeClient, role; profiles)
        name = string(role)
        ids = Tuple(string.(collect(profiles)))
        runtime_control_token(name)
        1 <= length(ids) <= 64 && allunique(ids) || throw(ArgumentError("expected 1–64 unique profile IDs"))
        foreach(runtime_control_token, ids)
        return new{typeof(ids)}(client, name, ids)
    end
end

"""
    PreparationStatus(client, role; parameters=Dict{String,Any}())

Display preparation evidence separately from worker connection and assignment.
Offer explicit preparation and cancellation using bounded passive `parameters`.
Without fresh executor evidence the component displays unknown, never ready.
Input values are sent only on preparation, not during status polling or to X-ray.
"""
struct PreparationStatus <: AbstractRuntimeControl
    "Same-origin owned-run context."
    client::RuntimeClient
    "Registered application role."
    role::String
    "Passive preparation inputs; excluded from diagnostic metadata."
    parameters::Dict{String,Any}
    function PreparationStatus(client::RuntimeClient, role; parameters::AbstractDict=Dict{String,Any}())
        name = string(role)
        runtime_control_token(name)
        inputs = normalize_wire(parameters)
        ncodeunits(JSON3.write(inputs)) <= 65536 || throw(ArgumentError("preparation inputs exceed 64 KiB"))
        return new(client, name, inputs)
    end
end

"""
    ScientificResult

Hold one successful display value with its existing assigned-result provenance.
`current` means the owned browser view still matches its input revision and
prepared target. This view projection never authorizes a job or a download.
"""
struct ScientificResult
    "Validated wire provenance of the displayed result."
    provenance::AssignedResult
    "Complete passive scientific value, excluded from X-ray metadata."
    value::Dict{String,Any}
    "Whether the display still matches current inputs and executor evidence."
    current::Bool
end
Base.show(io::IO, result::ScientificResult) = print(io, "ScientificResult(", result.provenance.result.job_id, ", <private view>)")

"""
    ScientificJob(client, role, operation; parameters=Dict{String,Any}())

Compose explicit Run/Cancel controls and last-successful result provenance using
the shared, same-origin runtime client. `parameters` is a passive input object or
an Observable containing one. Changing it marks the retained result outdated;
construction, rendering and input changes never submit or prepare work.

`result` is an Observable containing `nothing` or a `ScientificResult`. Compose
plots and tables from this binding without replacing their canvas or surrounding
DOM. Late results from superseded inputs/assignments are not projected. Create
the component and its input Observable in the owning Bonito session. No scientific
packages, broker credentials or user-code callbacks are installed by this control.

For browser-edited fields, put this control and its fields inside one
`data-runtime-input-scope` container and mark the fields' container with
`data-runtime-input-fields`. Run then waits for the latest field round trip;
late echoes cannot clear a newer edit. Display-only controls stay outside the
marked fields container. No scientific validation is duplicated in JavaScript.
"""
struct ScientificJob{I<:Observable} <: AbstractRuntimeControl
    "Same-origin owned-run context."
    client::RuntimeClient
    "Registered application role."
    role::String
    "Registered scientific operation, not a Julia expression."
    operation::String
    "Session-owned passive input object or projection from ordinary field values."
    parameters::I
    "Last successful display value and provenance; never an execution authority."
    result::Observable{Union{Nothing,ScientificResult}}
end

function scientific_job_parameters(value)
    inputs = normalize_wire(value)
    inputs isa Dict{String,Any} || throw(ArgumentError("scientific parameters require an object"))
    ncodeunits(JSON3.write(inputs)) <= 65536 || throw(ArgumentError("scientific inputs exceed 64 KiB"))
    return inputs
end

function ScientificJob(client::RuntimeClient, role, operation; parameters=Dict{String,Any}())
    name, op = string(role), string(operation)
    runtime_control_token(name); operation_subject_token(op)
    ncodeunits(op) <= 128 || throw(ArgumentError("operation exceeds 128 bytes"))
    inputs = parameters isa Observable ? parameters : Observable(scientific_job_parameters(parameters))
    scientific_job_parameters(inputs[])
    return ScientificJob(client, name, op, inputs, Observable{Union{Nothing,ScientificResult}}(nothing))
end
Base.show(io::IO, job::ScientificJob) = print(io, "ScientificJob(", job.role, ", ", job.operation, ", <private inputs>)")

function invalidate_job_projection!(job::ScientificJob)
    previous = job.result[]
    previous === nothing || !previous.current ||
        (job.result[] = ScientificResult(previous.provenance, previous.value, false))
    return nothing
end

# Browser projection is view data, not server authorization. Validate it again
# at the Bonito boundary and fence queued messages against the current draft.
function apply_job_projection!(job::ScientificJob, text::AbstractString, draft::String)
    next = try
        ncodeunits(text) <= 4 * 1024^2 + 262144 || return false
        packet = JSON3.read(text)
        Set(keys(packet)) == Set((:draft_id,:current,:provenance,:value)) || return false
        packet.draft_id == draft && packet.current isa Bool || return false
        envelope = decode_runtime_message(AssignedResult, JSON3.write(packet.provenance))
        envelope.fence.run_id == string(job.client.run_id) && envelope.fence.role == job.role &&
            envelope.result.operation == job.operation && envelope.result.failure === nothing || return false
        if !packet.current
            previous = job.result[]
            previous !== nothing && previous.provenance == envelope || return false
            invalidate_job_projection!(job)
            return true
        end
        envelope.result.input_hash == input_hash(job.operation, scientific_job_parameters(job.parameters[])) || return false
        value = normalize_wire(packet.value)
        value isa Dict{String,Any} || return false
        ScientificResult(envelope, value, true)
    catch
        return false
    end
    job.result[] = next
    return true
end

"""
    WorkerDiagnostics(client)

Display worker approval, presence, capacity, and bounded structured control events.
Events are authorized by the gateway and fetched independently of inventory;
arbitrary worker logs, terminal input, and credentials are not inspected.
"""
struct WorkerDiagnostics <: AbstractRuntimeControl
    "Same-origin owned-run or inventory-only context."
    client::RuntimeClient
end

"""
    WorkerControlPanel(client, selectors...)

Compose the shared selectors, preparation status, administrator registration
controls, and diagnostics. All selectors must refer to `client`'s run and have
distinct roles. An empty selector list provides inventory and administration only.
Registration controls are hidden for non-administrators and server-authorized.
"""
struct WorkerControlPanel{S<:Tuple} <: AbstractRuntimeControl
    "Same-origin owned-run or inventory-only context."
    client::RuntimeClient
    "Role selectors reused by the panel renderer."
    selectors::S
    function WorkerControlPanel(client::RuntimeClient, selectors::WorkerSelector...)
        length(selectors) <= 64 && allunique(s.role for s in selectors) ||
            throw(ArgumentError("expected at most 64 distinct worker roles"))
        all(s.client.run_id == client.run_id for s in selectors) ||
            throw(ArgumentError("all worker selectors must belong to the panel's run"))
        return new{typeof(selectors)}(client, selectors)
    end
end

function runtime_control_token(value::String)
    occursin(r"^[a-z0-9][a-z0-9_-]{0,63}$", value) || throw(ArgumentError("invalid runtime control identity"))
    return value
end

runtime_control_options(c::WorkerSelector) = (kind="selector", role=c.role, profiles=c.profiles)
runtime_control_options(c::PreparationStatus) = (kind="preparation", role=c.role, profiles=(), parameters=c.parameters)
runtime_control_options(c::ScientificJob) = (kind="execution", role=c.role, operation=c.operation,
    parameters=scientific_job_parameters(c.parameters[]))
runtime_control_options(::WorkerDiagnostics) = (kind="diagnostics",)
runtime_control_options(c::WorkerControlPanel) = (kind="panel",
    roles=Tuple((role=s.role, profiles=s.profiles) for s in c.selectors))

runtime_control_actions(::AbstractRuntimeControl) = (:refresh,)
runtime_control_actions(::WorkerSelector) = (:refresh, :assign, :release)
runtime_control_actions(::PreparationStatus) = (:refresh, :prepare, :cancelScientific)
runtime_control_actions(::ScientificJob) = (:submitJob, :jobResult, :jobArtifact, :cancelJob)
runtime_control_actions(::WorkerControlPanel) = (:refresh, :assign, :release, :prepare, :cancelScientific, :enroll, :registration)

runtime_inspection_options(c::AbstractRuntimeControl) = runtime_control_options(c)
runtime_inspection_options(c::ScientificJob) = (kind="execution", role=c.role, operation=c.operation)

function ComponentXRay.inspection(component::AbstractRuntimeControl)
    options = runtime_inspection_options(component)
    return ComponentXRay.ComponentInspection(component;
        name=string(nameof(typeof(component))),
        source=ComponentXRay.source_reference(@__MODULE__, @__FILE__, @__LINE__),
        parameters=[ComponentXRay.PropertyInspection(name, value) for (name, value) in pairs(options) if name != :parameters],
        actions=[ComponentXRay.ActionInspection(action, "Same-origin RuntimeClient.$action", nothing)
            for action in runtime_control_actions(component)],
        css_scopes=[".lc-runtime-controls", ".lc-runtime-section", ".lc-runtime-fields",
            ".lc-runtime-records", ".lc-runtime-events", ".lc-runtime-event-log", ".lc-runtime-note"],
        notes=["Shared browser renderer; gateway authorization is authoritative.",
            "Credentials, event payloads, and other owners' state are excluded from X-ray."])
end

runtime_control_binding(::Session, ::AbstractRuntimeControl, target) = js"(element, control) => {}"

function runtime_control_binding(session::Session, component::ScientificJob, target)
    draft = Ref(string(uuid4()))
    previous = Ref{Union{Nothing,String}}(nothing)
    function project_inputs(raw)
        inputs = try scientific_job_parameters(raw) catch; nothing end
        key = inputs === nothing ? "invalid" : input_hash(component.operation, inputs)
        if key != previous[]
            draft[] = string(uuid4())
            previous[] = key
            invalidate_job_projection!(component)
        end
        return (draft_id=draft[], parameters=inputs)
    end
    wire = Observable{@NamedTuple{draft_id::String, parameters::Union{Nothing,Dict{String,Any}}}}(
        project_inputs(component.parameters[]))
    map!(project_inputs, session, wire, component.parameters; update=false)
    input_epoch = Observable(0)
    input_ack = Observable((epoch=0, inputs=wire[]))
    on(session, input_epoch) do epoch
        epoch isa Integer && !(epoch isa Bool) && 0 < epoch <= 9_007_199_254_740_991 || return
        # The scoped bubbling event follows the native fields' input messages
        # on this session's ordered socket. Echo the resulting canonical draft.
        input_ack[] = (epoch=Int(epoch), inputs=wire[])
        return nothing
    end
    projected = Observable("")
    projection_status = Observable("idle")
    on(session, projected) do value
        projection_status[] = apply_job_projection!(component, value, draft[]) ? "accepted" : "rejected"
        return nothing
    end
    onjs(session, projection_status, js"value => { $(target).dataset.resultProjection = value; }")
    onjs(session, wire, js"value => $(target).__lcmJobInputs?.(value)")
    onjs(session, input_ack, js"value => $(target).__lcmJobInputAck?.(value)")
    return js"""
        (element, control) => {
            if (!control.job || element.__lcmJobInputs) return;
            element.dataset.resultProjection = $(projection_status).value;
            let draftId, previous, epoch = 0, confirmedEpoch = 0;
            element.__lcmJobInputs = value => {
                if (epoch !== confirmedEpoch) return;
                draftId = value.draft_id;
                try { control.job.setInputs(value.parameters); } catch {}
            };
            element.__lcmJobInputAck = value => {
                if (value.epoch !== epoch || value.epoch <= confirmedEpoch) return;
                confirmedEpoch = epoch;
                element.__lcmJobInputs(value.inputs);
            };
            element.__lcmJobInputs($(wire).value);
            const scope = element.closest('[data-runtime-input-scope]');
            if (scope) {
                const ownedInput = event => {
                    const input = event.target;
                    if (!(input instanceof Element) || !input.closest('[data-runtime-input-fields]') ||
                        input.closest('[data-runtime-input-scope]') !== scope) return false;
                    const onChange = input.tagName === 'SELECT' || ['checkbox', 'radio'].includes(input.type);
                    return event.type === (onChange ? 'change' : 'input');
                };
                const begin = event => {
                    if (!ownedInput(event) || control.job.closed) return;
                    epoch++;
                    control.job.markInputsPending();
                };
                const acknowledge = event => {
                    if (ownedInput(event) && !control.job.closed) $(input_epoch).notify(epoch);
                };
                for (const kind of ['input', 'change']) {
                    scope.addEventListener(kind, begin, true);
                    scope.addEventListener(kind, acknowledge);
                }
                const destroy = control.destroy.bind(control);
                control.destroy = () => {
                    for (const kind of ['input', 'change']) {
                        scope.removeEventListener(kind, begin, true);
                        scope.removeEventListener(kind, acknowledge);
                    }
                    destroy();
                };
            }
            control.job.subscribe(state => {
                const good = state.lastGood;
                if (!good) return;
                const key = good.receipt.id + ':' + state.current;
                if (key === previous) return;
                previous = key;
                $(projected).notify(JSON.stringify({draft_id:draftId, current:state.current,
                    provenance:good.provenance, value:good.value}));
            });
        }
    """
end

function Bonito.jsrender(session::Session, component::AbstractRuntimeControl)
    client_script = DOM.script(src=RUNTIME_CLIENT_SCRIPT)
    controls_script = DOM.script(src=RUNTIME_CONTROLS_SCRIPT)
    configuration = JSON3.write(merge(runtime_control_options(component),
        (run_id=component.client.run_id === nothing ? nothing : string(component.client.run_id),)))
    target = DOM.div("Loading runtime controls…";
        class="lc-runtime-controls", var"data-lcm-runtime-controls"=configuration)
    binding = runtime_control_binding(session, component, target)
    Bonito.onload(session, target, js"""
        (element) => {
            const clientScript = $(client_script), controlsScript = $(controls_script);
            const start = () => {
                if (!element.isConnected || !globalThis.LineCableModelsRuntimeClient ||
                    !globalThis.LineCableModelsRuntimeControls) return;
                try {
                    const control = globalThis.LineCableModelsRuntimeControls.mount(element);
                    ($(binding))(element, control);
                }
                catch { element.textContent = "Runtime controls could not be initialized."; }
            };
            clientScript.addEventListener("load", start, {once:true});
            controlsScript.addEventListener("load", start, {once:true});
            const failed = () => { element.textContent = "Runtime control library is unavailable."; };
            clientScript.addEventListener("error", failed, {once:true});
            controlsScript.addEventListener("error", failed, {once:true});
            start();
        }
    """)
    return Bonito.jsrender(session, DOM.div(
        DOM.style(RUNTIME_CONTROLS_CSS; var"data-lcm-css-source"="assets/runtime-controls.css"),
        client_script, controls_script, ComponentXRay.instrument(session, target, component);
        style="display: contents;"))
end
