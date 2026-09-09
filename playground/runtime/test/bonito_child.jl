# Real Bonito boundary fixture. State belongs to this process, not the gateway.
ccall(:getppid, Cint, ()) == parse(Int, ENV["LCM_UI_PARENT_PID"]) || exit(70)
using LineCableModelsPlayground, LineCableModelsPlaygroundProtocol, Bonito, JSON3, UUIDs
const LCM = LineCableModelsPlayground
const count = Observable(0)
const run_id = ENV["LCM_RUN_ID"]
counter_app = App() do session
    content = DOM.div(
        DOM.button("Increment"; id="increment", onclick=js"event => $(count).notify($(count).value + 1)"),
        DOM.output(count; id="count"), DOM.span(run_id; id="run-id"))
    Bonito.onload(session, content, js"element => document.documentElement.dataset.fixtureReady = 'true'")
    LCM.widget_shell("RUNTIME CONTRACT", "Owned counter", content)
end
workbench_app = LCM.TemplateWorkbench.app(; xray=true)
const WB = LCM.WorkbenchUI
struct RuntimeControlWorkbench <: WB.AbstractWorkbench end
WB.initialize(::RuntimeControlWorkbench, session) = (active=Observable(:controls),)
WB.handle!(::RuntimeControlWorkbench, state, action) = nothing
function WB.compose(::RuntimeControlWorkbench, state)
    client = RuntimeClient(LCM.UUID(run_id))
    selector = WorkerSelector(client, :parameters; profiles=("line-parameters",))
    return WB.Workbench(namespace=:runtime_controls_fixture,
        identity=WB.Identity("Runtime controls", "Browser integration fixture · no scientific roles"),
        navigation=WB.Sidebar(WB.NavGroup("Fixture", WB.NavItem(:controls, "Runtime controls")); active=state.active),
        workspace=WB.ViewStack(WB.View(:controls, "Runtime controls", DOM.div(selector,
            PreparationStatus(client, :parameters), WorkerDiagnostics(client),
            JuliaTerminal(client, :terminal))); active=state.active))
end
runtime_workbench = WB.workbench_app(RuntimeControlWorkbench(); xray=LCM.ComponentXRay.XRayPolicy(true))
# Actual Bonito input/result binding fixture. Only the browser test replaces this
# control's HTTP transport; this child does not grant leases or execute science.
job_view = App() do session
    client = RuntimeClient(UUID(run_id))
    parameters = Observable{Any}(Dict("value"=>3))
    job = ScientificJob(client, :parameters, "system.echo"; parameters)
    fence = AssignmentFence(string(uuid4()), run_id, "fixture-owner", "parameters",
        "worker-a", string(uuid4()), string(uuid4()), "fixture", "1.0.0", repeat("a",64), 1)
    execution = PreparedExecution(string(uuid4()),1,repeat("b",64))
    samples = Dict{String,Any}()
    for value in (3,4,5)
        request = new_job_request(job.operation, Dict("value"=>value); session_id=run_id)
        result = JobResult("1.0", request.job_id, job.operation, "1.0", request.input_hash,
            "fixture", fence.fingerprint, fence.worker_id, "miss", utc_timestamp(), utc_timestamp(),
            Dict{String,Any}("value"=>value), nothing, nothing, String[])
        samples[string(value)] = JSON3.read(encode_message(AssignedResult("2.0",fence,result,execution)))
    end
    display = map(session, job.result) do result
        result === nothing ? "No value" : string(result.value["value"], result.current ? " · current" : " · outdated")
    end
    content = DOM.div(
        DOM.div(; id="job-fixture-data", var"data-job-fixture"=JSON3.write(samples)),
        DOM.button("Input 3"; id="input-three", onclick=js"event => $(parameters).notify({value:3})"),
        DOM.button("Input 4"; id="input-four", onclick=js"event => $(parameters).notify({value:4})"),
        DOM.button("Input 5"; id="input-five", onclick=js"event => $(parameters).notify({value:5})"),
        DOM.button("Invalid draft"; id="input-invalid", onclick=js"event => $(parameters).notify([])"),
        DOM.div(DOM.input(;type="number",id="draft-field",value="3",
            oninput=js"event => $(parameters).notify({value:event.currentTarget.valueAsNumber})"),
            DOM.input(;type="checkbox",id="draft-checkbox",
                onchange=js"event => $(parameters).notify({value:event.currentTarget.checked ? 4 : 3})"),
            DOM.select(DOM.option("Three";value="3"),DOM.option("Five";value="5");id="draft-choice",
                onchange=js"event => $(parameters).notify({value:Number(event.currentTarget.value)})");
            var"data-runtime-input-fields"=""),
        DOM.output(display; id="scientific-display"), job; var"data-runtime-input-scope"="")
    LCM.widget_shell("BINDING CONTRACT", "Scientific display", content)
end
LCM.serve_owned_ui(("/counter"=>counter_app, "/workbench"=>workbench_app,
    "/runtime-controls"=>LCM.runtime_controls_widget(), "/runtime-workbench"=>runtime_workbench,
    "/job-view"=>job_view, "/terminal-view"=>LCM.runtime_terminal_widget(RuntimeClient(UUID(run_id)))))
