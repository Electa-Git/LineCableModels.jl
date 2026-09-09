import assert from "node:assert/strict";
import vm from "node:vm";
import {readFile} from "node:fs/promises";
import {setTimeout as delay} from "node:timers/promises";
import {webcrypto} from "node:crypto";

const script = await readFile(new URL("../../assets/runtime-client.js", import.meta.url), "utf8");
const run = webcrypto.randomUUID(), lease = webcrypto.randomUUID();
const initialControl = () => ({schema_version: 1, enabled: true, broker: "online",
  administrator: true, workers: [], profiles: [], provisioned: []});
const json = (body, status = 200) => new Response(JSON.stringify(body), {
  status, headers: {"Content-Type": "application/json"}});
const assigned = () => ({id: lease, run_id: run, role: "main", profile: "line-parameters",
  state: "reserving", usable: false});
async function until(check, ms = 2000) {
  const deadline = Date.now() + ms;
  while (!check() && Date.now() < deadline) await delay(5);
  assert.ok(check(), "condition did not complete within its bound");
}
function environment(fetcher) {
  const document = new EventTarget();
  document.visibilityState = "visible";
  const context = vm.createContext({fetch: fetcher, document, crypto: webcrypto,
    AbortController, Response, TextDecoder, TextEncoder, Uint8Array, setTimeout, clearTimeout});
  vm.runInContext(script, context);
  return {api: context.LineCableModelsRuntimeClient, document, context};
}

// An owned UI run can exist on a publisher whose worker control is disabled.
{
  const fixture = environment(async url => url.endsWith("/control") ?
    json({...initialControl(), enabled:false, broker:"disabled"}) : json({error:"Worker control is not configured"},503));
  const client = fixture.api.acquire(run);
  await client.refresh();
  assert.equal(client.state.stale,false);
  assert.equal(client.state.control.enabled,false);
  assert.equal(client.state.assignments.length,0);
  assert.throws(() => client.assign("main","line-parameters",{mode:"automatic"}));
  client.close();
}

// No network on load/construction, one client per run, one coalesced poll.
{
  let calls = 0, down = false;
  const fixture = environment(async (url) => {
    calls++;
    if (down) throw new Error("private broker password must never be displayed");
    await delay(5);
    return json(url.endsWith("/control") ? initialControl() : []);
  });
  assert.equal(calls, 0);
  const client = fixture.api.acquire(run);
  assert.equal(fixture.api.acquire(run), client);
  assert.equal(calls, 0);
  fixture.document.dispatchEvent(new Event("visibilitychange"));
  await delay(10);
  assert.equal(calls, 0);
  await Promise.all([client.refresh(), client.refresh(), client.refresh()]);
  assert.equal(calls, 2);
  assert.equal(client.state.stale, false);
  assert.ok(Object.isFrozen(client.state));
  const good = client.state.control;
  down = true;
  await client.refresh();
  assert.equal(client.state.stale, true);
  assert.equal(client.state.control, good);
  assert.ok(!client.state.error.includes("password"));
  await assert.rejects(() => client.assign("main", "line-parameters", {mode: "automatic"}));
  client.close();
  assert.notEqual(fixture.api.acquire(run), client);
  fixture.api.acquire(run).close();
  assert.throws(() => new fixture.api.RuntimeClient("../foreign"));
  assert.throws(() => new fixture.api.RuntimeClient(null, {timeoutMs: 0}));
}

// Log transport is independent of inventory; cursor, gap and bounded history
// survive a log failure without creating a second discovery implementation.
{
  let controls = 0, stalled = true;
  let epoch = webcrypto.randomUUID(), sequence = 1;
  const queries = [];
  const fixture = environment(async (url, options) => {
    if (url.endsWith("/control")) { controls++; return json(initialControl()); }
    assert.ok(url.includes("/control/events"));
    queries.push(url);
    if (stalled) return new Promise((resolve, reject) => {
      options.signal.addEventListener("abort", () => reject(new Error("cancelled fixture")), {once: true});
    });
    const query = new URL(url, "http://fixture.test").searchParams;
    const changed = query.has("epoch") && query.get("epoch") !== epoch;
    const after = changed ? 0 : Number(query.get("after") || 0);
    return json({epoch, cursor: sequence, gap: changed, evicted: 0,
      records: Array.from({length: sequence - after}, (_, i) => ({sequence: after + i + 1, code: "control_connected"}))});
  });
  const client = new fixture.api.RuntimeClient(null, {pollMs: 50, timeoutMs: 120});
  let delivered = 0;
  const listener = () => { delivered++; };
  const remove = client.subscribe(listener, {withEvents: true});
  assert.throws(() => client.subscribe(listener));
  await until(() => client.state.eventsStale);
  assert.ok(controls >= 3, "log timeout must not delay inventory polls");
  assert.equal(client.state.stale, false);
  stalled = false;
  await until(() => client.state.events?.records.length === 1 && !client.state.eventsStale);
  await until(() => queries.some(q => q.includes("after=1&epoch=" + epoch)));
  sequence = 600;
  await until(() => client.state.events?.cursor === 600);
  assert.equal(client.state.events.records.length, 512);
  assert.equal(client.state.events.localDropped, 88);
  epoch = webcrypto.randomUUID();
  sequence = 601;
  await until(() => client.state.events?.epoch === epoch);
  assert.equal(client.state.events.gap, true);
  assert.equal(client.state.events.records.length, 512);
  const beforeHidden = controls;
  fixture.document.visibilityState = "hidden";
  await delay(160);
  assert.ok(controls <= beforeHidden + 1);
  fixture.document.visibilityState = "visible";
  fixture.document.dispatchEvent(new Event("visibilitychange"));
  await until(() => controls > beforeHidden + 1);
  const beforeClose = delivered;
  remove(); remove();
  await delay(150);
  assert.equal(delivered, beforeClose);
  assert.equal(client.closed, true);
  assert.equal(client.controllers.size, 0);
}

// A lost mutation reply is never replayed. Its explicit retry retains the
// caller's exact idempotency identity and discovers the already-owned result.
{
  let posts = 0, effects = 0, owned = [], loseReply = true;
  const seen = new Set(), requests = [];
  const fixture = environment(async (url, options) => {
    requests.push({url, ...options});
    if (options.method === "GET") return json(url.endsWith("/control") ? initialControl() : owned);
    posts++;
    const input = JSON.parse(options.body);
    if (!seen.has(input.request_id)) { seen.add(input.request_id); effects++; owned = [assigned()]; }
    if (loseReply) { loseReply = false; throw new Error("reply was lost after commit"); }
    return json(owned[0], 202);
  });
  const client = new fixture.api.RuntimeClient(run);
  await client.refresh();
  const requestId = webcrypto.randomUUID();
  await assert.rejects(() => client.assign("main", "line-parameters", {mode: "pinned", worker_id: "worker-a"}, {requestId}),
    error => error.uncertain && error.requestId === requestId);
  assert.equal(posts, 1);
  assert.equal(client.state.assignments[0].id, lease);
  assert.equal(client.state.pending, false);
  await client.assign("main", "line-parameters", {mode: "pinned", worker_id: "worker-a"}, {requestId});
  assert.equal(posts, 2);
  assert.equal(effects, 1);
  for (const request of requests) {
    assert.ok(request.url.startsWith("/runtime/api/"));
    assert.equal(request.credentials, "same-origin");
    assert.equal(request.redirect, "error");
    if (request.method !== "GET") assert.equal(request.headers["X-LCM-Request"], "1");
  }
  assert.throws(() => client.assign("main", "line-parameters", {mode: "pinned", worker_id: "../bad"}));
  client.close();
}

// Reconciliation cannot use a pre-action in-flight snapshot.
{
  let owned = [], block = false, resume = null;
  const fixture = environment(async (url, options) => {
    if (url.endsWith("/control")) return json(initialControl());
    if (options.method !== "GET") { owned = [assigned()]; return json(owned[0], 202); }
    if (block) {
      block = false;
      const previous = structuredClone(owned);
      return new Promise(resolve => { resume = () => resolve(json(previous)); });
    }
    return json(owned);
  });
  const client = new fixture.api.RuntimeClient(run);
  await client.refresh();
  block = true;
  const old = client.refresh();
  await until(() => resume !== null);
  const mutation = client.assign("main", "line-parameters", {mode: "automatic"});
  await until(() => owned.length === 1);
  assert.equal(client.state.pending, true);
  resume();
  await old;
  await mutation;
  assert.equal(client.state.assignments[0].id, lease);
  assert.equal(client.state.pending, false);
  client.close();
}

// Bound responses before parsing, and retain no active request after teardown.
{
  const fixture = environment(async () => new Response("{}", {
    headers: {"Content-Length": String(4 * 1024 * 1024 + 1)}}));
  const client = new fixture.api.RuntimeClient();
  await client.refresh();
  assert.equal(client.state.stale, true);
  assert.equal(client.state.control, null);
  assert.equal(client.controllers.size, 0);
  client.close();
}
// Preparation polling is independent; finite evidence cannot survive latency,
// a mutation, a hidden page, a lost lease, or a stale inventory response.
{
  let controls=0, reads=0, fail=false, usable=true, block=false, resume=null;
  let ttl=180, slow=0;
  const mutations=[];
  const report=() => ({channel:"online",phase:"idle",preparation:"ready",valid_for_ms:ttl,
    pending:false,accepted:true,progress:1,elapsed_seconds:1,output_lines:0,
    executor_id:lease,executor_generation:1,current_request_id:null,preparation_key:"b".repeat(64)});
  const fixture=environment(async (url, options) => {
    if (url.endsWith("/control")) {controls++; return json({...initialControl(),preparation_control:true});}
    if (url.endsWith("/assignments")) return json([{...assigned(),state:"active",usable}]);
    assert.ok(url.endsWith("/assignments/"+lease+"/science"));
    if (options.method!=="GET") {mutations.push(JSON.parse(options.body));return json({accepted:true},202);}
    reads++;
    if (fail) throw new Error("private child output");
    if (block) {block=false;await new Promise(resolve=>{resume=resolve;});}
    if (slow) await delay(slow);
    return json(report());
  });
  const client=new fixture.api.RuntimeClient(run,{pollMs:1000,timeoutMs:1000});
  const remove=client.subscribe(()=>{}, {withScience:true});
  await until(()=>client.state.science[lease]?.preparation==="ready");
  await until(()=>client.state.science[lease]?.preparation==="unknown");
  assert.equal(client.state.science[lease].preparation_key,null,"finite evidence expires without a poll");
  assert.equal(reads,1);
  ttl=30;slow=50;await client.refreshScience();slow=0;
  assert.equal(client.state.science[lease].preparation,"unknown","HTTP latency cannot extend ready time");
  ttl=500;fail=true;await client.refreshScience();fail=false;
  assert.equal(client.state.science[lease].channel,"offline");
  assert.equal(client.state.stale,false,"science outage is not an inventory outage");
  block=true;const previous=client.refreshScience();await until(()=>resume!==null);
  const before=controls;await client.refresh();assert.ok(controls>before,"stalled science must not stall inventory");
  const requestId=webcrypto.randomUUID();
  await client.prepare(lease,{length:2},{requestId});
  assert.equal(mutations.length,1);assert.equal(mutations[0].request_id,requestId);
  assert.equal(mutations[0].action,"prepare");
  resume();await previous;
  assert.equal(client.state.science[lease],undefined,"pre-action report must not restore readiness");
  await client.refreshScience();assert.equal(client.state.science[lease].preparation,"ready");
  fixture.document.visibilityState="hidden";fixture.document.dispatchEvent(new Event("visibilitychange"));
  assert.equal(Object.keys(client.state.science).length,0);
  fixture.document.visibilityState="visible";await client.refreshScience();
  assert.equal(client.state.science[lease].preparation,"ready");
  usable=false;await client.refresh();assert.equal(Object.keys(client.state.science).length,0);
  await client.cancelScientific(lease,requestId);
  assert.equal(mutations.at(-1).target_id,requestId);
  assert.equal(mutations.at(-1).action,"cancel");
  assert.throws(()=>client.prepare(lease,[]));
  assert.throws(()=>client.prepare(lease,{large:"x".repeat(65536)}));
  assert.throws(()=>client.cancelScientific(lease,"foreign"));
  remove();assert.equal(client.closed,true);
  assert.equal(client.controllers.size,0);
}
// Job requests use the same protected client, immutable retry identity and exact
// result provenance. Failed/revoked jobs remain readable but are never replayed.
{
  const jobId=webcrypto.randomUUID(), boot=webcrypto.randomUUID(), executor=webcrypto.randomUUID();
  const receipt={id:jobId,request_id:webcrypto.randomUUID(),run_id:run,lease_id:lease,role:"main",
    operation:"fixture.echo",input_hash:"a".repeat(64),worker_id:"worker-a",worker_boot:boot,generation:1,
    executor_id:executor,executor_generation:1,preparation_key:"b".repeat(64),
    submitted_at:"2026-09-08T00:00:00.000Z",deadline:"2026-09-08T00:01:00.000Z",state:"queued",
    current_assignment:true,channel:"online",cancel_requested:false,cancel_acknowledged:false};
  const outcome={protocol_version:"2.0",fence:{lease_id:lease,run_id:run,role:"main",worker_id:"worker-a",
    worker_boot:boot,generation:1,fingerprint:"c".repeat(64)},
    execution:{executor_id:executor,executor_generation:1,preparation_key:"b".repeat(64)},
    result:{protocol_version:"1.0",job_id:jobId,operation:"fixture.echo",schema_version:"1.2.0",
      input_hash:"a".repeat(64),environment_fingerprint:"c".repeat(64),worker_id:"worker-a",
      inline_result:{value:3},artifact:null,failure:null,warnings:[]}};
  let lose=true, badResult=false, badReceipt=false, foreign=false, reads=0, badArtifact=false;
  const submissions=[],cancellations=[];
  const fixture=environment(async(url,options)=>{
    if(url.endsWith("/control"))return json(initialControl());
    if(url.endsWith("/assignments"))return json([assigned()]);
    if(url.endsWith("/artifact")){
      assert.equal(url,"/runtime/api/jobs/"+jobId+"/artifact");
      assert.equal(options.method,"GET");
      return json(badArtifact?[]:{values:[1,2,3]});
    }
    if(options.method==="POST"){
      const body=JSON.parse(options.body);
      if(url.endsWith("/cancel")){
        cancellations.push(body);receipt.cancel_requested=true;
        return json(receipt,202);
      }
      submissions.push(body);
      receipt.request_id=body.request_id;
      if(lose){lose=false;throw new Error("lost after durable receipt");}
      return json({...receipt,operation:badReceipt?"foreign.operation":receipt.operation},202);
    }
    const value={...receipt,run_id:foreign?webcrypto.randomUUID():run};
    if(url.endsWith("/result")){
      reads++;
      const result=structuredClone(outcome);
      if(badResult)result.execution.executor_generation++;
      return json({schema_version:1,job:value,result});
    }
    return json(url.endsWith("/jobs")?[value]:value);
  });
  const client=new fixture.api.RuntimeClient(run);
  await client.refresh();
  const requestId=webcrypto.randomUUID();
  await assert.rejects(()=>client.submitJob(lease,"fixture.echo",{value:3},{requestId}),
    error=>error.uncertain&&error.requestId===requestId);
  assert.equal(submissions.length,1);
  const saved=await client.submitJob(lease,"fixture.echo",{value:3},{requestId});
  assert.equal(saved.id,jobId);assert.equal(submissions.length,2);
  assert.equal(submissions[0].request_id,submissions[1].request_id);
  assert.deepEqual(Object.keys(submissions[0]).sort(),["operation","parameters","request_id"]);
  assert.equal((await client.listJobs())[0].id,jobId);
  assert.equal((await client.job(jobId)).state,"queued");
  assert.equal((await client.jobResult(jobId)).result.result.inline_result.value,3);
  assert.deepEqual(JSON.parse(JSON.stringify(await client.jobArtifact(jobId))),{values:[1,2,3]});
  badArtifact=true;await assert.rejects(()=>client.jobArtifact(jobId));badArtifact=false;
  badResult=true;await assert.rejects(()=>client.jobResult(jobId));badResult=false;
  foreign=true;await assert.rejects(()=>client.job(jobId));await assert.rejects(()=>client.listJobs());foreign=false;
  const cancelId=webcrypto.randomUUID();
  const canceled=await client.cancelJob(jobId,{requestId:cancelId});
  assert.equal(canceled.state,"queued");assert.equal(canceled.cancel_requested,true);
  assert.equal(canceled.cancel_acknowledged,false,"requested must not pretend acknowledged or terminated");
  assert.deepEqual(cancellations,[{request_id:cancelId}]);
  receipt.state="revoked";receipt.current_assignment=false;
  assert.equal((await client.jobResult(jobId)).job.current_assignment,false);
  assert.equal(submissions.length,2,"historical reads cannot submit another operation");
  badReceipt=true;
  await assert.rejects(()=>client.submitJob(lease,"fixture.echo",{value:3},{requestId}),
    error=>error.uncertain&&error.requestId===requestId);
  badReceipt=false;
  for(const parameters of [[],null,{value:NaN},{value:Infinity},{value:()=>3},{value:undefined},
      {value:new Date()},{large:"x".repeat(65536)}]){
    assert.throws(()=>client.submitJob(lease,"fixture.echo",parameters));
  }
  const getter={get value(){throw new Error("must not invoke an input getter");}};
  assert.throws(()=>client.submitJob(lease,"fixture.echo",getter),/passive/);
  assert.throws(()=>client.submitJob(lease,"eval",{}));
  assert.throws(()=>client.job("../private"));
  assert.throws(()=>client.cancelJob("../private"));
  assert.throws(()=>client.jobArtifact("https://foreign.test/private"));
  assert.ok(reads>=3);
  client.close();assert.equal(client.controllers.size,0);
}
// A shared control panel must not poll private terminal roles as scientific executors.
{
  const terminalLease=webcrypto.randomUUID(), calls=[];
  const control={...initialControl(),preparation_control:true,profiles:[
    {id:"line-parameters",kind:"scientific"},{id:"julia-terminal",kind:"terminal"}]};
  const owned=[{...assigned(),usable:true},{...assigned(),id:terminalLease,role:"terminal",profile:"julia-terminal",usable:true}];
  const fixture=environment(async url=>{
    calls.push(url);
    return json(url.endsWith("/control")?control:url.endsWith("/assignments")?owned:{});
  });
  const client=fixture.api.acquire(run);
  try {
    await client.refresh();await client.refreshScience();
    assert.equal(calls.some(url=>url.endsWith("/assignments/"+terminalLease+"/science")),false);
    assert.equal(calls.some(url=>url.endsWith("/assignments/"+lease+"/science")),true);
    client.publishScience({[terminalLease]:{preparation:"cold"}});
    assert.equal(Object.hasOwn(client.state.science,terminalLease),false);
  } finally {client.close();}
}
console.log("Runtime browser client: passive ownership, coalescing, stale state, independent logs/preparation, terminal-role separation, expiring evidence, job provenance, mutation fencing and teardown passed.");
