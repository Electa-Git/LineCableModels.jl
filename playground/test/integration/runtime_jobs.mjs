import assert from "node:assert/strict";
import vm from "node:vm";
import {readFile} from "node:fs/promises";
import {webcrypto, createHash} from "node:crypto";
import {setTimeout as delay} from "node:timers/promises";

const script = await readFile(new URL("../../assets/runtime-client.js",import.meta.url),"utf8");
const uuid = () => webcrypto.randomUUID();
const json = value => new Response(JSON.stringify(value),{headers:{"Content-Type":"application/json"}});
async function until(check) {
  for (let i=0;i<400;i++) { if (check()) return; await delay(5); }
  assert.fail("Job view did not reach its finite expected state");
}
function fixture() {
  const run = uuid();
  const lease = {id:uuid(),run_id:run,role:"main",profile:"fixture",worker_id:"worker-a",
    worker_boot:uuid(),generation:1,usable:true,state:"active"};
  const report = {channel:"online",phase:"idle",preparation:"ready",valid_for_ms:5000,pending:false,accepted:true,
    progress:0,elapsed_seconds:0,output_lines:0,current_request_id:null,executor_id:uuid(),executor_generation:1,
    preparation_key:"b".repeat(64),failure:null};
  const inventory = {schema_version:1,enabled:true,preparation_control:true,assigned_execution:true,broker:"online",
    workers:[],provisioned:[],profiles:[{id:"fixture",operations:["fixture.echo"]}]};
  const rows = new Map(), requests = [], submissions = [], cancellations = [];
  const knobs = {loseSubmit:false,loseCancel:false,artifact:false,gate:null,artifactWaiting:false,offline:false,resultOffline:false};
  const fetcher = async (url, options) => {
    requests.push([url,options.method]);
    if (knobs.offline) throw new Error("fixture offline");
    if (url.endsWith("/control")) return json(inventory);
    if (url.endsWith("/assignments")) return json([lease]);
    if (url.endsWith("/science")) return json(report);
    if (options.method === "POST") {
      const body = JSON.parse(options.body);
      if (url.endsWith("/jobs")) {
        submissions.push(body);
        if (!rows.has(body.request_id)) {
          const receipt = {id:uuid(),request_id:body.request_id,run_id:run,lease_id:lease.id,role:lease.role,
            operation:body.operation,input_hash:createHash("sha256").update(JSON.stringify(body.parameters)).digest("hex"),
            worker_id:lease.worker_id,worker_boot:lease.worker_boot,generation:lease.generation,
            executor_id:report.executor_id,executor_generation:report.executor_generation,preparation_key:report.preparation_key,
            submitted_at:"2026-09-08T00:00:00.000Z",deadline:"2026-09-08T00:01:00.000Z",state:"queued",
            current_assignment:true,channel:"online",cancel_requested:false,cancel_acknowledged:false};
          rows.set(body.request_id,{receipt,parameters:body.parameters,artifact:knobs.artifact});
        }
        if (knobs.loseSubmit) { knobs.loseSubmit=false; throw new Error("lost after saving receipt"); }
        return json(rows.get(body.request_id).receipt);
      }
      cancellations.push({url,body});
      const row = [...rows.values()].find(r => url.includes(r.receipt.id));
      row.receipt.cancel_requested=true;
      if (knobs.loseCancel) {knobs.loseCancel=false;throw new Error("lost cancellation reply");}
      row.receipt.cancel_acknowledged=true;
      return json(row.receipt);
    }
    const row = [...rows.values()].find(r => url.includes(r.receipt.id));
    assert.ok(row,"unknown job URL");
    if (url.endsWith("/artifact")) {
      knobs.artifactWaiting=true;
      if (knobs.gate) await knobs.gate;
      return json(row.parameters);
    }
    const receipt = row.receipt;
    if (!url.endsWith("/result")) return json(receipt);
    if (knobs.resultOffline) throw new Error("fixture result channel unavailable");
    const result = ["succeeded","failed","canceled"].includes(receipt.state) ? {
      protocol_version:"2.0",fence:{...receipt,fingerprint:"c".repeat(64)},
      execution:{executor_id:receipt.executor_id,executor_generation:receipt.executor_generation,preparation_key:receipt.preparation_key},
      result:{protocol_version:"1.0",job_id:receipt.id,operation:receipt.operation,input_hash:receipt.input_hash,
        worker_id:receipt.worker_id,environment_fingerprint:"c".repeat(64),schema_version:"fixture.v1",warnings:[],
        inline_result:receipt.state === "succeeded" && !row.artifact ? row.parameters : null,
        artifact:receipt.state === "succeeded" && row.artifact ? {retrieval_reference:"https://must-not-follow.invalid/"} : null,
        failure:receipt.state !== "succeeded" ? {category:receipt.state} : null}
    } : null;
    return json({schema_version:1,job:receipt,result});
  };
  const document = new EventTarget(); document.visibilityState="visible";
  const context = vm.createContext({fetch:fetcher,document,crypto:webcrypto,AbortController,Response,
    TextEncoder,TextDecoder,Uint8Array,setTimeout,clearTimeout});
  vm.runInContext(script,context);
  const api = context.LineCableModelsRuntimeClient;
  const client = new api.RuntimeClient(run,{pollMs:50});
  const job = new api.RuntimeJob(client,{role:"main",operation:"fixture.echo",parameters:{value:1}});
  const row = () => rows.get(job.active?.requestId);
  async function ready() {await client.refresh();await client.refreshScience();await until(()=>job.state.canRun);}
  async function complete(state="succeeded") {
    row().receipt.state=state;
    await until(()=>!job.refreshing);
    await job.refresh();
  }
  return {client,job,lease,report,inventory,rows,requests,submissions,cancellations,knobs,row,ready,complete};
}

const f = fixture();
const {client,job,knobs,submissions,cancellations} = f;
assert.equal(f.requests.length,0,"construction must be inert");
const detach = job.subscribe(()=>{});
try {
  await f.ready();
  assert.equal(submissions.length,0,"status and mount never run a calculation");
  await job.run();await f.complete();await until(()=>job.state.current);
  assert.equal(job.state.lastGood.value.value,1);
  const first = job.state.lastGood;
  job.setInputs({value:1});assert.equal(job.state.current,true,"identical input echoes are not changes");
  job.markInputsPending();
  assert.equal(job.state.inputsPending,true);
  assert.equal(job.state.current,false,"a local edit immediately invalidates current data");
  assert.equal(job.state.canRun,false,"Run cannot overtake an unacknowledged field edit");
  const beforeEdit=submissions.length;
  await assert.rejects(()=>job.run(),/Select and prepare/);
  assert.equal(submissions.length,beforeEdit,"pending input never submits the previous value");
  job.setInputs({value:1});
  assert.equal(job.state.inputsPending,false,"an acknowledged same-value edit also clears the barrier");
  assert.equal(job.state.canRun,true);
  assert.equal(job.state.lastGood,first);
  job.setInputs({value:2});assert.equal(job.state.current,false);assert.equal(job.state.lastGood,first);

  // Changing inputs while a private artifact is in flight cannot replace data.
  knobs.artifact=true;
  let release;knobs.gate=new Promise(resolve=>release=resolve);
  await f.ready();await job.run();f.row().receipt.state="succeeded";
  await until(()=>!job.refreshing);const pending=job.refresh();
  await until(()=>knobs.artifactWaiting);
  job.setInputs({value:3});release();await pending;
  assert.equal(job.state.superseded,true);assert.equal(job.state.lastGood,first);
  assert.ok(f.requests.every(([url])=>url.startsWith("/runtime/api/")),"never follow a supplied artifact URL");
  knobs.gate=null;knobs.artifact=false;knobs.artifactWaiting=false;
  await f.ready();await job.run();await f.complete();await until(()=>job.state.current);
  const third=job.state.lastGood;
  assert.equal(third.value.value,3);
  assert.throws(()=>job.setInputs({value:NaN}));
  assert.equal(job.state.current,false);assert.equal(job.state.canRun,false);
  job.setInputs({value:3});await f.ready();

  // Explicit retry preserves the original input and UUID, even after a new draft.
  knobs.loseSubmit=true;
  await job.run();assert.equal(job.state.phase,"unconfirmed");
  const count=submissions.length, original=submissions.at(-1);
  job.setInputs({value:4});await delay(100);
  assert.equal(submissions.length,count,"unconfirmed input must never replay automatically");
  assert.equal(job.state.canRun,false);
  await job.retry();assert.deepEqual(submissions.at(-1),original);
  await f.complete();assert.equal(job.state.lastGood,third);assert.equal(job.state.superseded,true);

  // Cancellation is an acknowledgement, not a fabricated canceled result.
  await f.ready();await job.run();knobs.loseCancel=true;
  await job.cancel();assert.equal(job.pending.kind,"cancel");
  await job.retry();assert.deepEqual(cancellations[0],cancellations[1]);
  assert.equal(job.state.receipt.state,"queued");assert.equal(job.state.receipt.cancel_acknowledged,true);
  await f.complete("canceled");assert.equal(job.state.phase,"canceled");assert.equal(job.state.lastGood,third);

  // A newly prepared process cannot inherit the previous process's completion.
  await f.ready();await job.run();
  f.report.executor_id=uuid();f.report.executor_generation++;
  await client.refresh();await client.refreshScience();
  await f.complete();assert.equal(job.state.superseded,true);assert.equal(job.state.lastGood,third);

  // Frozen/expired readiness disables starts and retains old data honestly.
  knobs.offline=true;await client.refresh();
  assert.equal(job.state.canRun,false);assert.equal(job.state.current,false);assert.equal(job.state.lastGood,third);
  knobs.offline=false;await f.ready();
  f.inventory.assigned_execution=false;await client.refresh();assert.equal(job.state.canRun,false);
  f.inventory.assigned_execution=true;await f.ready();
  f.inventory.broker="unavailable";await client.refresh();assert.equal(job.state.canRun,false);
  f.inventory.broker="online";await f.ready();

  // Readiness may expire between publishing button state and capturing intent.
  const context=job.context.bind(job);let checks=0;
  job.context=()=>++checks===2 ? null : context();
  const beforeExpiry=submissions.length;
  await assert.rejects(()=>job.run(),/Select and prepare/);
  assert.equal(submissions.length,beforeExpiry);
  job.context=context;

  // Teardown must not close a sibling component's shared inventory connection.
  const sibling=client.subscribe(()=>{});
  const before=submissions.length;detach();
  assert.equal(job.closed,true);assert.equal(job.state,null);assert.equal(client.closed,false);
  await assert.rejects(()=>job.run(),/closed/);
  await delay(80);assert.equal(submissions.length,before);
  sibling();assert.equal(client.closed,true);assert.equal(client.controllers.size,0);
} finally {job.close();client.close();}
const revoked=fixture();revoked.job.subscribe(()=>{});
try {
  await revoked.ready();await revoked.job.run();await revoked.complete();
  assert.equal(revoked.job.state.current,true);
  const retained=revoked.job.state.lastGood;
  revoked.row().receipt.current_assignment=false;await revoked.job.refresh();
  assert.equal(revoked.job.state.current,false,"explicit negative assignment evidence overrides cached readiness");
  revoked.row().receipt.current_assignment=true;await revoked.job.refresh();
  assert.equal(revoked.job.state.current,false,"an older status cannot revive a revoked display");
  assert.equal(revoked.job.state.lastGood,retained);
  await revoked.ready();await revoked.job.run();await revoked.complete();
  assert.equal(revoked.job.state.current,true,"a new explicit successful calculation can replace historical data");
} finally {revoked.job.close();revoked.client.close();}
console.log("Scientific job view: explicit run/retry/cancel, private artifact projection, input/executor/revocation fencing, last-good retention and shared-client teardown passed.");

const independent=fixture();independent.job.subscribe(()=>{});
try {
  await independent.ready();await independent.job.run();await independent.complete();
  const retained=independent.job.state.lastGood;
  await independent.ready();await independent.job.run();
  independent.knobs.resultOffline=true;
  const before=independent.requests.filter(([url])=>url.endsWith('/result')).length;
  await independent.complete('revoked');
  assert.equal(independent.job.state.phase,'revoked');
  assert.equal(independent.job.state.lastGood,retained);
  assert.equal(independent.job.state.error,null);
  assert.equal(independent.requests.filter(([url])=>url.endsWith('/result')).length,before,
    'Revoked receipt must remain visible without a broker result');
  independent.knobs.resultOffline=false;
  await independent.ready();await independent.job.run();
  independent.knobs.resultOffline=true;
  await independent.complete();
  assert.equal(independent.job.state.phase,'succeeded');
  assert.equal(independent.job.state.lastGood,retained);
  independent.knobs.resultOffline=false;
  await until(()=>independent.job.state.current && independent.job.state.lastGood!==retained);
  assert.notEqual(independent.job.state.lastGood,retained,'Delayed successful data must still recover');
} finally {independent.job.close();independent.client.close();}
console.log('Durable job lifecycle remains visible during result-channel loss; delayed successful data recovers.');
