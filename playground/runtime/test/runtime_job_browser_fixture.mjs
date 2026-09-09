// Explicit HTTP fixture for the real Bonito browser binding test. Production
// gateway/worker/artifact authorization is covered by the separate TLS suite.
export function installJobFixture() {
  const samples=JSON.parse(document.querySelector('#job-fixture-data').dataset.jobFixture);
  const seed=samples['3'], fence=seed.fence, execution=seed.execution;
  const api=globalThis.LineCableModelsRuntimeClient;
  const client=api.acquire(fence.run_id), original=client.fetcher;
  const lease={id:fence.lease_id,run_id:fence.run_id,role:fence.role,profile:fence.profile_id,
    worker_id:fence.worker_id,worker_boot:fence.worker_boot,generation:fence.generation,usable:true,state:'active'};
  const inventory={schema_version:1,enabled:true,preparation_control:true,assigned_execution:true,broker:'online',
    workers:[],provisioned:[],profiles:[{id:'fixture',operations:['system.echo']}]};
  const report={channel:'online',phase:'idle',preparation:'ready',valid_for_ms:5000,pending:false,accepted:true,
    progress:0,elapsed_seconds:0,output_lines:0,current_request_id:null,...execution,failure:null};
  const fixture={submissions:[],rows:new Map(),loseSubmit:false,inventory,report};
  const json=value=>new Response(JSON.stringify(value),{headers:{'Content-Type':'application/json'}});
  client.fetcher=async (url,options) => {
    if (!url.startsWith('/runtime/api/')) throw Error('Fixture URL escaped the runtime API');
    if (url.endsWith('/control')) return json(inventory);
    if (url.endsWith('/assignments')) return json([lease]);
    if (url.endsWith('/science')) return json(report);
    if (url.endsWith('/jobs') && options.method==='POST') {
      const body=JSON.parse(options.body);fixture.submissions.push(body);
      if (!fixture.rows.has(body.request_id)) {
        const envelope=structuredClone(samples[String(body.parameters.value)]);
        if (!envelope) throw Error('No scientific fixture for these inputs');
        const receipt={...lease,id:crypto.randomUUID(),request_id:body.request_id,lease_id:lease.id,
          operation:body.operation,input_hash:envelope.result.input_hash,...execution,
          submitted_at:'2026-09-08T00:00:00.000Z',deadline:'2026-09-08T00:01:00.000Z',state:'queued',
          current_assignment:true,channel:'online',cancel_requested:false,cancel_acknowledged:false};
        envelope.result.job_id=receipt.id;
        fixture.rows.set(body.request_id,{receipt,envelope});
      }
      if (fixture.loseSubmit) {fixture.loseSubmit=false;throw Error('Fixture lost acknowledgement');}
      return json(fixture.rows.get(body.request_id).receipt);
    }
    const row=[...fixture.rows.values()].find(row=>url.includes(row.receipt.id));
    if (!row) throw Error('Unknown fixture job');
    if (!url.endsWith('/result')) return json(row.receipt);
    return json({schema_version:1,job:row.receipt,result:row.receipt.state==='succeeded'?row.envelope:null});
  };
  fixture.job=globalThis.LineCableModelsRuntimeControls.mount(document.querySelector('[data-runtime-kind="execution"]')).job;
  fixture.complete=async () => {
    fixture.rows.get(fixture.job.active.requestId).receipt.state='succeeded';
    if (fixture.job.refreshing) await fixture.job.refreshing;
    await fixture.job.refresh();
  };
  fixture.ready=async () => {await client.refresh();await client.refreshScience();};
  fixture.restore=()=>{client.fetcher=original;};
  globalThis.__scientificFixture=fixture;
  return fixture.ready();
}
