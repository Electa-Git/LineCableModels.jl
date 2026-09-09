// Actual registered Bonito/Reveal consumers, actual runtime HTTP, no fetch/state fixtures.
import assert from 'node:assert/strict';
import {spawn} from 'node:child_process';
import {readFile,writeFile,access} from 'node:fs/promises';
import {join} from 'node:path';
import {setTimeout as delay} from 'node:timers/promises';

const [base,directory] = process.argv.slice(2);
assert.equal(new URL(base).hostname,'127.0.0.1');
let browser,socket,sequence=0,checks=0,chromeLog='';
const pending=new Map(),errors=[],traffic=[],evidence={deck:{},workbench:{},timings:{}};
// A recovery failure must not discard already measured numerical results.
// This checkpoint is deliberately not the successful acceptance artifact.
const checkpoint=()=>writeFile(join(directory,'scientific-partial-results.json'),JSON.stringify(evidence,null,2));
const stop=()=>{
  for(const item of pending.values()){clearTimeout(item.timer);item.reject(Error('Browser fixture stopped'));}
  pending.clear();socket?.close();
  if(browser && browser.exitCode===null)browser.kill('SIGTERM');
};
process.once('SIGTERM',()=>{process.exitCode=1;stop();});
const command=(method,params={})=>new Promise((resolve,reject)=>{
  const id=++sequence,timer=setTimeout(()=>{pending.delete(id);reject(Error('CDP deadline: '+method));},30000);
  pending.set(id,{resolve,reject,timer});socket.send(JSON.stringify({id,method,params}));
});
const evaluate=async expression=>{
  const result=await command('Runtime.evaluate',{expression,returnByValue:true,awaitPromise:true});
  if(result.exceptionDetails)throw Error(result.exceptionDetails.exception?.description ?? JSON.stringify(result.exceptionDetails));
  return result.result.value;
};
const wait=async(expression,label,seconds=90)=>{
  for(const end=Date.now()+seconds*1000;Date.now()<end;){
    if(await evaluate(expression))return;
    await delay(150);
  }
  throw Error(label);
};
const check=async(expression,label)=>{assert(await evaluate(expression),label);checks++;};
const screenshot=async name=>{
  const shot=await command('Page.captureScreenshot',{format:'png'});
  await writeFile(join(directory,name+'.png'),Buffer.from(shot.data,'base64'));
};
const button=(kind,label)=>`btn(view.querySelector('[data-runtime-kind=${kind}]'),${JSON.stringify(label)})`;
const click=async(kind,label)=>{
  const expression=button(kind,label);
  await wait(`${expression} && !${expression}.disabled`,'Action did not become available: '+label);
  await evaluate(`${expression}.click()`);
};
const current=()=>wait("job.state.current && view.querySelector('[data-runtime-kind=execution] [data-result-current]').dataset.resultCurrent==='true' && view.querySelector('[data-runtime-kind=execution]').dataset.resultProjection==='accepted' && !!view.querySelector('.lc-study-curve').getAttribute('points')",'Successful data was not projected',30);
const runCalculation=async()=>{
  const previous=await evaluate('job.state.lastGood?.receipt.id ?? null');
  await click('execution','Run calculation');
  await wait('(()=>{if(["revoked","failed","canceled","rejected","uncertain"].includes(job.state.phase))throw Error("Calculation ended: "+job.state.phase+" · "+job.state.error);return job.state.lastGood?.receipt.id && job.state.lastGood.receipt.id!=='+JSON.stringify(previous)+';})()','Calculation did not produce a new successful result',150);
  await current();
};
async function bind(consumer,name){
  if(consumer==='deck'){
    await evaluate("studyFrame.style.visibility='hidden';Reveal.slide(Reveal.getSlides().findIndex(s=>s.querySelector('iframe[data-lcm-src=\"/science/"+(name==='parameters'?'line-parameters':'corridor')+"\"]')))");
    const selector='iframe[data-lcm-src="/science/'+(name==='parameters'?'line-parameters':'corridor')+'"]';
    await wait('document.querySelector('+JSON.stringify(selector)+")?.dataset.lcmState==='ready'",'Scientific deck frame did not mount');
    await evaluate('window.doc=document.querySelector('+JSON.stringify(selector)+').contentDocument;window.view=doc.querySelector(".lc-study-view")');
  }else{
    await evaluate("studyFrame.style.visibility='visible'");
    await evaluate("window.doc=studyFrame.contentDocument;[...doc.querySelectorAll('.lc-wb-nav-item')].find(b=>b.textContent.includes("+JSON.stringify(name==='parameters'?'Line parameters':'OHL / UGC')+")).click()");
    await wait("doc.querySelector('[data-view="+name+"]').classList.contains('is-active')",'Workbench view did not activate');
    await evaluate("window.view=doc.querySelector('[data-view="+name+"] .lc-study-view')");
  }
  await wait("doc.defaultView.WEBSOCKET?.isopen() && !!view?.querySelector('[data-runtime-kind=execution]') && typeof view.querySelector('[name=frequency_points]')?.oninput==='function'",'Scientific control bindings did not initialize');
  await evaluate("window.job=doc.defaultView.LineCableModelsRuntimeControls.mount(view.querySelector('[data-runtime-kind=execution]')).job;window.plot=view.querySelector('.lc-study-plot'); for(const d of view.querySelectorAll('details'))d.open=true");
}
async function calculate(consumer,name){
  console.log('Scientific browser: '+consumer+' / '+name+' cold preparation');
  await bind(consumer,name);
  await check('!job.state.lastGood && !job.state.canRun','Fresh view falsely implies a calculation or ready worker');
  await evaluate("change('minimum_frequency_hz','100');change('maximum_frequency_hz','150');change('frequency_points','2')");
  await wait("job.inputsValid && (job.parameters.frequency_points===2 || job.parameters.frequencies_hz?.length===2)",'Browser inputs did not reach the shared job');
  const profile=name==='parameters'?'line-parameters':'power-flow', worker=name==='parameters'?'worker-a':'worker-b';
  await wait("[...view.querySelector('[data-runtime-kind=selector] select').options].some(o=>o.value==="+JSON.stringify(profile)+")",'Registered profile unavailable');
  await evaluate("(()=>{const [profile,placement,worker]=view.querySelectorAll('[data-runtime-kind=selector] select');setSelect(profile,"+JSON.stringify(profile)+");setSelect(placement,'pinned');setSelect(worker,"+JSON.stringify(worker)+");})()");
  await click('selector','Assign worker');
  await wait(`(() => {
    const assignment = job.assignment();
    if (assignment && ['released','releasing','reconciling','failed','expired'].includes(assignment.state))
      throw Error('Cold assignment lost before preparation: ' + assignment.state + ' / revision ' + assignment.revision);
    return view.querySelector('[data-preparation]')?.dataset.preparation === 'cold';
  })()`, 'Acknowledged assignment was not cold');
  await check('!job.state.canRun && !job.state.lastGood','Assignment alone allowed execution');
  const started=Date.now();
  await click('preparation','Prepare executor');
  // UI navigation/health must not wait for preparation to finish.
  assert.equal((await fetch(base+'/health')).status,200);checks++;
  await check("document.querySelector('.lcm-deck-status')!==null",'Presentation shell disappeared during preparation');
  await wait('(()=>{if(["released","reconciling","releasing"].includes(job.assignment()?.state))throw Error("Assignment lost during preparation");const status=job.client.state.science[job.assignment()?.id];if(status?.phase==="failed")throw Error("Preparation failed: "+status.failure);return job.state.canRun;})()',
    'Explicit cold preparation did not become ready',650);
  evidence.timings[consumer+'-'+name+'-cold-seconds']=(Date.now()-started)/1000;
  await check('!job.state.lastGood','Preparation implicitly submitted the scientific calculation');
  const warmExecutor=await evaluate('job.context().executor_id');
  const warmStarted=Date.now();
  await click('preparation','Prepare executor');
  await wait('job.state.canRun','Repeated explicit preparation did not recover',30);
  evidence.timings[consumer+'-'+name+'-warm-seconds']=(Date.now()-warmStarted)/1000;
  await check('job.context().executor_id==='+JSON.stringify(warmExecutor),'Warm preparation replaced the executor');
  await evaluate("for(const d of view.querySelectorAll('details'))d.open=false");
  await runCalculation();
  evidence[consumer][name]=await evaluate('({parameters:job.parameters,receipt:job.state.lastGood.receipt,value:job.state.lastGood.value,provenance:job.state.lastGood.provenance})');
  await checkpoint();
  await check("plot===view.querySelector('.lc-study-plot') && [...plot.querySelectorAll('polyline')].filter(p=>p.getAttribute('points')).length>0",'Numerical result did not update the persistent SVG');
  await check("view.textContent.includes('Current result') && view.querySelector('[aria-label=\"Result provenance\"]').textContent.includes(job.state.lastGood.receipt.input_hash)",'Visible result provenance missing');
  if(consumer==='workbench' && name==='parameters'){
    await evaluate("window.savedGood=job.state.lastGood;window.savedCurve=view.querySelector('.lc-study-curve').getAttribute('points');change('separation_m','0.8')");
    await wait("!job.state.current && view.textContent.includes('Outdated result')",'Edited inputs did not mark last-good data outdated');
    await check("job.state.lastGood===savedGood && view.querySelector('.lc-study-curve').getAttribute('points')===savedCurve",'Editing erased the last successful plot');
    await check("(()=>{change('separation_m','0.5');return job.state.inputsPending && !job.state.canRun;})()",
      'Run overtook the latest visible input edit');
    await runCalculation();
  }
  await screenshot(consumer+'-'+name+'-calculated');
  console.log('Scientific browser: '+consumer+' / '+name+' real result passed');
}
async function cancelAndRecover(){
  console.log('Scientific browser: cancel cold corridor preparation and explicitly retry');
  await bind('workbench','corridor');
  await evaluate("window.beforeCancel=job.state.lastGood;window.beforeCancelPlot=view.querySelector('.lc-study-curve').getAttribute('points');window.beforeCancelLease=job.assignment().id;for(const d of view.querySelectorAll('details'))d.open=true");
  // Warm scientific jobs can finish before a status poll. Exercise genuinely
  // cold preparation; the separate TLS suite covers running-job cancellation.
  await click('selector','Release assignment');
  await wait('job.assignment()?.state==="released"','Owned assignment was not released',30);
  await click('selector','Assign worker');
  await wait('job.assignment()?.id!==beforeCancelLease && view.querySelector("[data-preparation]")?.dataset.preparation==="cold"',
    'Replacement assignment was not acknowledged cold',30);
  await click('preparation','Prepare executor');
  await wait('job.client.state.science[job.assignment()?.id]?.phase==="preparing"',
    'Cold preparation did not reach a cancellable phase',60);
  await evaluate('window.beforeCancelExecutor=job.client.state.science[job.assignment().id].executor_id');
  await click('preparation','Cancel preparation');
  await wait('job.client.state.science[job.assignment()?.id]?.failure==="canceled"',
    'Worker did not confirm preparation cancellation',30);
  await check('!job.state.current && job.state.lastGood===beforeCancel && view.querySelector(".lc-study-curve").getAttribute("points")===beforeCancelPlot',
    'Cancellation erased the previous successful result');
  await wait('!job.state.canRun','Canceled executor still advertised ready',10);
  await check('job.assignment().usable','Canceling preparation incorrectly released its assignment');
  await evaluate('window.beforeRecoveryRequest=job.client.state.science[job.assignment().id].request_id');
  const started=Date.now();
  await click('preparation','Prepare executor');
  await wait('(()=>{if(["released","reconciling","releasing"].includes(job.assignment()?.state))throw Error("Assignment lost during recovery");const status=job.client.state.science[job.assignment()?.id];if(status?.phase==="failed" && status.request_id!==beforeRecoveryRequest)throw Error("Recovery preparation failed: "+status.failure);return job.state.canRun;})()',
    'Explicit post-cancellation preparation did not become ready',650);
  evidence.timings['workbench-corridor-recovery-seconds']=(Date.now()-started)/1000;
  await checkpoint();
  await check('job.context().executor_id!==beforeCancelExecutor','Cancellation recovery reused the canceled process identity');
  await evaluate("for(const d of view.querySelectorAll('details'))d.open=false");
  await runCalculation();
  await check('job.state.lastGood!==beforeCancel && job.state.current','Explicit retry did not publish a new current result');
  await screenshot('corridor-cancellation-recovered');
}
try{
  browser=spawn(process.env.LCM_BROWSER||'google-chrome',['--headless=new','--disable-gpu','--no-first-run','--no-default-browser-check',
    '--remote-debugging-address=127.0.0.1','--remote-debugging-port=0','--user-data-dir='+join(directory,'chrome'),'about:blank'],{stdio:['ignore','ignore','pipe']});
  browser.stderr.on('data',chunk=>{chromeLog=(chromeLog+chunk.toString()).slice(-65536);});
  browser.on('error',error=>errors.push(String(error)));
  let debugPort;
  for(let i=0;i<150;i++){
    try{debugPort=(await readFile(join(directory,'chrome','DevToolsActivePort'),'utf8')).split('\n')[0];break;}catch{}
    if(browser.exitCode!==null)throw Error('Owned Chrome exited: '+chromeLog);
    await delay(100);
  }
  assert(debugPort,'Owned Chrome did not start');
  const pages=await(await fetch('http://127.0.0.1:'+debugPort+'/json')).json();
  socket=new WebSocket(pages.find(page=>page.type==='page').webSocketDebuggerUrl);
  await new Promise((resolve,reject)=>{socket.addEventListener('open',resolve,{once:true});socket.addEventListener('error',reject,{once:true});});
  socket.addEventListener('message',event=>{
    const message=JSON.parse(event.data);
    if(message.method==='Runtime.exceptionThrown')errors.push(message.params.exceptionDetails);
    if(message.method==='Network.responseReceived' && message.params.response.url.includes('/runtime/api/')){
      traffic.push({url:message.params.response.url,status:message.params.response.status});
      if(traffic.length>200)traffic.shift();
    }
    const item=pending.get(message.id);if(!item)return;
    pending.delete(message.id);clearTimeout(item.timer);
    message.error?item.reject(Error(message.error.message)):item.resolve(message.result);
  });
  await command('Page.enable');await command('Runtime.enable');await command('Network.enable');
  await command('Emulation.setDeviceMetricsOverride',{width:1440,height:900,deviceScaleFactor:1,mobile:false});
  let runs;
  for(const end=Date.now()+180000;Date.now()<end;){
    runs=await(await fetch(base+'/runtime/api/runs')).json();
    assert(!runs.some(run=>run.state==='failed'),JSON.stringify(runs));
    if(runs.length===2 && runs.every(run=>run.state==='running'))break;
    await delay(250);
  }
  assert(runs.length===2 && runs.every(run=>run.state==='running'),'Registered UI hosts did not start');
  const deck=runs.find(run=>run.application==='ichqp-showcase'),study=runs.find(run=>run.application==='cable-study');
  await command('Page.navigate',{url:base+'/presentations/showcase.html?lcm-run='+deck.id});
  await wait('window.Reveal?.isReady()','Reveal did not initialize');
  await evaluate("window.btn=(root,label)=>[...root.querySelectorAll('button')].find(b=>b.textContent===label);window.setSelect=(input,value)=>{input.value=value;input.dispatchEvent(new input.ownerDocument.defaultView.Event('change',{bubbles:true}));};window.change=(name,value)=>{const input=view.querySelector('[name='+name+']');input.value=value;input.dispatchEvent(new doc.defaultView.Event('input',{bubbles:true}));}");
  await evaluate("window.studyFrame=document.createElement('iframe');studyFrame.id='study-frame';studyFrame.src="+JSON.stringify(base+'/applications/runs/'+study.id+'/workbenches/cable-study')+";studyFrame.style.cssText='position:fixed;inset:0;width:100vw;height:100vh;border:0;z-index:10000';document.body.append(studyFrame)");
  await wait("studyFrame.contentWindow.WEBSOCKET?.isopen() && typeof [...studyFrame.contentDocument.querySelectorAll('.lc-wb-nav-item')].find(b=>b.textContent.includes('Line parameters'))?.onclick==='function'",'Workbench actions did not initialize');
  await calculate('deck','parameters');
  await calculate('deck','corridor');
  await calculate('workbench','parameters');
  await calculate('workbench','corridor');
  await cancelAndRecover();
  await evaluate('window.retained=job.state.lastGood;window.retainedPlot=view.querySelector(".lc-study-curve").getAttribute("points")');
  await writeFile(join(directory,'stop-power-worker'),'stop only the owned worker-b fixture\n');
  await wait("!job.state.current && !job.state.canRun && job.client.state.control.workers.find(w=>w.registration.worker_id==='worker-b')?.liveness!=='online'",'Worker loss did not invalidate scientific readiness',35);
  await check("job.state.lastGood===retained && view.querySelector('.lc-study-curve').getAttribute('points')===retainedPlot && view.textContent.includes('Outdated result')",'Worker loss erased retained scientific data');
  assert.equal((await fetch(base+'/health')).status,200);checks++;
  await bind('workbench','parameters');
  await runCalculation();
  await check('job.state.current','Power worker loss contaminated the line-parameter role');
  await screenshot('line-survives-power-worker-loss');
  if(process.env.LCM_PHYSICAL_ARTIFACT_CHECK==='1'){
    console.log('Scientific browser: verify a real cross-host private S3 result');
    await evaluate("change('frequency_points','200')");
    await wait('job.parameters.frequencies_hz?.length===200 && job.inputsValid','Large-result inputs were not accepted');
    await runCalculation();
    const item=await evaluate('({parameters:job.parameters,receipt:job.state.lastGood.receipt,value:job.state.lastGood.value})');
    const result=await(await fetch(base+'/runtime/api/jobs/'+item.receipt.id+'/result')).json();
    assert.equal(result.result.result.inline_result,null);checks++;
    assert.equal(result.result.result.artifact.storage_backend,'s3');checks++;
    assert(result.result.result.artifact.size>65536);checks++;
    const download=await fetch(base+'/runtime/api/jobs/'+item.receipt.id+'/artifact');
    assert.equal(download.status,200);checks++;
    assert.deepEqual(await download.json(),item.value);checks++;
    assert.equal((await fetch(base+'/artifacts/sha256/'+result.result.result.artifact.sha256)).status,404);checks++;
    evidence.artifact={...item,reference:result.result.result.artifact};
    await checkpoint();
  }
  await evaluate("studyFrame.remove();Reveal.slide(1);Reveal.next()");
  await check('Reveal.getIndices().h>=1','Static presentation navigation failed after worker loss');
  assert.equal(errors.length,0,JSON.stringify(errors));checks++;
  await writeFile(join(directory,'scientific-results.json'),JSON.stringify(evidence,null,2));
  console.log('PASS: '+checks+' real scientific consumer browser checks; physical isolation evidence belongs to the calling host gate');
  if(process.env.LCM_TEST_KEEP_CONSUMERS==='1'){
    // Retain the real presentation's UI socket while the separate terminal
    // browser uses that application run. This is a fixture handshake, not a
    // longer production disconnect grace or a synthetic runtime heartbeat.
    await writeFile(join(directory,'scientific-browser-ready'),'ready\n');
    const done=()=>access(join(directory,'terminal-browser-finished')).then(()=>true,()=>false);
    const deadline=Date.now()+400_000;
    while(!await done()){
      if(Date.now()>=deadline)throw Error('Physical terminal companion did not finish');
      await delay(100);
    }
  }
}catch(error){
  await screenshot('scientific-live-failure').catch(()=>{});
  const state=await evaluate("({text:window.view?.textContent,projection:window.view?.querySelector('[data-runtime-kind=execution]')?.dataset.resultProjection,job:window.job?{state:job.state,inputs:job.parameters,context:job.context(),control:job.client.state}:null})").catch(()=>null);
  await writeFile(join(directory,'scientific-live-failure.json'),JSON.stringify({error:String(error),partial:evidence,state,errors,traffic},null,2));
  throw error;
}finally{
  stop();
  if(browser && browser.exitCode===null){
    await Promise.race([new Promise(resolve=>browser.once('exit',resolve)),delay(5000)]);
    if(browser.exitCode===null){browser.kill('SIGKILL');await new Promise(resolve=>browser.once('exit',resolve));}
  }
  await writeFile(join(directory,'chrome.log'),chromeLog);
}
