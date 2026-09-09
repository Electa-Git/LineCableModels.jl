import assert from 'node:assert/strict';
import {installJobFixture} from './runtime_job_browser_fixture.mjs';
const [base, debug] = process.argv.slice(2);
const deadline = Date.now() + 150000;
let runs;
while (Date.now() < deadline) {
  runs = await (await fetch(base + '/runtime/api/runs')).json();
  if (runs.length === 2 && runs.every(run => run.state === 'running')) break;
  assert(!runs.some(run => run.state === 'failed'), JSON.stringify(runs));
  await new Promise(resolve => setTimeout(resolve, 200));
}
assert(runs?.length === 2 && runs.every(run => run.state === 'running'), 'Bonito children did not become ready');
const pages = await (await fetch(debug + '/json')).json();
const socket = new WebSocket(pages.find(page => page.type === 'page').webSocketDebuggerUrl);
await new Promise(resolve => socket.addEventListener('open', resolve, {once: true}));
let sequence = 0;
const pending = new Map(), errors = [], failedRequests = [];
socket.addEventListener('message', event => {
  const message = JSON.parse(event.data);
  if (message.method === 'Runtime.exceptionThrown') errors.push(message.params.exceptionDetails);
  if (message.method === 'Network.loadingFailed') failedRequests.push(message.params);
  if (message.method === 'Network.responseReceived' && message.params.response.status >= 400)
    failedRequests.push({url:message.params.response.url,status:message.params.response.status});
  const request = pending.get(message.id);
  if (!request) return;
  pending.delete(message.id);
  message.error ? request.reject(Error(message.error.message)) : request.resolve(message.result);
});
const command = (method, params = {}) => new Promise((resolve, reject) => {
  const id = ++sequence;
  pending.set(id, {resolve, reject}); socket.send(JSON.stringify({id, method, params}));
});
const evaluate = async expression => {
  const result = await command('Runtime.evaluate', {expression, returnByValue: true, awaitPromise: true});
  if (result.exceptionDetails) throw Error(result.exceptionDetails.exception?.description);
  return result.result.value;
};
const wait = async (expression, message) => {
  const until = Date.now() + 90000;
  while (Date.now() < until) {
    if (await evaluate(expression)) return;
    await new Promise(resolve => setTimeout(resolve, 100));
  }
  const documents = await evaluate(`(() => [document, ...[...document.querySelectorAll('iframe')].map(f=>f.contentDocument)].map(d=>({
    url:d.URL, ready:d.readyState, fixture:d.documentElement.dataset.fixtureReady,
    bonito:typeof d.defaultView.Bonito, text:d.body?.innerText.slice(0,1200),
    scripts:[...d.scripts].map(s=>({src:s.src,type:s.type,inline:s.src?'':s.textContent.slice(0,300)}))
  })))()`);
  throw Error(message + '\n' + JSON.stringify({errors, failedRequests, documents}));
};
const check = async (expression, message) => assert(await evaluate(expression), message);
const first = '/applications/runs/' + runs[0].id + '/';
const second = '/applications/runs/' + runs[1].id + '/';
try {
  await command('Page.enable'); await command('Runtime.enable'); await command('Network.enable');
  await command('Emulation.setDeviceMetricsOverride', {width:1600, height:1000, deviceScaleFactor:1, mobile:false});
  await command('Page.navigate', {url:base + first + 'counter'});
  await wait("document.readyState === 'complete' && document.documentElement.dataset.fixtureReady === 'true' && typeof Bonito !== 'undefined'", 'Counter did not mount through namespaced proxy');
  console.log('Mounted first Bonito counter');
  await evaluate(`window.count = doc => doc.querySelector('#count')?.textContent;
    window.addFrame = (id, src) => { const frame = document.createElement('iframe');
      frame.id=id; frame.src=src; frame.width='1400'; frame.height='900'; document.body.append(frame); };
    addFrame('peer', ${JSON.stringify(first + 'counter')});
    addFrame('other', ${JSON.stringify(second + 'counter')});
    addFrame('workbench', ${JSON.stringify(first + 'workbench')});
    addFrame('runtime-gallery', ${JSON.stringify(first + 'runtime-controls')});
    addFrame('runtime-workbench', ${JSON.stringify(first + 'runtime-workbench')});
    addFrame('terminal-gallery', ${JSON.stringify(first + 'terminal-view')});
    window.doc = id => document.getElementById(id).contentDocument;`);
  await wait("doc('peer').documentElement.dataset.fixtureReady === 'true' && doc('other').documentElement.dataset.fixtureReady === 'true'", 'Deck-style frames did not mount');
  await wait("!!document.getElementById('workbench').contentWindow.lcmXRay", 'Workbench X-ray did not mount through proxy');
  await wait("['runtime-gallery','runtime-workbench'].every(id => doc(id).querySelectorAll('[data-runtime-kind]').length >= 3 && [...doc(id).querySelectorAll('[data-runtime-kind]')].every(n => n.dataset.runtimeStale === 'false'))", 'Shared runtime controls did not mount in the gallery and workbench');
  await check("doc('runtime-workbench').body.textContent.includes('Worker control is not configured') && !doc('runtime-workbench').querySelector('.lc-runtime-controls button:not(:disabled)').textContent.includes('Assign')", 'Disabled runtime configuration was reported as a connection failure');
  console.log('Mounted shared runtime controls in gallery and workbench frames');
  await wait("['terminal-gallery','runtime-workbench'].every(id => doc(id).querySelector('.lc-runtime-terminal .xterm'))", 'Private terminal did not mount through Bonito assets and onload');
  await check("['terminal-gallery','runtime-workbench'].every(id => [...doc(id).querySelectorAll('.lc-runtime-terminal button')].find(b=>b.textContent==='Connect')?.disabled)", 'Unassigned terminal allowed connection');
  console.log('Mounted the same private terminal in real Bonito gallery/deck-style and workbench hosts');
  await evaluate("doc('peer').querySelector('#increment').click()");
  await wait("count(document)==='1' && count(doc('peer'))==='1'", 'Bonito binary updates did not reach both same-run frames');
  await check("count(doc('other'))==='0'", 'Observable state leaked across application processes');
  await evaluate(`window.savedPeer=doc('peer'); window.savedOther=doc('other'); window.savedWorkbench=doc('workbench');
    window.wb=document.getElementById('workbench').contentWindow;
    window.rootIdentity=wb.lcmXRay.root;`);
  await evaluate(`window.peerWindow=document.getElementById('peer').contentWindow;
    window.peerTransport=peerWindow.WEBSOCKET;
    window.reconnected=false;
    peerTransport.on_open(()=>{window.reconnected=true;}); peerTransport.close();`);
  await wait('reconnected && peerTransport.isopen()', 'Surviving Bonito session did not reconnect');
  await evaluate("doc('peer').querySelector('#increment').click()");
  await wait("count(document)==='2' && count(doc('peer'))==='2'", 'Reconnected session lost its observable binding');
  await check("savedPeer===doc('peer') && count(doc('other'))==='0'", 'Reconnect remounted or crossed the run boundary');
  for (const theme of ['light', 'dark', 'light', 'dark']) {
    await evaluate(`(() => { const input = doc('workbench').querySelector('[data-lc-wb-theme-selector]');
      input.value=${JSON.stringify(theme)}; input.dispatchEvent(new wb.Event('change')); })()`);
    await wait(`[document,doc('peer'),doc('other'),doc('workbench'),doc('runtime-gallery'),doc('runtime-workbench'),doc('terminal-gallery')].every(d => d.documentElement.dataset.lcmResolvedTheme===${JSON.stringify(theme)})`, 'Theme did not propagate across mounted hosts');
    await check(`(() => {const colors=['terminal-gallery','runtime-workbench'].map(id=>{const d=doc(id),s=d.defaultView.getComputedStyle(d.querySelector('.lc-terminal-viewport'));return [s.color,s.backgroundColor,s.borderColor];});return JSON.stringify(colors[0])===JSON.stringify(colors[1]) && colors[0][0]!==colors[0][1];})()`, 'Private terminal styles differ between gallery and workbench');
    await check(`(() => {const styles=['runtime-gallery','runtime-workbench'].map(id => {const d=doc(id),s=d.defaultView.getComputedStyle(d.querySelector('.lc-runtime-fields select'));return [s.color,s.backgroundColor];});return JSON.stringify(styles[0])===JSON.stringify(styles[1]) && styles.every(s=>s[0]!==s[1]);})()`, 'Worker selector styling differs between gallery and workbench');
    await check("savedPeer===doc('peer') && savedOther===doc('other') && savedWorkbench===doc('workbench') && rootIdentity===wb.lcmXRay.root && count(document)==='2'", 'Theme switching remounted a session');
  }
  await evaluate(`wb.lcmXRay.enable(); doc('workbench').querySelector('.lc-wb-menubar').dispatchEvent(new wb.MouseEvent('click',{bubbles:true,composed:true}));
    window.xr=doc('workbench').querySelector('.lc-xray-host').shadowRoot;`);
  await wait("!xr.querySelector('.xray-panel').hidden && xr.querySelector('.xray-body').textContent.includes('WorkbenchUI.MenuBar')", 'Owned X-ray metadata did not survive proxy');
  await check("!!xr.querySelector('[data-css-property]') && !xr.querySelector('.xray-properties input')", 'X-ray CSS preview controls did not mount');
  await evaluate(`window.rw=document.getElementById('runtime-workbench').contentWindow;rw.lcmXRay.enable();
    doc('runtime-workbench').querySelector('.lc-runtime-controls').dispatchEvent(new rw.MouseEvent('click',{bubbles:true,composed:true}));`);
  await wait("doc('runtime-workbench').querySelector('.lc-xray-host').shadowRoot.querySelector('.xray-body').textContent.includes('WorkerSelector')", 'Worker selector X-ray metadata did not mount');
  await evaluate("doc('runtime-workbench').querySelector('.lc-runtime-terminal').dispatchEvent(new rw.MouseEvent('click',{bubbles:true,composed:true}))");
  await wait("doc('runtime-workbench').querySelector('.lc-xray-host').shadowRoot.querySelector('.xray-body').textContent.includes('JuliaTerminal')", 'Terminal X-ray metadata did not mount');
  await check("!doc('runtime-workbench').querySelector('.lc-xray-host').shadowRoot.querySelector('.xray-body').textContent.includes('writer_id')", 'Private terminal identity leaked into X-ray');
  await evaluate(`addFrame('scientific', ${JSON.stringify(first+'job-view')})`);
  await wait("!!doc('scientific').querySelector('[data-runtime-kind=execution]') && typeof document.getElementById('scientific').contentWindow.Bonito !== 'undefined'", 'Scientific component did not mount through real Bonito');
  await evaluate(`(async () => { window.jw=document.getElementById('scientific').contentWindow;
    window.displayIdentity=doc('scientific').querySelector('#scientific-display');
    await jw.eval(${JSON.stringify('('+installJobFixture.toString()+')()')});
    window.jf=jw.__scientificFixture;
    window.jobButton=label=>[...doc('scientific').querySelectorAll('button')].find(b=>b.textContent===label); })()`);
  await wait("jf.job.state.canRun", 'Fixture prepared state did not reach the shared job control');
  await check("jf.submissions.length===0 && displayIdentity.textContent==='No value'", 'Mount or readiness submitted work');
  await evaluate("jobButton('Run calculation').click()");
  await wait("jf.job.state.receipt?.state==='queued'", 'Run button did not submit the explicit intent');
  await evaluate("jf.complete()");
  await wait("displayIdentity.textContent==='3 · current'", 'Successful browser result did not reach the Julia Observable');
  await evaluate("doc('scientific').querySelector('#input-three').click()");
  await check("displayIdentity.textContent==='3 · current'", 'Identical input echo marked the result outdated');
  await evaluate("doc('scientific').querySelector('#input-four').click()");
  await wait("displayIdentity.textContent==='3 · outdated' && jf.job.parameters.value===4", 'Julia input change did not invalidate the displayed value');
  await check("jf.submissions.length===1", 'Changing an input automatically submitted work');
  await evaluate("jobButton('Run calculation').click()");
  await wait("jf.submissions.length===2 && jf.job.state.receipt?.state==='queued'", 'Second explicit job did not submit');
  await evaluate("doc('scientific').querySelector('#input-five').click()");
  await wait("jf.job.parameters.value===5", 'New Julia input did not reach the browser draft');
  await evaluate("jf.complete()");
  await check("displayIdentity.textContent==='3 · outdated' && jf.job.state.superseded", 'Late completion replaced the Julia view');
  await evaluate("jf.loseSubmit=true;jobButton('Run calculation').click()");
  await wait("jf.job.state.phase==='unconfirmed'", 'Lost acknowledgement was not visible');
  await check("!jobButton('Retry acknowledgement').hidden && jobButton('Run calculation').disabled", 'Unconfirmed submission allows duplicate Run');
  await evaluate("jobButton('Retry acknowledgement').click()");
  await wait("jf.submissions.length===4 && jf.job.state.receipt?.state==='queued'", 'Explicit retry did not recover the same receipt');
  await check("JSON.stringify(jf.submissions[2])===JSON.stringify(jf.submissions[3])", 'Browser retry changed request identity or inputs');
  await evaluate("jf.complete()");
  await wait("displayIdentity.textContent==='5 · current'", 'Retried successful result did not reach the Julia view');
  for (const theme of ['light','dark']) {
    await evaluate(`(() => {const input=doc('workbench').querySelector('[data-lc-wb-theme-selector]');
      input.value=${JSON.stringify(theme)};input.dispatchEvent(new wb.Event('change'));})()`);
    await wait(`doc('scientific').documentElement.dataset.lcmResolvedTheme===${JSON.stringify(theme)}`, 'Scientific control theme did not propagate');
    await check(`(() => { const controls=['runtime-gallery','scientific'].map(id=>{
      const d=doc(id),button=[...d.querySelectorAll('button')].find(b=>b.textContent==='Run calculation');
      const disabled=button.disabled;button.disabled=true;
      const s=d.defaultView.getComputedStyle(button),result=[s.color,s.backgroundColor,s.borderColor];
      button.disabled=disabled;return result;
    });return JSON.stringify(controls[0])===JSON.stringify(controls[1]) && controls[1][0]!==controls[1][1];})()`, 'Scientific control styling drifted from the gallery');
    await check("displayIdentity===doc('scientific').querySelector('#scientific-display') && displayIdentity.textContent==='5 · current'", 'Theme switch remounted the scientific display');
  }
  await check(`(() => {
    window.beforeDraftSubmissions=jf.submissions.length;
    const input=doc('scientific').querySelector('#draft-field');
    input.value='4';input.dispatchEvent(new jw.Event('input',{bubbles:true}));
    jobButton('Run calculation').click();
    input.value='5';input.dispatchEvent(new jw.Event('input',{bubbles:true}));
    const root=doc('scientific').querySelector('[data-runtime-kind=execution]');
    root.__lcmJobInputAck({epoch:1,inputs:{draft_id:'stale-fixture-ack',parameters:{value:4}}});
    return jf.job.state.inputsPending && !jf.job.state.canRun && jobButton('Run calculation').disabled &&
      jf.submissions.length===beforeDraftSubmissions;
  })()`, 'A rapid edit or stale acknowledgement let Run submit the previous field value');
  await wait("!jf.job.state.inputsPending && jf.job.parameters.value===5 && jf.job.state.canRun", 'Latest field acknowledgement did not enable the canonical draft');
  await check(`(() => {doc('scientific').querySelector('[data-runtime-kind=execution]').__lcmJobInputAck(
    {epoch:2,inputs:{draft_id:'duplicate-fixture-ack',parameters:{value:4}}});
    return jf.job.parameters.value===5 && !jf.job.state.inputsPending;})()`, 'Duplicate input acknowledgement replaced the canonical draft');
  await check("jf.submissions.length===beforeDraftSubmissions", 'Input acknowledgement implicitly submitted a calculation');
  await evaluate("jobButton('Run calculation').click()");
  await wait("jf.submissions.length===beforeDraftSubmissions+1 && jf.job.state.receipt?.state==='queued'", 'Explicit Run did not capture the acknowledged draft');
  await check("jf.submissions.at(-1).parameters.value===5", 'Run captured an older field value');
  await evaluate("jf.complete()");
  await wait("displayIdentity.textContent==='5 · current'", 'Acknowledged field result did not reach the shared view');
  await check(`(() => {const input=doc('scientific').querySelector('#draft-checkbox');input.click();
    return jf.job.state.inputsPending && !jf.job.state.canRun;})()`, 'Checkbox change escaped the input fence');
  await wait("!jf.job.state.inputsPending && jf.job.parameters.value===4", 'Checkbox acknowledgement preceded its actual value');
  await check(`(() => {const input=doc('scientific').querySelector('#draft-choice');input.value='5';
    input.dispatchEvent(new jw.Event('change',{bubbles:true}));
    return jf.job.state.inputsPending && !jf.job.state.canRun;})()`, 'Dropdown change escaped the input fence');
  await wait("!jf.job.state.inputsPending && jf.job.parameters.value===5", 'Dropdown acknowledgement preceded its actual value');
  await check("jf.submissions.length===beforeDraftSubmissions+1", 'Choice acknowledgements submitted implicit jobs');
  await evaluate("doc('scientific').querySelector('#input-invalid').click()");
  await wait("displayIdentity.textContent==='5 · outdated' && !jf.job.inputsValid && jobButton('Run calculation').disabled", 'Invalid Julia draft did not fail safely');
  await check(`(() => {const root=doc('scientific').createElement('div');root.textContent='unchanged';doc('scientific').body.append(root);
    const api=jw.LineCableModelsRuntimeClient, run=crypto.randomUUID(), client=api.acquire(run);
    let rejected=false;try {jw.LineCableModelsRuntimeControls.mount(root,{kind:'execution',run_id:run,role:'parameters',operation:'Core.eval'});} catch {rejected=true;}
    const clean=rejected && client.closed && root.textContent==='unchanged';root.remove();return clean;
  })()`, 'Invalid execution configuration leaked an acquired client or altered the document');
  await evaluate("jf.restore();document.getElementById('scientific').remove()");
  await wait("jf.job.closed && jf.job.client.closed", 'Scientific frame teardown leaked its controller/client');
  console.log('PASS: real Bonito scientific input/result bindings with explicit mocked job HTTP; last-good view, late result fencing, explicit retry, theme parity and teardown');
  assert.equal(errors.length, 0, JSON.stringify(errors));
  const stopped = await fetch(base + '/runtime/api/runs/' + runs[0].id, {
    method:'DELETE', headers:{Origin:base, 'X-LCM-Request':'1'}});
  assert.equal(stopped.status, 200);
  assert.equal((await fetch(base + first + 'counter')).status, 503);
  assert.equal((await fetch(base + '/health')).status, 200);
  await evaluate("doc('other').querySelector('#increment').click()");
  await wait("count(doc('other'))==='1'", 'Stopping one UI host disabled the surviving run');
  await command('Page.navigate', {url:base + '/runtime/runs/' + runs[0].id});
  await wait("document.querySelector('#runtime-status')?.textContent==='stopped' && !document.querySelector('#runtime-restart').hidden", 'Stopped run did not expose a clean restart');
  for (const theme of ['light','dark']) {
    await evaluate(`LineCableModelsTheme.select(${JSON.stringify(theme)})`);
    await check(`document.documentElement.dataset.lcmResolvedTheme===${JSON.stringify(theme)} &&
      getComputedStyle(document.body).backgroundColor !== getComputedStyle(document.querySelector('h1')).color &&
      getComputedStyle(document.documentElement).caretColor==='rgba(0, 0, 0, 0)'`, 'Runtime recovery surface lost the shared theme/caret contract');
  }
  await evaluate("document.querySelector('#runtime-restart').click()");
  await wait(`location.pathname.startsWith('/applications/runs/') && !location.pathname.startsWith(${JSON.stringify(first)}) &&
    document.documentElement.dataset.fixtureReady==='true'`, 'Explicit clean restart did not open a new owned run');
  await check("document.querySelector('#count').textContent==='0'", 'Clean restart falsely restored lost state');
  console.log('PASS: real Bonito assets/binary sockets, shared run frames, independent state, reconnect, theme identity, X-ray, owned stop, surviving host and clean restart');
} finally {
  socket.close();
}
