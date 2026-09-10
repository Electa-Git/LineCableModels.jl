#!/usr/bin/env node
// Real Chrome + real shared browser client/renderer, with an explicit local HTTP
// state fixture. This is not a claim of broker execution or Bonito integration.
import assert from "node:assert/strict";
import {createServer} from "node:http";
import {spawn} from "node:child_process";
import {mkdtemp, readFile, writeFile} from "node:fs/promises";
import {tmpdir} from "node:os";
import {join} from "node:path";
import {randomUUID} from "node:crypto";
import {setTimeout as delay} from "node:timers/promises";

const scratch = await mkdtemp(join(tmpdir(), "lcm-runtime-controls."));
const root = new URL("../../", import.meta.url);
const assets = new Map(await Promise.all(["brand.css", "control-contract.css", "forms.css", "data-views.css", "runtime-controls.css",
  "runtime-client.js", "runtime-controls.js"].map(async name => [name, await readFile(new URL("assets/" + name, root))])));
assets.set("control.js", await readFile(new URL("runtime/ui/control.js", root)));
const themeInit = await readFile(new URL("assets/theme-init.html", root), "utf8");
const run = randomUUID(), lease = randomUUID(), epoch = randomUUID();
let down = false, uncertain = false, eventFailure = false, runState = "running", slowRefresh = false;
const actions = [], requests = [];
const control = {schema_version:1, enabled:true, preparation_control:true, broker:"online", administrator:true,
  profiles:[{id:"line-parameters", version:"1.0.0"}], provisioned:[{worker_id:"worker-a"}, {worker_id:"worker-b"}],
  workers:[{registration:{worker_id:"worker-a", profiles:["line-parameters"], capacity:2, state:"approved", revision:1},
    liveness:"online", occupied:0, preparation:"unknown", report:{capacity:2,
      profiles:[{profile_id:"line-parameters", version:"1.0.0"}]}}]};
let assignments = [];
let scientific = {channel:"online",phase:"idle",preparation:"cold",valid_for_ms:0,pending:true,accepted:true,
  progress:0,elapsed_seconds:0,output_lines:0,current_request_id:null,executor_id:null,executor_generation:0,
  preparation_key:null,failure:null,reason:"accepted"};
const eventBatch = {epoch, cursor:3, gap:false, records:[
  {sequence:1, at:"2026-09-07T12:00:00", code:"worker_reported", worker_id:"worker-a"},
  {sequence:2, at:"2026-09-07T12:00:01", code:"fixture_owned_run", run_id:run},
  {sequence:3, at:"2026-09-07T12:00:02", code:"fixture_other_run", run_id:randomUUID()}
]};
const config = {kind:"panel", run_id:run, roles:[{role:"parameters", profiles:["line-parameters"]}]};
const html = `<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
  ${themeInit}${["brand.css", "control-contract.css", "forms.css", "data-views.css", "runtime-controls.css"].map(name => `<link rel="stylesheet" href="/${name}">`).join("")}
  <style>body{margin:12px;background:var(--lc-bg);color:var(--lc-text)}main{max-width:1000px;margin:auto;min-width:0}</style></head>
  <body><main><label>Theme<select data-lcm-theme-selector class="lc-control-select lc-form-control"><option value="dark">Dark</option><option value="light">Light</option></select></label>
  <div id="control" data-lcm-runtime-controls='${JSON.stringify(config)}'></div></main>
  <script src="/runtime-client.js"></script><script src="/runtime-controls.js"></script><script src="/control.js"></script></body></html>`;
const server = createServer(async (req, res) => {
  const url = new URL(req.url, "http://fixture");
  if (url.pathname.startsWith("/runtime/api/")) {
    requests.push([req.method, url.pathname]);
    res.setHeader("Content-Type", "application/json");
    const reply = (value, status=200) => {res.statusCode=status; res.end(JSON.stringify(value));};
    if (down) return reply({error:"Runtime status is unavailable."}, 503);
    if (req.method === "GET") {
      if (slowRefresh) await delay(300);
      if (url.pathname === "/runtime/api/runs/" + run) return reply({id:run,application:"cable-study",state:runState,reason:runState === "stopped" ? "UI host stopped" : ""});
      if (url.pathname.endsWith("/control")) return reply(control);
      if (url.pathname.endsWith("/assignments")) return reply(assignments);
      if (url.pathname.endsWith("/science")) return reply(scientific);
      if (url.pathname.endsWith("/events")) return eventFailure ? reply({error:"Events unavailable"},503) :
        reply({...eventBatch, records:Number(url.searchParams.get("after")) >= eventBatch.cursor ? [] : eventBatch.records});
    }
    const chunks = []; for await (const chunk of req) chunks.push(chunk);
    const body = JSON.parse(Buffer.concat(chunks).toString());
    actions.push({method:req.method, path:url.pathname, body, headers:req.headers});
    if (req.method === "POST" && url.pathname.endsWith("/science")) {
      scientific = {...scientific, preparation_key:null,valid_for_ms:0};
      if (body.action === "prepare") scientific = {...scientific,phase:"preparing",preparation:"preparing",
        current_request_id:body.request_id,progress:0.25,elapsed_seconds:1.2,output_lines:2};
      else scientific = {...scientific,phase:"failed",preparation:"failed",failure:"canceled"};
      return reply(scientific,202);
    }
    if (req.method === "POST" && url.pathname.endsWith("/assignments")) {
      assignments = [{id:lease, run_id:run, role:body.role, profile:body.profile, worker_id:"worker-a", generation:1, state:"active", usable:true, preparation:"unknown"}];
      control.workers[0].occupied = 1;
      if (uncertain) { uncertain=false; return reply({error:"Control acknowledgement unavailable"},503); }
      return reply(assignments[0],202);
    }
    if (req.method === "DELETE") { assignments[0].state="releasing"; assignments[0].usable=false; return reply(assignments[0]); }
    if (req.method === "PATCH") {
      control.workers[0].registration.state=body.state; control.workers[0].registration.revision++; return reply(control.workers[0].registration);
    }
    return reply({error:"Fixture endpoint unavailable"},400);
  }
  const asset = assets.get(url.pathname.slice(1));
  res.setHeader("Content-Type", asset ? url.pathname.endsWith(".css") ? "text/css" : "text/javascript" : "text/html");
  res.end(asset ?? html);
});
await new Promise(resolve => server.listen(0, "127.0.0.1", resolve));
const base = "http://127.0.0.1:" + server.address().port;
let browser, socket;
const errors = [], timers = new Map();
try {
  browser = spawn(process.env.LCM_BROWSER || "google-chrome", ["--headless=new", "--disable-gpu", "--no-first-run", "--no-default-browser-check",
    "--remote-debugging-address=127.0.0.1", "--remote-debugging-port=0", "--user-data-dir=" + join(scratch,"chrome"), "about:blank"], {stdio:["ignore","ignore","pipe"]});
  let chromeLog = "";
  browser.stderr.on("data", chunk => {chromeLog = (chromeLog + chunk.toString()).slice(-65536);});
  browser.on("error", error => errors.push(error.message));
  let debugPort;
  for (let i=0; i<150; i++) {
    try { debugPort=(await readFile(join(scratch,"chrome","DevToolsActivePort"),"utf8")).split("\n")[0]; break; } catch {}
    if (browser.exitCode !== null) throw Error("Chrome exited before debugging was available: " + chromeLog);
    await delay(100);
  }
  assert.ok(debugPort, "Chrome startup deadline exceeded");
  const pages = await (await fetch("http://127.0.0.1:" + debugPort + "/json")).json();
  socket = new WebSocket(pages.find(page => page.type === "page").webSocketDebuggerUrl);
  await new Promise(resolve => socket.addEventListener("open", resolve, {once:true}));
  let sequence=0;
  socket.addEventListener("message", event => {
    const message=JSON.parse(event.data);
    if (message.method === "Runtime.exceptionThrown") errors.push(JSON.stringify(message.params.exceptionDetails));
    const pending=timers.get(message.id); if (!pending) return;
    timers.delete(message.id); clearTimeout(pending.timer);
    message.error ? pending.reject(Error(message.error.message)) : pending.resolve(message.result);
  });
  const command=(method,params={}) => new Promise((resolve,reject) => {
    const id=++sequence; const timer=setTimeout(() => {timers.delete(id);reject(Error("CDP timeout: " + method));},15000);
    timers.set(id,{resolve,reject,timer});socket.send(JSON.stringify({id,method,params}));
  });
  const evaluate=async expression => {
    const result=await command("Runtime.evaluate",{expression,awaitPromise:true,returnByValue:true});
    if (result.exceptionDetails) throw Error(result.exceptionDetails.exception?.description);
    return result.result.value;
  };
  const wait=async expression => {
    for (let i=0;i<100;i++) {if(await evaluate(expression))return; await delay(50);}
    throw Error("Timed out: " + expression);
  };
  await command("Runtime.enable"); await command("Page.enable");
  await command("Emulation.setDeviceMetricsOverride",{width:1280,height:1100,deviceScaleFactor:1,mobile:false});
  await command("Page.navigate",{url:base});
  await wait("document.querySelector('#control')?.dataset.runtimeStale === 'false'");
  await evaluate(`window.client=LineCableModelsRuntimeClient.acquire(${JSON.stringify(run)});
    window.button = text => [...document.querySelectorAll('#control button')].find(b => b.textContent === text);
    window.field = text => [...document.querySelectorAll('#control label')].find(l => l.firstChild.textContent === text).querySelector('select');
    window.choose=(name,value) => {const f=field(name); f.value=value;f.dispatchEvent(new Event('change'));};`);
  assert.equal(actions.length,0,"mount must allocate nothing");
  assert.equal(await evaluate("button('Assign worker').disabled"),true);
  await evaluate("choose('Profile','line-parameters');choose('Placement','pinned');choose('Worker','worker-a')");
  assert.equal(await evaluate("button('Assign worker').disabled"),false);
  slowRefresh = true;
  await evaluate("button('Refresh status').click()");
  assert.equal(await evaluate("button('Refreshing…').disabled && button('Refreshing…').getAttribute('aria-busy') === 'true'"),true);
  assert.equal(await evaluate("client.state.activity.some(e=>e.code==='refresh_started') && document.querySelector('[aria-label=\"Client action history\"]').textContent.includes('Refresh status requested')"),true);
  await wait("!client.state.refreshing"); slowRefresh = false;
  assert.equal(await evaluate("client.state.activity.at(-1).code === 'refresh_completed' && button('Refresh status').dataset.busy === 'false'"),true);
  assert.equal(await evaluate(`(() => {
    const log=document.querySelector('[aria-label="Structured control event history"]').textContent;
    return log.includes('worker_reported') && log.includes('fixture_owned_run') && !log.includes('fixture_other_run');
  })()`),true,"run diagnostics contain shared worker events and owned-run events, not another run's events");
  await evaluate(`window.statusChanges=0;window.statusObserver=new MutationObserver(r=>statusChanges+=r.length);
    statusObserver.observe(document.querySelector('.lc-runtime-connection [role="status"]'),{childList:true,subtree:true});`);
  await evaluate("client.refresh()");
  assert.equal(await evaluate("statusChanges"),0,'unchanged status must not repeat live announcements');
  await evaluate("statusObserver.disconnect()");
  runState = "stopped"; await evaluate("client.refresh()");
  assert.equal(await evaluate("button('Assign worker').disabled && document.querySelector('.lc-runtime-run-status').textContent.includes('start a new run')"),true);
  assert.equal(await evaluate(`document.querySelector('.lc-runtime-run-status a').getAttribute('href')`), '/runtime/runs/' + run);
  assert.equal(await evaluate("client.state.control.broker"),'online');
  runState = "running"; await evaluate("client.refresh()");
  await evaluate("field('Worker').focus(); client.refresh()");
  assert.equal(await evaluate("document.activeElement === field('Worker') && field('Worker').value === 'worker-a'"),true);
  control.workers[0].liveness="offline";
  await evaluate("client.refresh()");
  assert.equal(await evaluate("field('Worker').value === 'worker-a' && field('Worker').selectedOptions[0].textContent.includes('offline')"),true);
  assert.equal(await evaluate("button('Assign worker').disabled"),true,"offline pinned worker must remain visible but unavailable");
  // All selection controls keep the shared theme contract through repeated changes.
  const colors=[];
  for (const theme of ["light","dark","light"]) {
    await evaluate(`(() => {const themeSelect=document.querySelector('[data-lcm-theme-selector]');themeSelect.value=${JSON.stringify(theme)};themeSelect.dispatchEvent(new Event('change'));})()`);
    const styles=await evaluate(`(() => {const s=getComputedStyle(field('Profile')),o=getComputedStyle(field('Profile').options[1]);
      return {color:s.color,bg:s.backgroundColor,option:o.color,optionBg:o.backgroundColor,theme:document.documentElement.dataset.lcmResolvedTheme};})()`);
    assert.equal(styles.theme,theme);assert.notEqual(styles.color,styles.bg);assert.notEqual(styles.option,styles.optionBg);colors.push(styles);
    assert.equal(await evaluate(`(() => {
      const broker=document.querySelector('.lc-runtime-connection .lc-status-indicator');
      const offline=[...document.querySelectorAll('[aria-label="Worker inventory"] .lc-status-indicator')].find(n=>n.textContent==='offline');
      const sample=(n,key)=>{const swatch=document.createElement('span');swatch.style.color='var(--lc-'+key+')';n.append(swatch);const c=getComputedStyle(swatch).color;swatch.remove();return c;};
      return broker.dataset.tone==='success' && offline.dataset.tone==='danger' &&
        getComputedStyle(broker).color===sample(broker,'success') && getComputedStyle(offline).color===sample(offline,'danger') &&
        Number(getComputedStyle(broker).fontWeight)>=700 && Number(getComputedStyle(offline).fontWeight)>=700;
    })()`),true,theme+' semantic status styling');
    await writeFile(join(scratch,'status-' + theme + '.png'),Buffer.from((await command('Page.captureScreenshot',{format:'png'})).data,'base64'));
  }
  assert.deepEqual(colors[0],colors[2]);assert.notDeepEqual(colors[0],colors[1]);
  for (const width of [1280,390]) {
    await command("Emulation.setDeviceMetricsOverride",{width,height:1100,deviceScaleFactor:1,mobile:false});
    await evaluate("new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)))");
    assert.equal(await evaluate("document.documentElement.scrollWidth <= innerWidth"),true,"page must not overflow horizontally");
    const shot=await command("Page.captureScreenshot",{format:"png"});
    await writeFile(join(scratch,"controls-" + width + ".png"),Buffer.from(shot.data,"base64"));
  }
  await command("Emulation.setDeviceMetricsOverride",{width:1280,height:1100,deviceScaleFactor:1,mobile:false});
  control.workers[0].liveness="online";await evaluate("client.refresh()");
  uncertain=true;
  await evaluate("button('Assign worker').click()");
  await wait("!button('Retry same action').hidden && !client.state.pending");
  assert.equal(actions.length,1);await delay(2200);assert.equal(actions.length,1,"uncertain mutation must not auto-replay");
  assert.equal(await evaluate("button('Assign worker').disabled && button('Release assignment').disabled"),true);
  assert.equal(await evaluate("document.querySelector('[data-preparation]').dataset.preparation"),"cold");
  await evaluate("button('Retry same action').click()");
  await wait("button('Retry same action').hidden && !client.state.pending");
  assert.equal(actions.length,2);assert.deepEqual(actions[0].body,actions[1].body);
  assert.equal(actions[0].headers["x-lcm-request"],"1");
  assert.equal(actions[0].headers.origin,base);
  await wait("!button('Prepare executor').disabled");
  await evaluate("button('Prepare executor').click()");
  await wait("document.querySelector('[data-preparation]').dataset.preparation === 'preparing' && !client.state.pending");
  assert.equal(actions.at(-1).body.action,"prepare");
  assert.deepEqual(actions.at(-1).body.parameters,{});
  assert.equal(await evaluate("button('Prepare executor').disabled && !button('Cancel preparation').disabled"),true);
  assert.equal(await evaluate("document.querySelector('[data-preparation]').textContent.includes('25%')"),true);
  for (const theme of ['dark', 'light']) {
    await evaluate(`LineCableModelsTheme.select('${theme}')`);
    assert.equal(await evaluate(`(() => {const s=document.querySelector('[data-preparation] .lc-activity-status');
      const css=getComputedStyle(s,'::before'); return s.dataset.busy==='true' && css.animationName==='lc-activity-spin' &&
        css.borderTopColor!==css.borderBottomColor && getComputedStyle(s).color!==getComputedStyle(s.closest('section')).backgroundColor;})()`),true);
    const shot=await command('Page.captureScreenshot',{format:'png'});
    await writeFile(join(scratch,'preparing-'+theme+'.png'),Buffer.from(shot.data,'base64'));
  }
  await command('Emulation.setEmulatedMedia',{features:[{name:'prefers-reduced-motion',value:'reduce'}]});
  assert.equal(await evaluate("getComputedStyle(document.querySelector('[data-preparation] .lc-activity-status'),'::before').animationName"),'none');
  await command('Emulation.setEmulatedMedia',{features:[]});
  scientific.channel='offline';await evaluate('client.refreshScience()');
  assert.equal(await evaluate("document.querySelector('[data-preparation] .lc-activity-status').dataset.busy"),'false');
  scientific.channel='online';await evaluate('client.refreshScience()');
  await wait("document.querySelector('[data-preparation] .lc-activity-status').dataset.busy === 'true'");
  const target=actions.at(-1).body.request_id;
  await evaluate("button('Cancel preparation').click()");
  await wait("document.querySelector('[data-preparation]').dataset.preparation === 'failed' && !client.state.pending");
  assert.equal(await evaluate("document.querySelector('[data-preparation] .lc-activity-status').dataset.busy"),'false');
  assert.equal(actions.at(-1).body.target_id,target);
  assert.equal(actions.at(-1).body.action,"cancel");
  scientific={...scientific,phase:"idle",preparation:"ready",executor_id:lease,executor_generation:1,
    current_request_id:null,preparation_key:"b".repeat(64),valid_for_ms:5000,failure:null};
  await evaluate("client.refreshScience()");
  assert.equal(await evaluate("document.querySelector('[data-preparation]').dataset.preparation"),"ready");
  assert.equal(await evaluate("document.querySelector('[data-preparation] .lc-activity-status').dataset.busy"),'false');
  control.profiles[0].kind="terminal";
  const boundary=requests.length;
  await evaluate("client.refresh();");await delay(100);await evaluate("client.refreshScience()");
  assert.equal(await evaluate("document.querySelector('[data-preparation]').hidden && button('Prepare executor').disabled"),true,
    "Terminal roles must not offer scientific preparation");
  assert.equal(requests.slice(boundary).some(([,path])=>path.endsWith('/science')),false,
    "Terminal assignment was polled as a scientific executor");
  control.profiles[0].kind="scientific";await evaluate("client.refresh()");
  assert.equal(await evaluate("document.querySelector('[data-preparation]').hidden"),false);
  await evaluate("button('Release assignment').click()");
  await wait("client.state.assignments[0].state === 'releasing' && !client.state.pending");
  assert.equal(await evaluate("document.querySelector('[data-preparation]').dataset.preparation"),"unknown");
  assert.equal(await evaluate("button('Release assignment').disabled && button('Assign worker').disabled"),true);
  eventFailure=true;await evaluate("client.refreshEvents()");
  assert.equal(await evaluate("client.state.eventsStale && !client.state.stale"),true);
  down=true;await evaluate("client.refresh()");
  assert.equal(await evaluate("document.querySelector('#control').textContent.includes('Showing last-known') && field('Profile').disabled"),true);
  down=false;eventFailure=false;await evaluate("client.refresh()");
  await evaluate("choose('Registered worker','worker-a');choose('Registration state','draining');button('Apply registration').click()");
  await wait("client.state.control.workers[0].registration.state === 'draining' && !client.state.pending");
  assert.equal(actions.at(-1).body.expected_revision,1);
  control.administrator=false; await evaluate("client.refresh()");
  assert.equal(await evaluate("button('Apply registration').closest('section').hidden"),true);
  // A second component shares one run client; removing one preserves the other.
  await evaluate(`window.extra=document.createElement('div');document.querySelector('main').append(extra);
    LineCableModelsRuntimeControls.mount(extra,{kind:'diagnostics',run_id:${JSON.stringify(run)}})`);
  await evaluate("client.refresh()");
  assert.equal(await evaluate("client.listeners.size"),2);
  await evaluate("document.querySelector('#control').remove()");
  await wait("client.listeners.size === 1 && !client.closed");
  await evaluate("extra.remove()");await wait("client.closed");
  const count=requests.length;await delay(2200);assert.equal(requests.length,count,"detached components leaked polling");
  assert.deepEqual(errors,[]);
  console.log("Runtime control renderer: draft retention, offline pin, themes, compact layout, explicit retry, release, registration, diagnostics and teardown passed.");
} finally {
  socket?.close();for(const p of timers.values())clearTimeout(p.timer);
  if(browser && browser.exitCode===null) {
    browser.kill("SIGTERM");
    for(let i=0;i<100 && browser.exitCode===null;i++) await delay(50);
    if(browser.exitCode===null) browser.kill("SIGKILL");
  }
  server.closeAllConnections();await new Promise(resolve=>server.close(resolve));
  console.log("Runtime control browser diagnostics: " + scratch);
}
