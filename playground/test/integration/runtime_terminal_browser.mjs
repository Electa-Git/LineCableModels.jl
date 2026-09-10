#!/usr/bin/env node
// Real Chrome/xterm with explicit HTTP and terminal fixtures, not a runtime/PTY claim.
import assert from "node:assert/strict";
import {createServer} from "node:http";
import {spawn} from "node:child_process";
import {mkdtemp,readFile,writeFile} from "node:fs/promises";
import {tmpdir} from "node:os";
import {join} from "node:path";
import {randomUUID} from "node:crypto";
import {setTimeout as delay} from "node:timers/promises";

const scratch=await mkdtemp(join(tmpdir(),"lcm-terminal-browser."));
const live=process.env.LCM_TERMINAL_TEST_ORIGIN;
const startupMs=live ? Number(process.env.LCM_TERMINAL_TEST_STARTUP_MS) : 7500;
assert.ok(Number.isFinite(startupMs) && startupMs>0 && startupMs<=120000,
  'Live terminal fixture requires its finite configured startup bound');
if(live) {const url=new URL(live);assert.equal(url.protocol,"http:");assert.equal(url.hostname,"127.0.0.1");assert.equal(url.pathname,"/");}
const root=new URL("../../",import.meta.url);
const names=["brand.css","control-contract.css","forms.css","runtime-controls.css","runtime-terminal.css",
  "runtime-client.js","runtime-terminal-client.js","runtime-terminal.js","vendor/runtime-terminal.bundle.js","vendor/runtime-terminal.bundle.css"];
const assets=new Map(await Promise.all(names.map(async name=>["/"+name,await readFile(new URL("assets/"+name,root))])));
assets.set("/fixture.js",await readFile(new URL("runtime_terminal_browser_fixture.js",import.meta.url)));
const themeInit=await readFile(new URL("assets/theme-init.html",root),"utf8");
const run=randomUUID(),lease={id:randomUUID(),run_id:run,role:"terminal",profile:"julia-terminal",generation:1,
  worker_id:"fixture",worker_boot:randomUUID(),usable:true,state:"active"};
const inventory={schema_version:1,enabled:true,preparation_control:false,broker:"online",administrator:false,
  workers:[],provisioned:[],profiles:[{id:"julia-terminal",kind:"terminal",isolation:"container"}]};
const config=JSON.stringify({run_id:run,role:"terminal",title:"Julia REPL",rows:18});
const html=`<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
  <title>Terminal browser fixture</title>${themeInit}${names.filter(n=>n.endsWith("css")).map(n=>`<link rel="stylesheet" href="/${n}">`).join("")}
  <style>body{margin:12px;background:var(--lc-bg);color:var(--lc-text)}main{max-width:1000px;margin:auto;min-width:0}</style>
  </head><body><main><div id="terminal" data-lcm-runtime-terminal='${config}'></div></main>
  ${names.filter(n=>n.endsWith("js")).map(n=>`<script src="/${n}"></script>`).join("")}<script src="/fixture.js"></script></body></html>`;
const server=createServer((req,res)=>{
  if(req.url===`/runtime/api/runs/${run}`) {
    res.setHeader("Content-Type","application/json");res.end(JSON.stringify({id:run,application:"cable-study",state:"running",reason:""}));return;
  }
  if(req.url==="/runtime/api/control" || req.url===`/runtime/api/runs/${run}/assignments`) {
    res.setHeader("Content-Type","application/json");res.end(JSON.stringify(req.url.endsWith("control")?inventory:[lease]));return;
  }
  const asset=assets.get(req.url);
  res.setHeader("Content-Type",asset?req.url.endsWith("css")?"text/css":"text/javascript":"text/html");res.end(asset??html);
});
await new Promise(resolve=>server.listen(0,"127.0.0.1",resolve));
let browser,socket;const pending=new Map(),errors=[];let checks=0;
process.once("SIGTERM",()=>{
  process.exitCode=1;
  for(const item of pending.values()){clearTimeout(item.timer);item.reject(Error("Terminal browser fixture interrupted"));}
  pending.clear();socket?.close();if(browser && browser.exitCode===null)browser.kill("SIGTERM");
});
try {
  browser=spawn(process.env.LCM_BROWSER||"google-chrome",["--headless=new","--disable-gpu","--no-first-run","--no-default-browser-check",
    "--remote-debugging-address=127.0.0.1","--remote-debugging-port=0","--user-data-dir="+join(scratch,"chrome"),"about:blank"],{stdio:["ignore","ignore","pipe"]});
  let chromeLog="";browser.stderr.on("data",chunk=>{chromeLog=(chromeLog+chunk.toString()).slice(-65536);});
  browser.on("error",error=>errors.push(error.message));
  let port;
  for(let i=0;i<150;i++) {
    try {port=(await readFile(join(scratch,"chrome","DevToolsActivePort"),"utf8")).split("\n")[0];break;} catch{}
    if(browser.exitCode!==null) throw Error("Fixture Chrome exited: "+chromeLog);
    await delay(100);
  }
  assert.ok(port);
  const pages=await(await fetch("http://127.0.0.1:"+port+"/json")).json();
  socket=new WebSocket(pages.find(page=>page.type==="page").webSocketDebuggerUrl);
  await new Promise(resolve=>socket.addEventListener("open",resolve,{once:true}));
  let sequence=0;
  socket.addEventListener("message",event=>{
    const message=JSON.parse(event.data);
    if(message.method==="Runtime.exceptionThrown") errors.push(JSON.stringify(message.params.exceptionDetails));
    const item=pending.get(message.id);if(!item)return;
    pending.delete(message.id);clearTimeout(item.timer);message.error?item.reject(Error(message.error.message)):item.resolve(message.result);
  });
  const command=(method,params={})=>new Promise((resolve,reject)=>{
    const id=++sequence,timer=setTimeout(()=>{pending.delete(id);reject(Error("CDP deadline: "+method));},15000);
    pending.set(id,{resolve,reject,timer});socket.send(JSON.stringify({id,method,params}));
  });
  const evaluate=async expression=>{
    const result=await command("Runtime.evaluate",{expression,awaitPromise:true,returnByValue:true});
    if(result.exceptionDetails)throw Error(result.exceptionDetails.exception?.description);
    return result.result.value;
  };
  const wait=async (expression,milliseconds=7500)=>{
    const deadline=performance.now()+milliseconds;
    while(performance.now()<deadline){if(await evaluate(expression))return;await delay(50);}
    const shot=await command("Page.captureScreenshot",{format:"png"});
    await writeFile(join(scratch,"failure.png"),Buffer.from(shot.data,"base64"));
    const context=await evaluate(`(() => {
      const root=document.querySelector('#terminal');
      const config=root ? JSON.parse(root.dataset.lcmRuntimeTerminal) : null;
      const client=config && LineCableModelsRuntimeClient.acquire(config.run_id);
      return {phase:root?.dataset.terminalPhase,status:document.querySelector('.lc-terminal-status')?.textContent,
        runtime:client?.state,output:typeof text==='function'?text().slice(-4000):null};
    })()`);
    context.browserErrors=errors;
    await writeFile(join(scratch,"failure.json"),JSON.stringify(context,null,2));
    console.error("Explicit terminal fixture diagnostics: "+scratch);
    throw Error("Browser state deadline: "+expression);
  };
  await command("Runtime.enable");await command("Page.enable");
  await command("Emulation.setDeviceMetricsOverride",{width:1280,height:800,deviceScaleFactor:1,mobile:false});
  await command("Page.navigate",{url:live || "http://127.0.0.1:"+server.address().port});
  if(live) {
    await wait("globalThis.__liveTerminal && [...document.querySelectorAll('button')].some(b=>b.textContent==='Connect' && !b.disabled)");
    await evaluate("window.term=__liveTerminal;window.button=label=>[...document.querySelectorAll('#terminal button')].find(b=>b.textContent===label);window.text=()=>Array.from({length:term.buffer.active.length},(_,i)=>term.buffer.active.getLine(i)?.translateToString(true)).join('\\n')");
    const ready=(milliseconds=startupMs)=>wait(`(() => {
      const phase=document.querySelector('#terminal').dataset.terminalPhase;
      if(['failed','exited','uncertain','disconnected'].includes(phase))
        throw Error('Terminal failed before readiness: '+phase+' / '+document.querySelector('.lc-terminal-status').textContent);
      return phase==='ready' && !term.options.disableStdin;
    })()`,milliseconds);
    const enter=async code=>{
      await evaluate("term.focus()");await command("Input.insertText",{text:code});
      await command("Input.dispatchKeyEvent",{type:"keyDown",key:"Enter",code:"Enter",text:"\r",windowsVirtualKeyCode:13});
      await command("Input.dispatchKeyEvent",{type:"keyUp",key:"Enter",code:"Enter",windowsVirtualKeyCode:13});
    };
    const started=performance.now();
    await evaluate("button('Connect').click()");await ready();
    console.log('Live terminal guarded startup: '+((performance.now()-started)/1000).toFixed(3)+' s');
    await enter('terminal_browser_value=41; println("LCM_BROWSER:", terminal_browser_value+1)');
    await wait("text().includes('LCM_BROWSER:42')");checks++;
    await enter('browser_λ="λ"; println("UNICODE_BROWSER:", browser_λ)');
    await wait("text().includes('UNICODE_BROWSER:λ')");checks++;
    await evaluate(`(() => {const data=new DataTransfer();data.setData('text/plain',${JSON.stringify('begin\nbrowser_sum=sum(1:3)\nprintln("MULTI_BROWSER:", browser_sum)\nend')});term.textarea.dispatchEvent(new ClipboardEvent('paste',{clipboardData:data,bubbles:true,cancelable:true}));})()`);
    await enter('');
    await wait("text().includes('MULTI_BROWSER:6')");checks++;
    await enter('println("HISTORY_BROWSER:",42)');await wait("text().includes('HISTORY_BROWSER:42')");
    await command("Input.dispatchKeyEvent",{type:"keyDown",key:"ArrowUp",code:"ArrowUp",windowsVirtualKeyCode:38});
    await command("Input.dispatchKeyEvent",{type:"keyUp",key:"ArrowUp",code:"ArrowUp",windowsVirtualKeyCode:38});
    await enter('');await wait("text().split('HISTORY_BROWSER:42').length>=3");checks++;
    await command("Input.insertText",{text:"printl"});
    await command("Input.dispatchKeyEvent",{type:"keyDown",key:"Tab",code:"Tab",text:"\t",windowsVirtualKeyCode:9});
    await command("Input.dispatchKeyEvent",{type:"keyUp",key:"Tab",code:"Tab",windowsVirtualKeyCode:9});
    await delay(200);await enter('("COMPLETE_BROWSER:",43)');await wait("text().includes('COMPLETE_BROWSER:43')");checks++;
    await enter('sleep(30)');await delay(300);await evaluate("button('Interrupt').click()");
    await wait("text().includes('InterruptException')");
    await enter('println("AFTER_INTERRUPT:",1)');await wait("text().includes('AFTER_INTERRUPT:1')");checks++;
    for(const selected of ["light","dark"]) {
      await evaluate(`LineCableModelsTheme.select(${JSON.stringify(selected)})`);await delay(150);
      assert.equal(await evaluate("term.options.theme.background === getComputedStyle(document.querySelector('#terminal')).getPropertyValue('--lc-console-bg').trim()"),true);checks++;
    }
    await command("Emulation.setDeviceMetricsOverride",{width:900,height:760,deviceScaleFactor:1,mobile:false});await delay(500);
    await enter('println("PTY_SIZE:", displaysize(stdout))');
    await wait("text().includes('PTY_SIZE:(' + term.rows + ', ' + term.cols + ')')");checks++;
    assert.equal(await evaluate("document.querySelector('.xterm-screen').getBoundingClientRect().bottom <= document.querySelector('.lc-terminal-screen').getBoundingClientRect().bottom + 0.5"),true);checks++;
    await evaluate("button('Disconnect').click()");await wait("!button('Reconnect').disabled");await delay(700);
    await evaluate("button('Reconnect').click()");await ready(7500);
    await enter('println("RETAINED_BROWSER:", terminal_browser_value)');await wait("text().includes('RETAINED_BROWSER:41')");checks++;
    await evaluate("button('Restart').click();button('Restart Julia').click()");
    await wait("document.querySelector('#terminal').dataset.terminalPhase!=='ready'");await ready();
    await enter('println("NEW_BROWSER:", isdefined(Main, :terminal_browser_value))');await wait("text().includes('NEW_BROWSER:false')");checks++;
    await enter('for i in 1:2000; println("row ",i); end; println("OUTPUT_FINISHED")');
    await wait("text().includes('row 2000') && text().includes('OUTPUT_FINISHED')");checks++;
    assert.equal((await fetch(live+"health")).status,200);checks++;
    const shot=await command("Page.captureScreenshot",{format:"png"});await writeFile(join(scratch,"terminal-live.png"),Buffer.from(shot.data,"base64"));
    await evaluate("document.querySelector('#terminal').remove()");
    await wait("LineCableModelsRuntimeClient.acquire("+JSON.stringify(JSON.parse(await evaluate("document.querySelector('#terminal-config').textContent")).run_id)+").listeners.size===0");
    assert.deepEqual(errors,[]);checks++;
    console.log(`Live terminal browser: ${checks} assertions passed across Chrome, gateway, TLS broker and real fixture Julia REPL. Screenshot: ${scratch}`);
  } else {
  await wait("globalThis.__terminalFixture?.transport.state.canConnect");
  await evaluate("window.f=__terminalFixture; window.term=f.renderers[0]; window.text=()=>Array.from({length:term.buffer.active.length},(_,i)=>term.buffer.active.getLine(i)?.translateToString(true)).join('\\n')");
  assert.equal(await evaluate("f.sent.length"),0);checks++;
  await evaluate("f.button('Connect').click()");await wait("f.transport.state.phase === 'starting'");
  assert.equal(await evaluate("term.options.disableStdin"),true);checks++;
  await evaluate("f.phase='ready';f.append('julia> ')");await wait("f.transport.state.canInput && text().includes('julia>')");
  assert.equal(await evaluate("term.options.disableStdin"),false);checks++;
  // Ordinary input/output polls are not a loading state. Sample every transport
  // publication (including its busy moment), not only settled screenshots.
  for (const selected of ["dark", "light"]) {
    await evaluate(`LineCableModelsTheme.select(${JSON.stringify(selected)})`);
    for (const width of [1280,520]) {
      await command("Emulation.setDeviceMetricsOverride",{width,height:800,deviceScaleFactor:1,mobile:false});
      await delay(350);
      await evaluate(`window.layoutSamples=[]; window.captureLayout=()=>{
        const box=n=>{const r=n.getBoundingClientRect();return [r.x,r.y,r.width,r.height];};
        return {busy:document.querySelector('.lc-terminal-phase').dataset.busy,
          boxes:[...document.querySelectorAll('.lc-terminal-heading, .lc-runtime-actions, .lc-terminal-viewport, .lc-runtime-actions button')].map(box)};
      }; window.stableLayout=captureLayout(); window.unwatchLayout=f.transport.subscribe(()=>layoutSamples.push(captureLayout())); term.focus();`);
      await command("Input.insertText",{text:"steady"});
      await wait("text().includes('steady')");await delay(400);
      assert.equal(await evaluate("layoutSamples.length>2 && layoutSamples.every(s=>s.busy==='false' && JSON.stringify(s.boxes)===JSON.stringify(stableLayout.boxes))"),true,
        `${selected}/${width}: typing or output polling shifted terminal controls or showed loading`);checks++;
      await evaluate("unwatchLayout()");
    }
  }
  await command("Emulation.setDeviceMetricsOverride",{width:1280,height:800,deviceScaleFactor:1,mobile:false});await delay(350);
  await evaluate("window.opensBeforeLoss=f.sent.filter(p=>p.action==='open').length;f.transport.client.notify({stale:true})");
  assert.equal(await evaluate("document.querySelector('.lc-terminal-status').textContent.includes('Assignment or live inventory changed')"),true);checks++;
  assert.equal(await evaluate("term.options.disableStdin && !term.options.cursorBlink"),true);checks++;
  await evaluate("f.transport.client.notify({stale:false})");await delay(150);
  assert.equal(await evaluate("f.sent.filter(p=>p.action==='open').length===opensBeforeLoss && !f.transport.state.connected"),true);checks++;
  await evaluate("f.button('Reconnect').click()");await wait("f.transport.state.canInput");
  await delay(200);
  // Nonzero round-trip time exposes busy states that an immediate mock reply
  // hides. Read/input acknowledgements must not animate or move the REPL shell.
  for (const selected of ['dark','light']) {
    await evaluate(`LineCableModelsTheme.select(${JSON.stringify(selected)});f.replyDelay=65`);
    await delay(250);
    await evaluate(`window.terminalSamples=[];window.watchTerminal=true;
      window.sampleTerminal=()=>{if(!watchTerminal)return;const r=document.querySelector('.lc-terminal-viewport').getBoundingClientRect();
        terminalSamples.push({top:r.top,height:r.height,width:r.width,busy:document.querySelector('.lc-terminal-phase').dataset.busy,
          cols:term.cols,rows:term.rows});requestAnimationFrame(sampleTerminal);};requestAnimationFrame(sampleTerminal);term.focus()`);
    for (const letter of 'stable_repl') { await command('Input.insertText',{text:letter}); await delay(25); }
    await delay(350);
    const samples=await evaluate('watchTerminal=false;terminalSamples');
    assert.ok(samples.length>10);
    assert.ok(samples.every(s=>s.busy==='false'),selected+' input/read transport requests must not toggle lifecycle activity');
    for (const key of ['top','height','width','cols','rows']) assert.ok(
      Math.max(...samples.map(s=>s[key]))-Math.min(...samples.map(s=>s[key]))<0.5,
      selected+' terminal '+key+' changed while typing: '+JSON.stringify(samples));
    checks+=6;
  }
  await evaluate('f.replyDelay=0');
  await evaluate("window.statusMutations=0;window.statusObserver=new MutationObserver(records=>statusMutations+=records.length);statusObserver.observe(document.querySelector('.lc-terminal-status'),{childList:true,characterData:true,subtree:true})");
  await delay(300);assert.equal(await evaluate("statusMutations"),0,"Unchanged readiness must not flood screen-reader live regions");checks++;
  await evaluate("statusObserver.disconnect();f.transport.enqueue('x'.repeat(32769))");await delay(250);
  assert.equal(await evaluate("document.querySelector('.lc-terminal-status').textContent.includes('32 KiB')"),true);checks++;
  await evaluate("term.focus()");
  await command("Input.dispatchKeyEvent",{type:"keyDown",key:"l",code:"KeyL",text:"l"});
  await command("Input.dispatchKeyEvent",{type:"keyUp",key:"l",code:"KeyL"});
  await wait("f.sent.some(p=>p.action==='input' && p.bytes.includes(108))");
  assert.equal(await evaluate("f.keys"),0);checks++;
  await evaluate("f.append('\\x1b]2;BAD TITLE\\x07\\x1b]11;#ffffff\\x07\\x1b]52;c;c2VjcmV0\\x07\\x1b]8;;https://invalid.example/\\x07link\\x1b]8;;\\x07')");
  await wait("text().includes('link')");
  assert.equal(await evaluate("document.title"),"Terminal browser fixture");
  assert.equal(await evaluate("document.querySelector('.lc-terminal-viewport a') === null && [...document.querySelectorAll('#terminal a')].every(a=>a.classList.contains('lc-terminal-recovery') && a.pathname.startsWith('/runtime/runs/'))"),true);checks+=2;
  const palettes=[];
  for(const selected of ["light","dark","light"]) {
    await evaluate(`LineCableModelsTheme.select(${JSON.stringify(selected)})`);await delay(150);
    const colors=await evaluate("({theme:document.documentElement.dataset.lcmResolvedTheme,bg:term.options.theme.background,fg:term.options.theme.foreground,token:getComputedStyle(document.querySelector('#terminal')).getPropertyValue('--lc-console-bg').trim()})");
    assert.equal(colors.theme,selected);assert.equal(colors.bg,colors.token);assert.notEqual(colors.fg,colors.bg);checks+=3;palettes.push(colors);
    const shot=await command("Page.captureScreenshot",{format:"png"});await writeFile(join(scratch,"terminal-"+selected+".png"),Buffer.from(shot.data,"base64"));
  }
  assert.deepEqual(palettes[0],palettes[2]);assert.notDeepEqual(palettes[0],palettes[1]);checks+=2;
  for(const width of [1280,390,920]) {
    await command("Emulation.setDeviceMetricsOverride",{width,height:800,deviceScaleFactor:1,mobile:false});await delay(400);
    const geometry=await evaluate("({overflow:document.documentElement.scrollWidth>innerWidth,cols:term.cols,rows:term.rows,remote:f.sent.filter(p=>p.action==='resize').at(-1)})");
    assert.equal(geometry.overflow,false);assert.equal(geometry.remote.columns,geometry.cols);assert.equal(geometry.remote.rows,geometry.rows);checks+=3;
    assert.equal(await evaluate("document.querySelector('.xterm-screen').getBoundingClientRect().bottom <= document.querySelector('.lc-terminal-screen').getBoundingClientRect().bottom + 0.5"),true);checks++;
    const shot=await command("Page.captureScreenshot",{format:"png"});await writeFile(join(scratch,"terminal-"+width+".png"),Buffer.from(shot.data,"base64"));
  }
  await evaluate("f.button('Stop').click()");assert.equal(await evaluate("f.sent.some(p=>p.action==='stop')"),false);checks++;
  await evaluate("f.button('Cancel').click();f.button('Restart').click();f.button('Restart Julia').click()");
  await wait("f.sent.some(p=>p.action==='restart') && f.transport.state.phase==='starting'");
  await evaluate("f.phase='ready'");await wait("f.transport.state.canInput");
  await evaluate("f.uncertain=true;f.transport.enqueue('not_replayed()\\r')");await wait("f.transport.state.uncertain && !f.transport.busy");
  const inputs=await evaluate("f.sent.filter(p=>p.action==='input').length");await delay(300);
  assert.equal(await evaluate("f.sent.filter(p=>p.action==='input').length"),inputs);
  assert.equal(await evaluate("term.options.disableStdin && f.button('Interrupt').disabled"),true);checks+=2;
  await evaluate("f.button('Reconnect').click()");await wait("f.transport.state.connected && f.transport.state.reviewRequired && !f.button('Resume input').disabled");
  assert.equal(await evaluate("term.options.disableStdin"),true);checks++;
  await evaluate("f.button('Resume input').click()");await wait("f.transport.state.canInput");
  assert.equal(await evaluate("f.sent.filter(p=>p.action==='input').length"),inputs);checks++;
  await evaluate("document.querySelector('#terminal').remove()");await wait("f.transport.closed && f.transport.client.closed");
  const count=await evaluate("f.sent.length");await delay(300);
  assert.equal(await evaluate("f.sent.length"),count);assert.equal(await evaluate("f.sockets.every(s=>s.readyState===3)"),true);checks+=2;
  assert.deepEqual(errors,[]);checks++;
  console.log(`Private terminal browser: ${checks} assertions passed; real Chrome/xterm with explicit transport fixture. Screenshots: ${scratch}`);
  }
} finally {
  socket?.close();for(const item of pending.values())clearTimeout(item.timer);
  if(browser && browser.exitCode===null){browser.kill("SIGTERM");for(let i=0;i<100 && browser.exitCode===null;i++)await delay(50);if(browser.exitCode===null)browser.kill("SIGKILL");}
  server.closeAllConnections();await new Promise(resolve=>server.close(resolve));
}
