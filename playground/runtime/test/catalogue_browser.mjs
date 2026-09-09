import assert from 'node:assert/strict';
import {assertPublishedShell} from '../../test/integration/published_shell_browser.mjs';
const [base, debug] = process.argv.slice(2);
const pages=await (await fetch(debug+'/json')).json();
const socket=new WebSocket(pages.find(p=>p.type==='page').webSocketDebuggerUrl);
await new Promise(resolve=>socket.addEventListener('open',resolve,{once:true}));
let sequence=0;
const pending=new Map(), errors=[];
socket.addEventListener('message',event=>{
  const message=JSON.parse(event.data);
  if(message.method==='Runtime.exceptionThrown') errors.push(message.params.exceptionDetails);
  const request=pending.get(message.id); if(!request) return; pending.delete(message.id);
  message.error ? request.reject(Error(message.error.message)) : request.resolve(message.result);
});
const command=(method,params={})=>new Promise((resolve,reject)=>{
  const id=++sequence; pending.set(id,{resolve,reject});socket.send(JSON.stringify({id,method,params}));
});
const evaluate=async expression=>{
  const result=await command('Runtime.evaluate',{expression,returnByValue:true,awaitPromise:true});
  if(result.exceptionDetails) throw Error(result.exceptionDetails.exception?.description);
  return result.result.value;
};
const wait=async(expression,message)=>{
  const deadline=Date.now()+150000;
  while(Date.now()<deadline){if(await evaluate(expression))return;await new Promise(r=>setTimeout(r,100));}
  throw Error(message+'\n'+JSON.stringify(errors));
};
const navigate=async url=>{
  await command('Page.navigate',{url});
  const expected=new URL(url);
  // Reveal initializes its own title-slide hash; that is not a navigation
  // failure when the requested URL did not name a particular slide.
  await wait(`location.pathname===${JSON.stringify(expected.pathname)} && location.search===${JSON.stringify(expected.search)} &&
    (${JSON.stringify(expected.hash)}==='' || location.hash===${JSON.stringify(expected.hash)}) && document.readyState==='complete'`, 'Page did not load: '+url);
};
const viewport=async(width,height)=>{
  await command('Emulation.setDeviceMetricsOverride',{width,height,deviceScaleFactor:1,mobile:false});
  await evaluate('new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))');
};
const ownedRuns=async()=>{const response=await fetch(base+'/runtime/api/runs');assert.equal(response.status,200);return response.json();};
const createRun=async application=>{
  const response=await fetch(base+'/runtime/api/runs',{method:'POST',headers:{Origin:base,'X-LCM-Request':'1','Content-Type':'application/json'},
    body:JSON.stringify({application,request_id:crypto.randomUUID()})});
  assert.equal(response.status,202); return response.json();
};
try {
  await command('Page.enable');await command('Runtime.enable');
  assert.deepEqual(await ownedRuns(),[]);
  const catalogue=await (await fetch(base+'/assets/application-catalogue.json')).json();
  assert.equal(new Set(catalogue.map(entry=>entry.id)).size,catalogue.length);
  assert(catalogue.every(entry=>!('ui' in entry)));
  const publishedPublic=catalogue.filter(entry=>entry.visibility==='public');
  const registeredPublic=await (await fetch(base+'/runtime/api/applications')).json();
  const identities=entries=>entries.map(({id,entrypoint,version})=>({id,entrypoint,version})).sort((a,b)=>a.id.localeCompare(b.id));
  assert.deepEqual(identities(publishedPublic),identities(registeredPublic),'Published and registered public catalogues diverged');
  assert(publishedPublic.some(entry=>entry.id==='ichqp-showcase'));
  assert(publishedPublic.some(entry=>entry.id==='cable-study'));
  await assertPublishedShell({devtools:{command},baseUrl:base,
    evaluate:(_,expression)=>evaluate(expression),navigate:(_,url)=>navigate(url),
    waitUntil:(_,expression,message)=>wait(expression,message),setViewport:(_,w,h)=>viewport(w,h),assert});
  await viewport(1600,1000);
  const links=new Set();
  for(const path of ['/','/dev/','/presentations/','/workbenches/','/dev/workbench.html','/dev/presentations.html']) {
    await navigate(base+path);
    const targets=await evaluate(`[...document.querySelectorAll('a[href]')].map(a=>new URL(a.href))
      .filter(u=>u.origin===location.origin && (u.pathname.endsWith('.html')||u.pathname.endsWith('/')))
      .map(u=>u.pathname)`);
    targets.forEach(target=>links.add(target));
  }
  for(const path of links) assert.equal((await fetch(base+path)).status,200,'Broken published link: '+path);
  assert.deepEqual(await ownedRuns(),[], 'Published pages allocated runtime resources');
  await navigate(base+'/presentations/');
  await wait("!!document.querySelector('[data-lcm-application-catalogue] select')",'Public deck selector missing');
  await wait("document.querySelector('[data-lcm-application-catalogue] select').value==='ichqp-showcase' && !document.querySelector('[data-lcm-application-catalogue] button').disabled",
    'Registered scientific UI must be launchable independently of worker availability');
  await navigate(base+'/presentations/showcase.html');
  await wait("document.documentElement.dataset.lcmDeckReady==='true'",'Scientific narrative is not available as a static deck');
  assert.deepEqual(await ownedRuns(),[],'Opening static slides allocated a run');
  await navigate(base+'/dev/');
  await wait("!!document.querySelector('[data-kind=\"presentation\"] select') && !document.querySelector('[data-kind=\"presentation\"] button').disabled",'Installed developer decks are not launchable');
  await evaluate(`const chooser=document.querySelector('[data-kind="presentation"]');
    chooser.querySelector('select').value='starter-deck'; chooser.querySelector('select').dispatchEvent(new Event('change'));
    chooser.querySelector('button').click();`);
  await wait("location.pathname==='/presentations/starter.html' && new URLSearchParams(location.search).has('lcm-run') && document.documentElement.dataset.lcmDeckReady==='true'",'Explicit launch did not open the selected deck');
  const runId=await evaluate("new URLSearchParams(location.search).get('lcm-run')");
  assert.equal((await ownedRuns()).length,1);
  await evaluate("window.liveIndex=[...document.querySelectorAll('.reveal .slides > section')].findIndex(s=>s.querySelector('iframe')); Reveal.slide(liveIndex)");
  await wait("document.querySelector('iframe[data-lcm-src]')?.dataset.lcmState==='ready'",'Owned deck frame did not become ready');
  const source=await evaluate("document.querySelector('iframe[data-lcm-src]').getAttribute('src')");
  assert.equal(source,'/applications/runs/'+runId+'/widgets/control-panel');
  await evaluate("window.keptFrame=document.querySelector('iframe[data-lcm-src]');window.keptDocument=keptFrame.contentDocument");
  for(const theme of ['light','dark']) {
    await evaluate(`LineCableModelsTheme.select(${JSON.stringify(theme)}); Reveal.slide(0); Reveal.slide(liveIndex)`);
    await wait(`keptFrame.contentDocument.documentElement.dataset.lcmResolvedTheme===${JSON.stringify(theme)}`,'Theme crossed the deck boundary incorrectly');
    assert(await evaluate("keptFrame===document.querySelector('iframe[data-lcm-src]') && keptDocument===keptFrame.contentDocument"),'Slide/theme switch remounted a live frame');
  }
  const workbench=await createRun('template-workbench');
  await navigate(base+'/runtime/runs/'+workbench.id);
  await wait(`location.pathname==='/applications/runs/${workbench.id}/workbenches/template' && !!window.lcmXRay`,'Actual catalogue workbench driver did not start');
  assert(await evaluate("!!document.querySelector('.lc-wb-shell') && !!document.querySelector('[data-lc-wb-theme-selector]')"),'Workbench lost its reusable shell');
  const gallery=await createRun('toolkit-gallery');
  await navigate(base+'/runtime/runs/'+gallery.id);
  await wait(`location.pathname==='/widgets/index.html' && new URLSearchParams(location.search).get('lcm-run')==='${gallery.id}'`,
    'Owned widget gallery did not open its published entry surface');
  await wait("document.readyState==='complete' && !!document.querySelector('iframe[data-lc-published-src=\"/widgets/file-upload\"]') && [...document.querySelectorAll('iframe[data-lc-published-src]')].filter(f=>!f.hidden).every(f=>f.getAttribute('src')?.startsWith('/applications/runs/'))",'Gallery frames bypassed their run namespace');
  await evaluate("window.uploadFrame=document.querySelector('iframe[data-lc-published-src=\"/widgets/file-upload\"]'); uploadFrame.scrollIntoView({block:'center'})");
  await wait("uploadFrame.contentDocument?.querySelector('.lc-upload-field')?.dataset.uploadReady==='true'",'Original upload widget did not initialize in its owned gallery');
  await evaluate(`window.uploadWindow=uploadFrame.contentWindow;
    window.uploadField=uploadFrame.contentDocument.querySelector('.lc-upload-field');
    window.uploadInput=uploadField.querySelector('input[type=file]');
    const transfer=new uploadWindow.DataTransfer();
    transfer.items.add(new uploadWindow.File(['owned upload fixture'], 'fixture.txt', {type:'text/plain'}));
    uploadInput.files=transfer.files; uploadInput.dispatchEvent(new uploadWindow.Event('change',{bubbles:true}));`);
  await wait(`(() => {
    if (uploadField.dataset.uploadState==='failed') throw Error(uploadField.querySelector('.lc-upload-error').textContent);
    return uploadField.dataset.uploadState==='ready';
  })()`, 'Owned upload POST failed through the proxy');
  const uploadURL=await evaluate('uploadField.dataset.uploadUrl');
  assert(uploadURL.startsWith('/applications/runs/'+gallery.id+'/uploads/'));
  assert.equal(await (await fetch(base+uploadURL)).text(),'owned upload fixture');
  await evaluate("uploadField.querySelector('.lc-upload-remove').click()");
  await wait("uploadField.dataset.uploadState==='empty'",'Owned upload removal failed');
  assert.equal((await fetch(base+uploadURL)).status,404);
  assert.equal(errors.length,0,JSON.stringify(errors));
  console.log('PASS: actual lcm runtime launcher, public/developer shell and links, passive catalogues, owned deck/workbench/gallery, native upload and removal');
  // Keep the application socket open for the harness shutdown test.
} finally {socket.close();}
