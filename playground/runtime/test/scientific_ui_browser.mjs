import assert from 'node:assert/strict';
import {writeFile} from 'node:fs/promises';

const [base, debug, directory] = process.argv.slice(2);
const delay = ms => new Promise(resolve => setTimeout(resolve, ms));
let runs;
for (const end = Date.now()+180000; Date.now()<end;) {
  runs = await (await fetch(base+'/runtime/api/runs')).json();
  assert(!runs.some(run=>run.state==='failed'), JSON.stringify(runs));
  if (runs.length===2 && runs.every(run=>run.state==='running')) break;
  await delay(250);
}
assert(runs?.length===2 && runs.every(run=>run.state==='running'), 'Registered UI hosts did not start');
const deckRun = runs.find(r=>r.application==='ichqp-showcase');
const studyRun = runs.find(r=>r.application==='cable-study');
assert(deckRun && studyRun);
const pages = await (await fetch(debug+'/json')).json();
const socket = new WebSocket(pages.find(p=>p.type==='page').webSocketDebuggerUrl);
await new Promise((resolve,reject)=>{socket.addEventListener('open',resolve,{once:true});socket.addEventListener('error',reject,{once:true});});
let sequence=0, checks=0;
const pending=new Map(), errors=[], requests=[];
socket.addEventListener('message',event=>{
  const message=JSON.parse(event.data);
  if (message.method==='Runtime.exceptionThrown') errors.push(message.params.exceptionDetails);
  if (message.method==='Network.requestWillBeSent') requests.push(message.params.request.url);
  const waiter=pending.get(message.id);
  if (!waiter) return;
  pending.delete(message.id); clearTimeout(waiter.timer);
  message.error ? waiter.reject(Error(message.error.message)) : waiter.resolve(message.result);
});
const command=(method,params={})=>new Promise((resolve,reject)=>{
  const id=++sequence;
  const timer=setTimeout(()=>{pending.delete(id);reject(Error('CDP timeout: '+method));},30000);
  pending.set(id,{resolve,reject,timer});
  socket.send(JSON.stringify({id,method,params}));
});
const evaluate=async expression=>{
  const result=await command('Runtime.evaluate',{expression,returnByValue:true,awaitPromise:true});
  if(result.exceptionDetails) throw Error(result.exceptionDetails.exception?.description ?? JSON.stringify(result.exceptionDetails));
  return result.result.value;
};
const check=async(expression,label)=>{assert(await evaluate(expression),label);checks++;};
const wait=async(expression,label)=>{
  for(const end=Date.now()+90000;Date.now()<end;){
    if(await evaluate(expression)) return;
    await delay(150);
  }
  throw Error(label);
};
const screenshot=async name=>{
  const result=await command('Page.captureScreenshot',{format:'png'});
  await writeFile(directory+'/'+name+'.png',Buffer.from(result.data,'base64'));
};
const ready=()=>wait("!!window.Reveal?.isReady() && !!document.querySelector('.lcm-deck-status')",'Reveal did not initialize');
const goto=route=>evaluate("Reveal.slide(Reveal.getSlides().findIndex(s=>s.querySelector('iframe[data-lcm-src=\""+route+"\"]')))");
try {
  await command('Page.enable'); await command('Runtime.enable'); await command('Network.enable');
  await command('Emulation.setDeviceMetricsOverride',{width:1440,height:900,deviceScaleFactor:1,mobile:false});
  await command('Page.navigate',{url:base+'/presentations/showcase.html'});
  await ready();
  await goto('/science/line-parameters');
  await wait("document.querySelector('iframe[data-lcm-src=\"/science/line-parameters\"]').dataset.lcmState==='unavailable'",'Public deck did not show immediate owned-run placeholder');
  await check("[...document.querySelectorAll('iframe[data-lcm-requires-run]')].every(f=>!f.hasAttribute('src'))",'Public deck allocated live frames');
  assert(!requests.some(url=>/\/science\//.test(url)),'Public deck made scientific UI requests'); checks++;
  await check("Reveal.getCurrentSlide().querySelector('.lcm-live-placeholder').textContent.includes('Launch this deck')",'Public fallback is not actionable');
  await screenshot('showcase-public');

  await command('Page.navigate',{url:base+'/presentations/showcase.html?lcm-run='+deckRun.id});
  await ready();
  await goto('/science/line-parameters');
  await wait("document.querySelector('iframe[data-lcm-src=\"/science/line-parameters\"]').dataset.lcmState==='ready'",'Owned scientific slide did not mount');
  await evaluate("window.live=document.querySelector('iframe[data-lcm-src=\"/science/line-parameters\"]'); window.sd=live.contentDocument; window.plot=sd.querySelector('.lc-study-plot');");
  await wait("!!sd.querySelector('[data-runtime-kind=execution]')",'Scientific job renderer did not mount');
  await check("sd.querySelectorAll('.lc-study-fields input').length===6 && !!plot",'Deck did not use the real typed scientific fields');
  await check("[...sd.querySelectorAll('.lc-study-fields input')].every(i=>i.checkValidity())",'A default numeric input violates its native min/step contract');
  await check("sd.body.textContent.includes('No calculation yet') && !sd.querySelector('.lc-study-curve').getAttribute('points')",'View implied a result before a calculation');
  await check("[...sd.querySelectorAll('button')].find(b=>b.textContent==='Run calculation')?.disabled",'Unassigned scientific view allowed execution');
  await check("(()=>{const r=[...sd.querySelectorAll('button')].find(b=>b.textContent==='Run calculation').getBoundingClientRect();return r.top>=0 && r.bottom<=live.clientHeight;})()",'Run action is outside the initial presentation viewport');
  await evaluate("window.change=(d,name,value)=>{const input=d.querySelector('[name=\"'+name+'\"]');input.value=value;input.dispatchEvent(new d.defaultView.Event('input',{bubbles:true}));}; change(sd,'frequency_points','2.5');");
  await wait("sd.body.textContent.includes('Complete all fields')",'Fractional sample count was not invalidated');
  await evaluate("change(sd,'frequency_points','');");
  await wait("sd.body.textContent.includes('Complete all fields')",'Blank numeric input silently reused its previous value');
  await evaluate("change(sd,'frequency_points','3');");
  await wait("!sd.body.textContent.includes('Complete all fields')",'Corrected numeric input did not recover');
  await check("plot===sd.querySelector('.lc-study-plot')",'Editing remounted the scientific plot');
  await check("plot.querySelectorAll('text').length===12 && [...plot.querySelectorAll('*')].every(n=>n.namespaceURI==='http://www.w3.org/2000/svg')",'Reactive labels escaped the native SVG namespace');
  await evaluate("Reveal.toggleOverview(true)");
  await check("plot===sd.querySelector('.lc-study-plot')",'Overview remounted a live scientific viewport');
  await evaluate("Reveal.toggleOverview(false)");
  await wait("live.dataset.lcmState==='ready'",'Leaving overview broke the live view');
  await screenshot('showcase-line');

  await evaluate("window.wbframe=document.createElement('iframe');wbframe.id='study-frame';wbframe.src="+JSON.stringify(base+'/applications/runs/'+studyRun.id+'/workbenches/cable-study')+";wbframe.style.cssText='position:fixed;inset:0;width:100vw;height:100vh;border:0;z-index:10000';document.body.append(wbframe)");
  await wait("!!wbframe.contentWindow.lcmXRay && !!wbframe.contentDocument.querySelector('.lc-study-runtime')",'Registered CableStudy did not mount');
  await wait("wbframe.contentWindow.WEBSOCKET?.isopen() && typeof [...wbframe.contentDocument.querySelectorAll('.lc-wb-nav-item')].find(b=>b.textContent.includes('Line parameters'))?.onclick==='function'",'Workbench action bindings were not ready');
  await evaluate("window.wd=wbframe.contentDocument;window.ww=wbframe.contentWindow;ww.lcmXRay.disable();[...wd.querySelectorAll('.lc-wb-nav-item')].find(b=>b.textContent.includes('Line parameters')).click();");
  await wait("wd.querySelector('[data-view=parameters]').classList.contains('is-active')",'Workbench navigation did not select the scientific view');
  await check("wd.querySelector('[data-view=parameters] .lc-study-fields [name=frequency_points]').valueAsNumber===40",'Draft inputs leaked between deck and workbench runs');
  await evaluate("window.wp=wd.querySelector('[data-view=parameters] .lc-study-plot');");
  for(const theme of ['light','dark','light','dark']){
    await evaluate("(()=>{const selector=wd.querySelector('[data-lc-wb-theme-selector]');selector.value="+JSON.stringify(theme)+";selector.dispatchEvent(new ww.Event('change'));})()");
    await wait("document.documentElement.dataset.lcmResolvedTheme==="+JSON.stringify(theme)+" && sd.documentElement.dataset.lcmResolvedTheme==="+JSON.stringify(theme)+" && wd.documentElement.dataset.lcmResolvedTheme==="+JSON.stringify(theme),'Theme did not propagate to all consumers');
    await check("(()=>{const colors=[sd,wd].map(d=>{const node=d===wd?d.querySelector('[data-view=parameters] .lc-study-view'):d.querySelector('.lc-study-view');const s=d.defaultView.getComputedStyle(node);return [s.color,s.backgroundColor];});return JSON.stringify(colors[0])===JSON.stringify(colors[1]) && colors[0][0]!==colors[0][1];})()",'Scientific styles drifted between deck and workbench');
    await check("plot===sd.querySelector('.lc-study-plot') && wp===wd.querySelector('[data-view=parameters] .lc-study-plot')",'Theme change remounted scientific view');
    await check("(()=>{const link=wd.querySelector('.lc-wb-sidebar a[href=\"/\"]');const probe=wd.createElement('span');probe.style.color='var(--lc-link)';link.append(probe);const expected=ww.getComputedStyle(probe).color;probe.remove();return ww.getComputedStyle(link).color===expected;})()",'Plain workbench footer link lost the shared theme colour');
  }
  await screenshot('cable-study-dark');
  await evaluate("ww.lcmXRay.enable(); wd.querySelector('[data-view=parameters] .lc-study-view').dispatchEvent(new ww.MouseEvent('click',{bubbles:true,composed:true}));");
  await wait("wd.querySelector('.lc-xray-host').shadowRoot.querySelector('.xray-body').textContent.includes('ScientificView')",'Scientific X-ray metadata did not render');
  await check("!wd.querySelector('.lc-xray-host').shadowRoot.querySelector('.xray-body').textContent.includes("+JSON.stringify(studyRun.id)+")",'Private run identity leaked into scientific metadata');
  await evaluate("ww.lcmXRay.disable();wbframe.remove();");
  await goto('/science/geometry');
  await wait("document.querySelector('iframe[data-lcm-src=\"/science/geometry\"]').dataset.lcmState==='ready'",'Geometry slide did not mount');
  await evaluate("window.gd=document.querySelector('iframe[data-lcm-src=\"/science/geometry\"]').contentDocument;");
  await wait("gd.defaultView.WEBSOCKET?.isopen() && typeof gd.querySelector('[name=core]').oninput==='function'",'Geometry input binding was not ready');
  await evaluate("window.circle=gd.querySelector('.lc-study-core');window.oldRadius=circle.getAttribute('r');change(gd,'core','20')");
  await wait("circle.getAttribute('r')!==oldRadius",'Geometry input did not update its existing SVG circle');
  await check("Math.abs(Number(circle.getAttribute('r'))-140*20/26)<1e-6",'Geometry display did not preserve radial proportions');
  await check("circle===gd.querySelector('.lc-study-core')",'Geometry change remounted the canvas');
  await goto('/science/corridor');
  await wait("document.querySelector('iframe[data-lcm-src=\"/science/corridor\"]').dataset.lcmState==='ready'",'Corridor slide did not mount');
  await check("document.querySelector('iframe[data-lcm-src=\"/science/corridor\"]').contentDocument.body.textContent.includes('not confidence intervals')",'Corridor view lost scientific assumptions');
  await goto('/science/terminal');
  await wait("document.querySelector('iframe[data-lcm-src=\"/science/terminal\"]').dataset.lcmState==='ready'",'Actual registered deck terminal did not mount');
  await check("!!document.querySelector('iframe[data-lcm-src=\"/science/terminal\"]').contentDocument.querySelector('.xterm')",'Registered deck is not using shared xterm renderer');
  await goto('/science/runtime');
  await wait("document.querySelector('iframe[data-lcm-src=\"/science/runtime\"]').dataset.lcmState==='ready'",'Preparation slide did not mount');
  await check("[...document.querySelectorAll('iframe[data-lcm-src]')].every(f=>f.getAttribute('src').startsWith('/applications/runs/"+deckRun.id+"/'))",'Live deck frames lost their single owned run');

  await command('Page.navigate',{url:base+'/presentations/showcase.html?lcm-print'});
  await ready();
  await check("[...document.querySelectorAll('iframe[data-lcm-src]')].every(f=>!f.hasAttribute('src'))",'Print view loaded live frames');
  const pdf=await command('Page.printToPDF',{printBackground:true,preferCSSPageSize:true,displayHeaderFooter:false});
  await writeFile(directory+'/showcase.pdf',Buffer.from(pdf.data,'base64'));
  assert.equal(errors.length,0,JSON.stringify(errors));checks++;
  console.log('PASS: '+checks+' actual registered consumer browser checks; no scientific execution or terminal-isolation claim');
} catch(error) {
  await screenshot('scientific-ui-failure').catch(()=>{});
  const state=await evaluate("(()=>{const d=window.wbframe?.contentDocument;return {workbenchReady:d?.readyState,socketOpen:window.wbframe?.contentWindow?.WEBSOCKET?.isopen(),navigation:d?[...d.querySelectorAll('.lc-wb-nav-item')].map(n=>({text:n.textContent,handler:typeof n.onclick,current:n.getAttribute('aria-current')})):[],geometry:window.gd?{value:gd.querySelector('[name=core]')?.value,handler:typeof gd.querySelector('[name=core]')?.oninput,radius:gd.querySelector('.lc-study-core')?.getAttribute('r')}:null};})()").catch(()=>null);
  await writeFile(directory+'/scientific-ui-failure.json',JSON.stringify({error:String(error),state,errors,requests},null,2));
  throw error;
} finally {
  for(const item of pending.values()){clearTimeout(item.timer);item.reject(Error('test closing'));}
  pending.clear();socket.close();
}
