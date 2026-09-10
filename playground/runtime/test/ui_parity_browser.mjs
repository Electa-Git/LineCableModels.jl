import assert from 'node:assert/strict';
import {writeFile} from 'node:fs/promises';
import {assertRuntimeActivity} from './runtime_activity_browser.mjs';
const [base, debug, directory] = process.argv.slice(2);
const delay = ms => new Promise(resolve => setTimeout(resolve, ms));
let runs;
// The fixture starts four hosts sequentially, each with a bounded 180 s cold
// startup. Its aggregate wait must cover those windows; a failed host still
// fails immediately below rather than consuming the full allowance.
for (const end = Date.now() + (4 * 180 + 30) * 1000; Date.now() < end;) {
  runs = await (await fetch(base + '/runtime/api/runs')).json();
  assert(runs.every(run => ['reserved', 'starting', 'running'].includes(run.state)), JSON.stringify(runs));
  if (runs.length === 4 && runs.every(run => run.state === 'running')) break;
  await delay(200);
}
assert(runs.length === 4 && runs.every(run => run.state === 'running'), 'UI hosts did not start');
const route = (app, path) => base + '/applications/runs/' + runs.find(run => run.application === app).id + path;
const pages = await (await fetch(debug + '/json')).json();
const socket = new WebSocket(pages.find(page => page.type === 'page').webSocketDebuggerUrl);
await new Promise(resolve => socket.addEventListener('open', resolve, {once: true}));
let sequence = 0, checks = 0;
const pending = new Map(), errors = [];
socket.addEventListener('message', event => {
  const message = JSON.parse(event.data), waiter = pending.get(message.id);
  if (message.method === 'Runtime.exceptionThrown') errors.push(message.params.exceptionDetails);
  if (!waiter) return;
  pending.delete(message.id); clearTimeout(waiter.timer);
  message.error ? waiter.reject(Error(message.error.message)) : waiter.resolve(message.result);
});
const command = (method, params = {}) => new Promise((resolve, reject) => {
  const id = ++sequence, timer = setTimeout(() => { pending.delete(id); reject(Error(method + ' timed out')); }, 30000);
  pending.set(id, {resolve, reject, timer}); socket.send(JSON.stringify({id, method, params}));
});
const read = async expression => {
  const result = await command('Runtime.evaluate', {expression, returnByValue: true, awaitPromise: true});
  if (result.exceptionDetails) throw Error(result.exceptionDetails.exception?.description ?? result.exceptionDetails.text);
  return result.result.value;
};
const wait = async (expression, message) => {
  for (const end = Date.now() + 90000; Date.now() < end;) {
    if (await read(expression)) return;
    await delay(100);
  }
  throw Error(message);
};
const viewport = (width, height) => command('Emulation.setDeviceMetricsOverride', {width, height, deviceScaleFactor: 1, mobile: false});
async function open(url, selector) {
  await command('Page.navigate', {url});
  await wait(`location.pathname === ${JSON.stringify(new URL(url).pathname)} && document.readyState === 'complete' &&
    !!document.querySelector(${JSON.stringify(selector)})`, 'Page did not mount: ' + url);
  await read('document.fonts.ready');
  await wait('window.WEBSOCKET?.isopen() === true', 'Bonito websocket not ready');
}
const theme = async value => {
  await read(`LineCableModelsTheme.select(${JSON.stringify(value)})`);
  await wait(`document.documentElement.dataset.lcmResolvedTheme === ${JSON.stringify(value)}`, 'Theme did not propagate');
};
const styles = selector => read(`(() => {
  const node = document.querySelector(${JSON.stringify(selector)}); if (!node) throw Error('Missing ' + ${JSON.stringify(selector)});
  const css = getComputedStyle(node);
  return Object.fromEntries(['fontFamily','fontSize','fontWeight','lineHeight','paddingTop','paddingRight','paddingBottom',
    'paddingLeft','borderTopWidth','borderRadius','color','backgroundColor'].map(key => [key, css[key]]));
})()`);
const equal = (actual, expected, label) => { assert.deepEqual(actual, expected, label); checks++; };
const check = async (expression, label) => { assert(await read(expression), label); checks++; };
const shot = async name => {
  const result = await command('Page.captureScreenshot', {format: 'png'});
  await writeFile(directory + '/' + name + '.png', Buffer.from(result.data, 'base64'));
};
async function selectView(label) {
  await wait(`typeof [...document.querySelectorAll('.lc-wb-nav-item')].find(node => node.textContent.includes(${JSON.stringify(label)}))?.onclick === 'function'`, 'Navigation action not ready');
  await read(`[...document.querySelectorAll('.lc-wb-nav-item')].find(node => node.textContent.includes(${JSON.stringify(label)})).click()`);
  await wait(`(() => {
    const item=[...document.querySelectorAll('.lc-wb-nav-item')].find(node => node.textContent.includes(${JSON.stringify(label)}));
    return item?.getAttribute('aria-current') === 'page' &&
      document.querySelector('.lc-wb-view.is-active')?.getAttribute('aria-label') === ${JSON.stringify(label)};
  })()`, 'View selection did not settle: ' + label);
}
const visible = '.lc-wb-view.is-active ';
async function dockInsets(mode) {
  const labels = await read(`[...document.querySelectorAll('.lc-wb-dock-tab')].map(n=>n.textContent)`);
  for (const label of labels) {
    await wait(`typeof [...document.querySelectorAll('.lc-wb-dock-tab')].find(n=>n.textContent===${JSON.stringify(label)})?.onclick==='function'`, 'Dock tab not ready');
    await read(`[...document.querySelectorAll('.lc-wb-dock-tab')].find(n=>n.textContent===${JSON.stringify(label)}).click()`);
    await wait(`document.querySelector('.lc-wb-dock-panel.is-active')?.getAttribute('aria-label')===${JSON.stringify(label)}`, 'Dock tab did not activate');
    await check(`(() => {
      const panel=document.querySelector('.lc-wb-dock-panel.is-active');
      const child=[...panel.querySelectorAll('*')].find(n=>n.getBoundingClientRect().width>0 && getComputedStyle(n).display!=='contents');
      const p=panel.getBoundingClientRect(), c=child.getBoundingClientRect(), style=getComputedStyle(panel);
      const inset=parseFloat(style.paddingTop);
      return inset>=10 && ['paddingLeft','paddingRight','paddingBottom'].every(k=>parseFloat(style[k])===inset) &&
        c.top>=p.top+inset-1 && c.left>=p.left+inset-1 && c.right<=p.right-inset+1;
    })()`, mode+' shared dock inset: '+label);
  }
}
try {
  await command('Page.enable'); await command('Runtime.enable'); await command('Network.enable');
  await command('Network.setCacheDisabled', {cacheDisabled: true});
  for (const mode of ['dark', 'light']) {
    await viewport(1440, 900);
    await open(route('toolkit-gallery', '/widgets/form-toolkit'), '.lc-form'); await theme(mode);
    const number = await styles('.lc-unit-number input');
    const choice = await styles('select[name=earth_model]');
    const button = await styles('.lc-button-secondary');
    const fieldLabel = await styles('.lc-field-label');
    await open(route('toolkit-gallery', '/widgets/overlay-toolkit'), '.lc-feedback-specimen'); await theme(mode);
    await check(`document.querySelector('[data-lcm-component="StatusIndicator"]') !== null || document.querySelectorAll('.lc-status-indicator').length >= 4`, mode+' shared status specimen');
    const statusSuccess = await read(`getComputedStyle(document.querySelector('.lc-status-indicator[data-tone="success"]')).color`);
    await wait(`typeof [...document.querySelectorAll('button')].find(b=>b.textContent==='Start feedback preview')?.onclick === 'function'`, 'Feedback example did not bind');
    await read(`[...document.querySelectorAll('button')].find(b=>b.textContent==='Start feedback preview').click()`);
    await wait(`[...document.querySelectorAll('button')].some(b=>b.textContent==='Preview active…' && b.disabled && b.getAttribute('aria-busy')==='true')`, 'ActionButton did not display owner-supplied activity');
    await check(`[...document.querySelectorAll('.lc-status-indicator')].some(n=>n.textContent==='Preview active · no worker contacted' && n.dataset.tone==='info' && n.dataset.busy==='true')`, mode+' StatusIndicator attributes follow activity');
    await read(`[...document.querySelectorAll('button')].find(b=>b.textContent==='End feedback preview').click()`);
    await wait(`[...document.querySelectorAll('button')].some(b=>b.textContent==='Start feedback preview' && !b.disabled)`, 'ActionButton did not reset');
    await check(`[...document.querySelectorAll('.lc-status-indicator')].some(n=>n.textContent==='Feedback preview complete' && n.dataset.tone==='success' && n.dataset.busy==='false')`, mode+' StatusIndicator attributes follow completion');
    await open(route('template-workbench', '/workbenches/template'), '.lc-wb-shell'); await theme(mode);
    await selectView('Cable geometry');
    equal(await styles(visible + '.lc-unit-number input'), number, mode + ' gallery/template numeric field');
    equal(await styles(visible + 'select[name=earth_model]'), choice, mode + ' gallery/template choice');
    equal(await styles(visible + '.lc-field-label'), fieldLabel, mode + ' gallery/template label');
    const heading = await styles(visible + '.lc-workspace-header h1');
    const page = await styles(visible + '.lc-workspace-page');
    const navigation = await styles('.lc-wb-nav-item[aria-current="page"]');
    await dockInsets(mode);
    await open(route('cable-study', '/workbenches/cable-study'), '.lc-wb-shell'); await theme(mode);
    await wait(`document.querySelectorAll('[data-runtime-kind="diagnostics"]').length === 1`, 'CableStudy must own exactly one diagnostics view');
    await wait(`document.querySelector('.lc-runtime-run-status .lc-status-indicator[data-tone="success"]') !== null`, 'Owned run status did not load');
    equal(await read(`getComputedStyle(document.querySelector('.lc-runtime-run-status .lc-status-indicator[data-tone="success"]')).color`), statusSuccess, mode+' gallery/live status colour');
    await selectView('Cable construction');
    await dockInsets(mode);
    equal(await styles(visible + '.lc-workspace-header h1'), heading, mode + ' template/live heading');
    equal(await styles(visible + '.lc-workspace-page'), page, mode + ' template/live page spacing');
    equal(await styles('.lc-wb-nav-item[aria-current="page"]'), navigation, mode + ' template/live navigation');
    await selectView('Line parameters');
    equal(await styles(visible + '.lc-unit-number input'), number, mode + ' gallery/live numeric field');
    equal(await styles(visible + 'select[name=quantity]'), choice, mode + ' gallery/live choice');
    equal(await styles(visible + '.lc-button-secondary'), button, mode + ' gallery/live action');
    equal(await styles(visible + '.lc-field-label'), fieldLabel, mode + ' gallery/live label');
    await check(`document.querySelector('${visible}.lc-wb-split') !== null`, 'Live view bypassed SplitPane');
    await check(`document.querySelector('${visible}.lc-study-plot') !== null`, 'Scientific plot did not mount');
    await shot('workbench-' + mode);
    await read(`window.keptPlot=document.querySelector('${visible}.lc-study-plot')`);
    const grip = await read(`(() => {const r=document.querySelector('${visible}.lc-wb-splitter-handle').getBoundingClientRect();return {x:r.x+r.width/2,y:r.y+r.height/2};})()`);
    const firstWidth = await read(`document.querySelector('${visible}.lc-wb-split-first').getBoundingClientRect().width`);
    await read(`window.splitDebug = {hit:document.elementFromPoint(${grip.x},${grip.y})?.outerHTML,
      before:document.querySelector('${visible}.lc-wb-split').style.cssText}`);
    await command('Input.dispatchMouseEvent', {type:'mouseMoved', ...grip});
    await command('Input.dispatchMouseEvent', {type:'mousePressed', ...grip, button:'left', buttons:1, clickCount:1});
    await read(`splitDebug.pressed=document.querySelector('${visible}.lc-wb-split').className`);
    for (const offset of [20,40,60,80]) {
      await command('Input.dispatchMouseEvent', {type:'mouseMoved', x:grip.x-offset, y:grip.y, button:'left', buttons:1});
      await delay(30);
    }
    await command('Input.dispatchMouseEvent', {type:'mouseReleased', x:grip.x-80, y:grip.y, button:'left', clickCount:1});
    await check(`Math.abs(document.querySelector('${visible}.lc-wb-split-first').getBoundingClientRect().width - ${firstWidth}) > 30`,
      'Scientific split does not resize: ' + JSON.stringify(await read(`({...splitDebug, after:document.querySelector('${visible}.lc-wb-split').style.cssText})`)) + JSON.stringify(errors));
    await check(`keptPlot === document.querySelector('${visible}.lc-study-plot')`, 'Resizing remounted the plot');
    for (const width of [1024,768,390]) {
      await viewport(width, 800); await delay(150);
      await check('document.documentElement.scrollWidth <= innerWidth + 1', mode + ' workbench page overflows at ' + width);
      await check(`(() => {const pane=document.querySelector('${visible}.lc-wb-split'); return pane.clientWidth >= pane.scrollWidth-1;})()`, mode + ' split overflows at ' + width);
    }
    await viewport(1440,900);
    await open(route('ichqp-showcase', '/science/line-parameters'), '.lc-study-view'); await theme(mode);
    equal(await styles('.lc-unit-number input'), number, mode + ' gallery/presentation numeric field');
    equal(await styles('select[name=quantity]'), choice, mode + ' gallery/presentation choice');
    equal(await styles('.lc-button-secondary'), button, mode + ' gallery/presentation action');
    await shot('scientific-frame-' + mode);
    await check('!!document.querySelector(".lc-workspace-page .lc-wb-split")', 'Presentation does not share the workspace composition');
    await open(route('toolkit-gallery', '/widgets/runtime-controls'), '.lc-runtime-controls'); await theme(mode);
    equal(await styles('.lc-button-secondary'), button, mode + ' runtime gallery/action');
    equal(await styles('.lc-field-label'), fieldLabel, mode + ' runtime gallery/label');
  }
  await assertRuntimeActivity({command,read,wait,base,run:runs[0],shot});
  console.log('Application parity: ' + checks + ' checks across actual gallery, template, workbench and scientific frame; both themes');
  assert.equal(errors.length, 0, JSON.stringify(errors));
} finally { socket.close(); }
