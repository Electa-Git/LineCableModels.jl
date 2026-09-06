#!/usr/bin/env node
import assert from 'node:assert/strict';

const base = process.argv[2] ?? 'http://127.0.0.1:18106';
const debug = process.argv[3] ?? 'http://127.0.0.1:19346';
const pages = await (await fetch(debug + '/json')).json();
const socket = new WebSocket(pages.find(page => page.type === 'page').webSocketDebuggerUrl);
await new Promise(resolve => socket.addEventListener('open', resolve, {once: true}));
let sequence = 0;
const pending = new Map(), errors = [];
socket.addEventListener('message', event => {
  const message = JSON.parse(event.data);
  if (message.method === 'Runtime.exceptionThrown') errors.push(message.params.exceptionDetails.text);
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
  const deadline = Date.now() + 45000;
  while (Date.now() < deadline) {
    if (await evaluate(expression)) return;
    await new Promise(resolve => setTimeout(resolve, 50));
  }
  throw Error(message);
};
const resize = async (width = 1600, height = 1000) => {
  await command('Emulation.setDeviceMetricsOverride', {width, height, deviceScaleFactor: 1, mobile: false});
  // CDP acknowledges metrics before the resize event and layout observers run.
  await evaluate('new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(() => resolve(true))))');
};
const move = (x, y) => command('Input.dispatchMouseEvent', {type: 'mouseMoved', x, y});
const point = expression => evaluate(`(() => { const r = (${expression}).getBoundingClientRect();
  return {x:r.x + Math.min(15,r.width/2), y:r.y + Math.min(15,r.height/2)}; })()`);
const hover = async expression => { const p = await point(expression); await move(p.x, p.y); };
const click = async expression => {
  const p = await point(expression); await move(p.x, p.y);
  await command('Input.dispatchMouseEvent', {type:'mousePressed', button:'left', clickCount:1, ...p});
  await command('Input.dispatchMouseEvent', {type:'mouseReleased', button:'left', clickCount:1, ...p});
};
const check = async (expression, message) => assert(await evaluate(expression), message);
const action = label => evaluate(`action(${JSON.stringify(label)}).click()`);
const theme = value => evaluate(`(() => { const input = document.querySelector('[data-lc-wb-theme-selector]');
  input.value = ${JSON.stringify(value)}; input.dispatchEvent(new Event('change')); })()`);
const open = async (route = '/fixture/xray') => {
  await command('Page.navigate', {url: base + route});
  await wait(`location.pathname === ${JSON.stringify(route)} && !!window.lcmXRay &&
    !!document.querySelector('[data-lcm-component]')`, 'X-ray did not mount: ' + route);
  if (route.startsWith('/fixture/')) await wait('!!window.__xrayFixtureReady', 'Bonito fixture handlers not ready');
  await evaluate(`window.s = document.querySelector('.lc-xray-host').shadowRoot;
    window.action = label => [...s.querySelectorAll('button')].find(n => n.textContent === label);
    window.rule = (selector, condition = '') => [...s.querySelectorAll('.xray-css-rule')].find(n =>
      n.querySelector('.xray-selector').textContent === selector &&
      (n.querySelector('.xray-css-condition')?.textContent || '') === condition);
    window.field = (property, selector = '.xp-card', condition = '') =>
      rule(selector, condition).querySelector('[data-css-property="' + property + '"]');
    window.enter = (property, value, selector = '.xp-card', condition = '') => {
      const row = field(property, selector, condition);
      const input = row.querySelector('[aria-label="' + property + ' preview"]');
      input.value = value; input.dispatchEvent(new Event(input.tagName === 'SELECT' ? 'change' : 'input'));
      return !row.classList.contains('is-invalid');
    };
    window.styleA = () => getComputedStyle(document.querySelector('#card-a'));
    window.styleB = () => getComputedStyle(document.querySelector('#card-b'));
    window.originalStyles = [...document.querySelectorAll('style[data-lcm-css-source]')].map(n => n.textContent);
    window.rootIdentity = lcmXRay.root;
  `);
};
const resetCheck = async () => check(`lcmXRay.preview.changes().length === 0 &&
  !document.querySelector('[data-lcm-preview-scope]') &&
  ![...document.styleSheets].some(sheet => { try { return [...sheet.cssRules].some(rule => rule.cssText.includes('data-lcm-preview-scope')); }
    catch (_) { return false; } })`, 'temporary declarations or scope markers leaked');

try {
  await command('Page.enable'); await command('Runtime.enable');
  await resize(); await open(); await theme('dark');
  await hover("document.querySelector('#card-a')");
  await wait("s.querySelector('.xray-outline').classList.contains('is-visible')", 'hover outline missing');
  await check("s.querySelector('.xray-panel').hidden", 'hover must not open inspection');
  await click("document.querySelector('#card-a')");
  await evaluate(`window.authoredRules = LineCableModelsCSSPreview.collect(document.querySelector('#card-a'), ['.xp-card'])
    .flatMap(record => (record.copies || []).map(copy => ({rule:copy.rule, text:copy.rule.style.cssText})))`);
  await check("!s.querySelector('.xray-panel').hidden && s.querySelector('.xray-body').textContent.includes('card-a')",
    'click must select the component');
  await check("!rule('.xp-cardinality') && [...s.querySelectorAll('.xray-css-rule')].length === 7",
    'ownership must use exact class tokens and deduplicate repeated stylesheets');
  await check(`field('padding').querySelector('input[type=number]') && field('display').querySelector('select') &&
    field('--application-value','element.style').classList.contains('is-readonly') &&
    ![...s.querySelectorAll('.xray-section')].filter(n=>!n.querySelector('.xray-css-rule')).some(n=>n.querySelector('input,select'))`,
    'typed CSS controls / read-only code boundary');
  await check("![...s.querySelectorAll('.xray-authored')].some(n=>/:\\s*$/.test(n.textContent))", 'unavailable authored values must not become blank editable defaults');
  await check("enter('padding', '28') && styleA().paddingTop === '28px' && styleB().paddingTop === '12px'",
    'preview must affect only selected instance, even with duplicate stylesheets');
  await evaluate("window.keptInput = field('padding').querySelector('input[type=number]'); keptInput.focus()");
  // Right-hand card can be partly behind the movable inspector: use its left edge.
  await hover("document.querySelector('#card-b')");
  await check(`field('padding').querySelector('input[type=number]') === keptInput && keptInput.value === '28' &&
    s.querySelector('.xray-body').textContent.includes('card-a')`, 'hover replaced the selected inspector or draft');

  // Picker consumes application clicks. Interact leaves the selection pinned.
  await click("document.querySelector('#card-a button')");
  await check("document.querySelector('#card-a output').textContent === '0'", 'pick accidentally invoked a callback');
  await action('Pick components');
  await click("document.querySelector('#card-a button')");
  await wait("document.querySelector('#card-a output').textContent === '1'", 'Interact blocked application callback');
  await wait("[...s.querySelectorAll('.xray-section')].find(n=>n.querySelector('h3').textContent==='Bindings').textContent.includes('1')",
    'read-only binding did not update');
  await check("field('padding').querySelector('input[type=number]') === keptInput && keptInput.value === '28'",
    'live binding update replaced draft controls');
  await action('Interact with application');

  await check("enter('gap', '34', '.xp-card', '@media (min-width: 900px)') && styleA().gap === '34px' && styleB().gap === '24px'",
    'conditional preview lost cascade/source order');
  await check("!enter('gap', '90', '.xp-card', '@media (min-width: 900px)') && styleA().gap === '34px'",
    'declared numeric maximum not enforced');
  await evaluate("enter('gap', '34', '.xp-card', '@media (min-width: 900px)')");
  await resize(800, 900);
  await check("styleA().gap === '10px' && styleB().gap === '10px'", 'media preview leaked outside its condition');
  await resize();
  await check("enter('border-left-width', '7') && styleA().borderLeftWidth === '7px' && styleB().borderLeftWidth === '2px'",
    '!important declaration priority not preserved');
  await check("enter('border-right-width', '6', '.xp-card', '@supports (display: grid)') && styleA().borderRightWidth === '6px'",
    '@supports grouping lost');

  await evaluate("enter('background-color', 'var(--lc-focus)', '.xp-card:hover')");
  await move(1, 999);
  await check("styleA().backgroundColor === styleB().backgroundColor", 'hover rule became unconditional');
  await hover("document.querySelector('#card-a')");
  await check("styleA().backgroundColor !== styleB().backgroundColor", 'hover CSS preview not applied');
  await evaluate("enter('color', 'var(--lc-focus)', '.xp-card::after')");
  await check("getComputedStyle(document.querySelector('#card-a'),'::after').color !== getComputedStyle(document.querySelector('#card-b'),'::after').color",
    'pseudo-element preview missing or leaking');
  await evaluate(`enter('color','var(--lc-focus)',':root[data-lcm-resolved-theme="light"] .xp-card')`);
  await check("styleA().color === styleB().color", 'light rule leaked into dark theme');
  for (const value of ['light', 'dark', 'light', 'dark']) {
    await theme(value);
    await check(`(styleA().color ${value === 'light' ? '!==' : '==='} styleB().color) &&
      rootIdentity === lcmXRay.root && field('padding').querySelector('input[type=number]') === keptInput`,
      'theme transition changed the component/editor identity or flattened a theme condition');
    await check(`getComputedStyle(s.querySelector('.xray-panel')).backgroundColor !== getComputedStyle(s.querySelector('.xray-panel')).color &&
      getComputedStyle(field('display').querySelector('select')).color !== getComputedStyle(field('display').querySelector('select')).backgroundColor`,
      'X-ray controls have indistinguishable foreground/background');
  }
  await action('Copy changes');
  await check(`s.querySelector('.xray-export').readOnly && s.querySelector('.xray-export').value.includes('test/integration/xray_fixture.css') &&
    s.querySelector('.xray-export').value.includes('@media') && s.querySelector('.xray-export').value.includes('!important') &&
    !s.querySelector('.xray-export').value.includes('data-lcm-preview-scope')`, 'source-aware export missing or leaking instance selectors');
  await check("JSON.stringify(originalStyles) === JSON.stringify([...document.querySelectorAll('style[data-lcm-css-source]')].map(n=>n.textContent))",
    'source stylesheet text was modified');
  await check("authoredRules.every(({rule,text})=>rule.style.cssText===text)", 'authored CSSOM declarations were modified');

  await evaluate("s.querySelector('[aria-label=\"Enable temporary CSS preview\"]').click()");
  await check("styleA().paddingTop === '12px' && lcmXRay.preview.changes().length > 0", 'preview toggle lost drafts or kept rules');
  await evaluate("s.querySelector('[aria-label=\"Enable temporary CSS preview\"]').click()");
  await check("styleA().paddingTop === '28px'", 'preview drafts did not restore');
  await evaluate("field('padding').querySelector('input[type=checkbox]').click()");
  await check("styleA().paddingTop === '12px'", 'per-property disable failed');
  await evaluate("field('padding').querySelector('input[type=checkbox]').click()");
  await evaluate("field('padding').querySelector('[aria-label=\"Reset padding\"]').click()");
  await check("styleA().paddingTop === '12px' && keptInput !== field('padding').querySelector('input[type=number]')", 'property reset failed');
  await action('Reset component'); await resetCheck();

  // Expression mode preserves authored expressions and rejects CSS injection.
  await evaluate(`const mode = field('padding').querySelector('select'); mode.value='expression'; mode.dispatchEvent(new Event('change'));`);
  await check("enter('padding','calc(1rem + 2px)') && styleA().paddingTop !== '12px'", 'expression preview failed');
  await check(`!enter('padding','1px; color: red') && !enter('padding','url(https://example.invalid/a)') &&
    !enter('padding','garbage')`, 'unsafe/invalid value accepted');
  await action('Reset all'); await resetCheck();

  // New selections preserve distinct per-instance drafts; reset-all clears both.
  await evaluate("enter('padding','20'); document.querySelector('#card-b').click(); enter('padding','30')");
  await check("styleA().paddingTop === '20px' && styleB().paddingTop === '30px'", 'instance drafts clobbered each other');
  await action('Reset component');
  await check("styleA().paddingTop === '20px' && styleB().paddingTop === '12px'", 'component reset affected sibling');
  await action('Reset all'); await resetCheck();

  // Move/resize the diagnostic window, then check it remains in a narrow viewport.
  const start = await point("s.querySelector('.xray-panel-header')");
  await move(start.x, start.y);
  await command('Input.dispatchMouseEvent', {type:'mousePressed',button:'left',buttons:1,clickCount:1,...start});
  await command('Input.dispatchMouseEvent', {type:'mouseMoved',button:'left',buttons:1,x:start.x-230,y:start.y+90});
  await command('Input.dispatchMouseEvent', {type:'mouseReleased',button:'left',buttons:0,clickCount:1,x:start.x-230,y:start.y+90});
  await check(`parseFloat(s.querySelector('.xray-panel').style.left) < ${start.x - 100}`, 'inspector drag failed');
  const handle = await point("s.querySelector('.xray-resize-handle')");
  await move(handle.x,handle.y);
  await command('Input.dispatchMouseEvent', {type:'mousePressed',button:'left',buttons:1,clickCount:1,...handle});
  await command('Input.dispatchMouseEvent', {type:'mouseMoved',button:'left',buttons:1,x:handle.x+90,y:handle.y-100});
  await command('Input.dispatchMouseEvent', {type:'mouseReleased',button:'left',buttons:0,clickCount:1,x:handle.x+90,y:handle.y-100});
  await check("!!s.querySelector('.xray-panel').style.height", 'inspector resize failed');
  await resize(640,800);
  await check("(()=>{const r=s.querySelector('.xray-panel').getBoundingClientRect(); return r.x>=0 && r.y>=0 && r.right<=640 && r.bottom<=800})()",
    'inspector escaped viewport');
  await resize();

  for (const cleanup of ["s.querySelector('.xray-close').click()", "lcmXRay.disable()",
    "document.dispatchEvent(new KeyboardEvent('keydown',{key:'Escape'}))"]) {
    await evaluate("lcmXRay.enable(); document.querySelector('#card-a').click(); enter('padding','31')");
    await evaluate(cleanup); await resetCheck();
  }
  await evaluate("lcmXRay.enable(); document.querySelector('#card-a').click(); enter('padding','31'); document.querySelector('#card-a').remove()");
  await wait("lcmXRay.preview.changes().length === 0", 'detached component preview was not released');
  await evaluate("document.querySelector('#card-b').click(); enter('padding','31'); lcmXRay.root.remove()");
  await wait("!document.querySelector('.lc-xray-host')", 'detached host controller leaked'); await resetCheck();

  await open('/fixture/readonly');
  await evaluate("document.querySelector('#card-a').click()");
  await check("!s.querySelector('.xray-css-controls') && s.querySelector('.xray-body').textContent.includes('CSS preview is disabled')",
    'host read-only policy ignored');

  // Actual workbench uses the same engine, typed catalogue and theme contract.
  await open('/workbenches/template');
  await wait("!!document.querySelector('[data-lcm-component=MenuBar]')", 'real workbench not registered');
  await evaluate("document.querySelector('[data-lcm-component=MenuBar]').click()");
  await check("[...s.querySelectorAll('.xray-css-source')].some(n=>n.textContent==='src/workbench/workbench.css')",
    'real workbench missing authored source identity');
  await check("!!s.querySelector('.xray-css-controls')", 'real workbench missing reusable CSS editors');
  await evaluate("window.originalMenuPadding = getComputedStyle(document.querySelector('.lc-wb-menubar')).padding");
  await check("enter('padding','0px 32px','.lc-wb-menubar') && getComputedStyle(document.querySelector('.lc-wb-menubar')).paddingRight === '32px'",
    'real workbench compound-value preview failed');
  for (const value of ['dark', 'light']) {
    await theme(value);
    await check(`(()=>{const field=s.querySelector('.xray-css-controls');const r=field.getBoundingClientRect();
      const body=s.querySelector('.xray-body'); return field.scrollWidth<=field.clientWidth+1 && body.scrollWidth<=body.clientWidth+1;})()`,
      'CSS editor horizontally overflows in ' + value);
  }
  await action('Reset component');
  await check("getComputedStyle(document.querySelector('.lc-wb-menubar')).padding === originalMenuPadding", 'real workbench reset failed');
  await resetCheck();
  // Even a root edit that hides the workbench must leave recovery reachable.
  await evaluate("document.querySelector('[data-lcm-component=Workbench]').click()");
  await check("enter('display','none','.lc-wb-shell') && getComputedStyle(lcmXRay.root).display === 'none'",
    'hostile root preview was not exercised');
  await check("!s.querySelector('.xray-panel').hidden && !action('Reset all').disabled", 'root preview hid its own recovery surface');
  await click("action('Reset all')");
  await check("getComputedStyle(lcmXRay.root).display === 'grid' && rootIdentity === lcmXRay.root", 'full reset failed to restore the original workbench without remounting');
  await resetCheck();
  // A nested ownership tree: scoped resets include every depth, not siblings.
  await open('/fixture/tree');
  await evaluate(`window.pick = id => document.getElementById(id).click();
    window.children = s.querySelector('[aria-label="Apply to children"]');
    window.padding = id => getComputedStyle(document.getElementById(id)).paddingTop;
    window.seedTree = () => {
      pick('card-a'); enter('padding','20');
      pick('card-b'); enter('padding','30');
      pick('card-c'); enter('padding','40');
      pick('group-inner'); enter('margin','10','.xp-group');
      pick('group-a'); enter('padding','24','.xp-group');
    };
    seedTree();
  `);
  await check(`!children.checked && lcmXRay.preview.changes().length === 5 &&
    s.querySelector('.xray-preview-status').textContent.includes('1 in selection') &&
    s.querySelector('.xray-preview-status').textContent.includes('3 in children')`, 'reset scope counts must include recursive children');
  await action('Reset component');
  await check(`padding('group-a') === '8px' && padding('card-a') === '20px' && padding('card-b') === '30px' &&
    padding('card-c') === '40px' && lcmXRay.preview.changes().length === 4 && action('Reset component').disabled`,
    'component-only reset changed child drafts or sibling branches');
  await evaluate("children.click()");
  await check("!action('Reset component').disabled && s.querySelector('.xray-preview-status').textContent.includes('3 in selection + children')",
    'checking children must expose nested-only recovery');
  await evaluate("document.querySelector('#group-inner').hidden = true");
  await action('Reset component');
  await evaluate("document.querySelector('#group-inner').hidden = false");
  await check(`padding('card-a') === '12px' && padding('card-b') === '12px' &&
    getComputedStyle(document.querySelector('#group-inner')).marginTop === '0px' &&
    padding('card-c') === '40px' && lcmXRay.preview.changes().length === 1`, 'recursive reset missed hidden/deep descendants or reset a sibling');

  await evaluate("seedTree(); pick('group-inner')");
  await action('Reset component');
  await check(`padding('card-b') === '12px' && padding('card-a') === '20px' && padding('card-c') === '40px' &&
    padding('group-a') === '24px'`, 'nested group reset crossed its ownership boundary');
  await evaluate("document.querySelector('[data-lcm-component=Workbench]').click()");
  await action('Reset component'); await resetCheck();
  await check("padding('group-a') === '8px' && padding('card-a') === '12px' && padding('card-c') === '12px'", 'workbench recursive recovery failed');

  // Changed highlights explain original → proposal, including invalid and paused edits.
  for (const value of ['dark', 'light']) {
    await theme(value);
    await evaluate("pick('card-a'); enter('padding','22')");
    await check(`field('padding').classList.contains('is-changed') && !field('padding').querySelector('.xray-css-diff').hidden &&
      field('padding').querySelector('.xray-original-value').textContent === '12px' &&
      field('padding').querySelector('.xray-proposed-value').textContent === '22px' &&
      getComputedStyle(field('padding')).borderInlineStartColor !== 'rgba(0, 0, 0, 0)'`, 'changed-value comparison missing in ' + value);
    await evaluate("pick('card-c'); pick('card-a')");
    await check("field('padding').classList.contains('is-changed') && field('padding').querySelector('.xray-proposed-value').textContent === '22px'",
      'changed highlight was lost when reselecting');
    await evaluate("field('padding').querySelector('input[type=checkbox]').click()");
    await check("padding('card-a') === '12px' && field('padding').textContent.includes('Override (disabled)') && field('padding').classList.contains('is-changed')",
      'disabled changes must remain visible and resettable');
    await evaluate("field('padding').querySelector('input[type=checkbox]').click(); enter('padding','-9')");
    await check(`padding('card-a') === '22px' && field('padding').classList.contains('is-invalid') &&
      field('padding').querySelector('.xray-proposed-value').textContent === '-9px' &&
      field('padding').textContent.includes('last valid override: 22px') && action('Copy changes').disabled`,
      'invalid draft must retain the last valid preview and explain the difference');
    await evaluate("pick('card-c'); pick('card-a')");
    await check(`field('padding').classList.contains('is-invalid') &&
      field('padding').querySelector('input[type=number]').value === '-9' && !action('Reset all').disabled`,
      'invalid draft or its reset was lost on selection change');
    await evaluate("field('padding').querySelector('[aria-label=\"Reset padding\"]').click()");
    await check(`padding('card-a') === '12px' && !field('padding').classList.contains('is-changed') &&
      field('padding').querySelector('.xray-css-diff').hidden`, 'per-property reset did not clear its comparison');
    await evaluate("enter('padding','-9')");
    await check("!action('Reset all').disabled && s.querySelector('.xray-preview-status').textContent.includes('1 invalid')", 'an invalid-only draft must not disable full reset');
    await action('Reset all'); await resetCheck();
    await check("field('padding').querySelector('input[type=number]').value === '12'", 'full reset did not discard invalid input');
  }
  await evaluate(`seedTree(); pick('card-b'); children.checked = false; children.dispatchEvent(new Event('change'));
    s.querySelector('[aria-label="Enable temporary CSS preview"]').click()`);
  await action('Reset all');
  await evaluate("s.querySelector('[aria-label=\"Enable temporary CSS preview\"]').click()");
  await resetCheck();
  await check("padding('card-a') === '12px' && padding('card-b') === '12px' && padding('card-c') === '12px'",
    'Reset all must include all branches even when children is unchecked and previews were paused');
  assert.deepEqual(errors, [], 'uncaught browser exceptions');
  console.log('X-ray preview passed: selection, isolated CSS, conditional rules, themes, bindings, subtree/global recovery, dirty-value comparisons, export and cleanup.');
} finally { socket.close(); }
