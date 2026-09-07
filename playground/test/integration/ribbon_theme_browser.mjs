#!/usr/bin/env node
import assert from 'node:assert/strict';
import { checkInteractionStates } from './interaction_state_contract.mjs';

const base = process.argv[2] ?? 'http://127.0.0.1:18102';
const debug = process.argv[3] ?? 'http://127.0.0.1:19342';
const pages = await (await fetch(debug + '/json')).json();
const target = pages.find(page => page.type === 'page');
const socket = new WebSocket(target.webSocketDebuggerUrl);
await new Promise(resolve => socket.addEventListener('open', resolve, {once: true}));
let sequence = 0;
const pending = new Map();
socket.addEventListener('message', event => {
  const message = JSON.parse(event.data);
  const request = pending.get(message.id);
  if (!request) return;
  pending.delete(message.id);
  message.error ? request.reject(Error(message.error.message)) : request.resolve(message.result);
});
const command = (method, params = {}, sessionId) => new Promise((resolve, reject) => {
  const id = ++sequence;
  pending.set(id, {resolve, reject});
  socket.send(JSON.stringify({id, method, params, sessionId}));
});
const evaluate = async (expression, sessionId) => {
  const result = await command('Runtime.evaluate', {expression, returnByValue: true, awaitPromise: true}, sessionId);
  if (result.exceptionDetails) throw Error(result.exceptionDetails.exception?.description);
  return result.result.value;
};
const wait = async (expression, message, sessionId) => {
  const deadline = Date.now() + 15000;
  while (Date.now() < deadline) {
    if (await evaluate(expression, sessionId)) return;
    await new Promise(resolve => setTimeout(resolve, 75));
  }
  throw Error(message);
};

let themeTarget;
try {
  // Change the actual publisher preference in a second tab. Real storage
  // events update already-mounted widget documents, iframe children and shells.
  themeTarget = (await command('Target.createTarget', {url: base + '/'})).targetId;
  const themeSession = (await command('Target.attachToTarget', {targetId: themeTarget, flatten: true})).sessionId;
  await wait('Boolean(window.LineCableModelsTheme)', 'publisher theme did not initialize', themeSession);
  await command('Target.activateTarget', {targetId: target.id});
  await command('Page.enable');
  await command('Network.enable');
  await command('Network.setCacheDisabled', {cacheDisabled: true});
  await command('Emulation.setDeviceMetricsOverride', {
    width: 1600, height: 1000, deviceScaleFactor: 1, mobile: false,
  });
  const report = [];
  for (const route of ['/widgets/ribbon', '/fixture/embedded', '/fixture/standalone', '/fixture/workbench']) {
    console.log('Checking ribbon theme transitions: ' + route);
    await command('Page.navigate', {url: base + route});
    const embedded = route === '/fixture/embedded';
    const documentExpression = embedded
      ? "document.querySelector('#ribbon-frame')?.contentDocument" : 'document';
    await wait(`location.pathname === ${JSON.stringify(route)} &&
      ${documentExpression}?.querySelector('.lc-ribbon')`, 'ribbon did not mount: ' + route);
    await evaluate(`window.auditDocument = ${documentExpression};
      window.auditRibbon = auditDocument.querySelector('.lc-ribbon');
      window.auditCounter = auditDocument.querySelector('#fixture-events');`);
    // Bonito may render markup before its session-local handlers arrive.
    await wait(`(() => {
      auditRibbon.querySelector('[data-ribbon-collapse]').click();
      return auditRibbon.dataset.collapsed === 'true';
    })()`, 'ribbon behavior did not initialize: ' + route);
    await evaluate("auditRibbon.querySelector('[data-ribbon-collapse]').click()");

    const inspect = async () => evaluate(`(() => {
      const doc = auditDocument, win = doc.defaultView;
      const token = (node, name) => {
        const probe = doc.createElement('span');
        probe.style.color = 'var(' + name + ')'; node.append(probe);
        const value = win.getComputedStyle(probe).color; probe.remove(); return value;
      };
      const buttons = [...doc.querySelectorAll('.lc-toolbar-button')];
      return buttons.map(button => {
        const s = win.getComputedStyle(button);
        const active = button.classList.contains('lc-toolbar-button-active');
        const hovered = button.matches(':hover') && !button.disabled;
        return {
          action: button.dataset.toolbarAction, active, disabled: button.disabled,
          color: s.color, bg: s.backgroundColor,
          icon: win.getComputedStyle(button.querySelector('svg')).stroke,
          expectedColor: token(button, active ? '--lc-accent-ink' : hovered ? '--lc-strong-text' : '--lc-widget-text'),
          expectedBg: token(button, active ? '--lc-widget-focus' : hovered ? '--lc-hover-bg-strong' : '--lc-toolbar-bg'),
          opacity: s.opacity, busy: button.getAttribute('aria-busy'),
          spinner: win.getComputedStyle(button, '::after').borderTopColor,
        };
      });
    })()`);
    const check = async label => {
      const buttons = await inspect();
      assert(buttons.some(button => button.active), label + ': no active button');
      for (const button of buttons) {
        assert.equal(button.color, button.expectedColor, label + ': text ' + button.action);
        assert.equal(button.bg, button.expectedBg, label + ': background ' + button.action);
        assert.equal(button.icon, button.color, label + ': icon ' + button.action);
        if (button.disabled) assert(Number(button.opacity) < 1, label + ': disabled cue');
        if (button.busy === 'true') assert.equal(button.spinner, button.color, label + ': busy spinner');
      }
      return buttons.length;
    };

    for (const theme of ['light', 'dark', 'light', 'dark']) {
      await evaluate(`LineCableModelsTheme.select(${JSON.stringify(theme)})`, themeSession);
      await wait(`auditDocument.documentElement.dataset.lcmResolvedTheme === ${JSON.stringify(theme)}`,
        route + ': child theme did not follow ' + theme);
      await command('Input.dispatchMouseEvent', {type: 'mouseMoved', x: 1599, y: 999});
      const count = await check(route + ' ' + theme);
      assert(await evaluate(`auditRibbon === auditDocument.querySelector('.lc-ribbon') &&
        auditCounter === auditDocument.querySelector('#fixture-events')`), 'theme change remounted controls');
      // Hover active, ordinary, disabled and busy controls in the visible panel.
      for (const state of ['.lc-toolbar-button-active', '.lc-toolbar-button:not(.lc-toolbar-button-active):not(:disabled)',
        '.lc-toolbar-button:disabled', '.lc-toolbar-button-busy']) {
        const point = await evaluate(`(() => {
          const node = auditRibbon.querySelector('.lc-ribbon-panel:not([hidden]) ' + ${JSON.stringify(state)});
          node.scrollIntoView({block:'nearest'});
          const r = node.getBoundingClientRect();
          const f = document.querySelector('#ribbon-frame')?.getBoundingClientRect();
          return {x:r.x + r.width/2 + (f?.x || 0), y:r.y + r.height/2 + (f?.y || 0)};
        })()`);
        await command('Input.dispatchMouseEvent', {type: 'mouseMoved', ...point});
        await check(route + ' ' + theme + ' hover ' + state);
      }
      report.push({route, theme, buttons: count});
    }
    // Moving groups into overflow or changing tabs cannot change their colors.
    await command('Emulation.setDeviceMetricsOverride', {
      width: 640, height: 900, deviceScaleFactor: 1, mobile: false,
    });
    if (embedded) await evaluate("document.querySelector('#ribbon-frame').width = '580'");
    await evaluate(`new Promise(resolve => setTimeout(resolve, 300))`);
    const tabs = await evaluate(`[...auditRibbon.querySelectorAll('[data-ribbon-tab]:not(:disabled)')]
      .map(tab => tab.dataset.ribbonTab)`);
    let overflowed = false;
    for (const tab of tabs) {
      await evaluate(`auditRibbon.querySelector('[data-ribbon-tab="' + ${JSON.stringify(tab)} + '"]').click()`);
      await evaluate(`new Promise(resolve => setTimeout(resolve, 150))`);
      overflowed ||= await evaluate(`Boolean(auditRibbon.querySelector('.lc-ribbon-panel:not([hidden]) .lc-ribbon-overflow-menu .lc-ribbon-group'))`);
      await check(route + ' narrow tab ' + tab);
    }
    if (route === '/widgets/ribbon' || embedded) assert(overflowed, 'overflow case was not exercised');
    // Restore a hovered, active control after the round trip and verify callbacks.
    if (!route.startsWith('/widgets/') && !embedded) {
      await evaluate(`auditDocument.querySelector('.lc-toolbar [data-toolbar-action="bus"]').click()`);
      await wait(`auditCounter.textContent === '1'`, route + ': callback failed after theme switching');
    }
    await command('Emulation.setDeviceMetricsOverride', {
      width: 1600, height: 1000, deviceScaleFactor: 1, mobile: false,
    });
  }
  const interactionStates = await checkInteractionStates({base, command, evaluate, wait,
    selectTheme: theme => evaluate(`LineCableModelsTheme.select(${JSON.stringify(theme)})`, themeSession)});
  console.log(JSON.stringify({ribbonThemeTransitions: report, interactionStates}, null, 2));
} finally {
  if (themeTarget) await command('Target.closeTarget', {targetId: themeTarget});
  socket.close();
}
