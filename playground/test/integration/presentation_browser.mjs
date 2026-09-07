#!/usr/bin/env node

import { writeFile } from "node:fs/promises";
import { assertPublishedText } from "./published_text_contract.mjs";
import { assertMathNotes } from "./math_notes_browser.mjs";
import { assertPublishedShell } from "./published_shell_browser.mjs";
import { assertIncrementalLists } from "./incremental_lists_browser.mjs";

const baseUrl = process.argv[2] ?? "http://127.0.0.1:8080";
const debuggingUrl = process.argv[3] ?? "http://127.0.0.1:9222";
const deckUrl = `${baseUrl}/presentations/specimen.html`;

class DevTools {
  constructor(url) {
    this.socket = new WebSocket(url);
    this.sequence = 0;
    this.pending = new Map();
    this.events = new Map();
  }

  async connect() {
    await new Promise((resolve, reject) => {
      this.socket.addEventListener("open", resolve, { once: true });
      this.socket.addEventListener("error", reject, { once: true });
    });
    this.socket.addEventListener("message", event => {
      const message = JSON.parse(event.data);
      if (message.id) {
        const pending = this.pending.get(message.id);
        if (!pending) return;
        this.pending.delete(message.id);
        if (message.error) pending.reject(new Error(message.error.message));
        else pending.resolve(message.result);
        return;
      }
      const listeners = this.events.get(message.method) ?? [];
      listeners.splice(0).forEach(resolve => resolve(message.params));
    });
  }

  command(method, params = {}) {
    const id = ++this.sequence;
    const result = new Promise((resolve, reject) => this.pending.set(id, { resolve, reject }));
    this.socket.send(JSON.stringify({ id, method, params }));
    return result;
  }

  once(method) {
    return new Promise(resolve => {
      const listeners = this.events.get(method) ?? [];
      listeners.push(resolve);
      this.events.set(method, listeners);
    });
  }

  close() {
    this.socket.close();
  }
}

function assert(condition, message) {
  if (!condition) throw new Error(message);
}

async function evaluate(devtools, expression) {
  const result = await devtools.command("Runtime.evaluate", {
    expression,
    awaitPromise: true,
    returnByValue: true,
  });
  if (result.exceptionDetails) {
    throw new Error(result.exceptionDetails.exception?.description ?? result.exceptionDetails.text);
  }
  return result.result.value;
}

async function navigate(devtools, url) {
  const loaded = devtools.once("Page.loadEventFired");
  await devtools.command("Page.navigate", { url });
  let timer;
  try {
    await Promise.race([
      loaded,
      new Promise((_, reject) => { timer = setTimeout(
        () => reject(new Error(`page load timed out: ${url}`)), 15_000,
      ); }),
    ]);
  } finally { clearTimeout(timer); }
}

async function waitUntil(devtools, expression, message, timeout = 15_000) {
  const deadline = Date.now() + timeout;
  while (Date.now() < deadline) {
    if (await evaluate(devtools, expression)) return;
    await new Promise(resolve => setTimeout(resolve, 75));
  }
  throw new Error(message);
}

async function setViewport(devtools, width, height) {
  await devtools.command("Emulation.setDeviceMetricsOverride", {
    width,
    height,
    deviceScaleFactor: 1,
    mobile: false,
  });
  await evaluate(devtools, 'new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)))');
  await waitUntil(devtools,
    `!document.documentElement.classList.contains('lcm-viewport-settling')`,
    `viewport did not settle at ${width} × ${height}`);
}

async function inspectGeometry(devtools) {
  return evaluate(devtools, `(() => {
    const stage = document.querySelector('.reveal .slides');
    const box = stage.getBoundingClientRect();
    const liveFrames = Array.from(document.querySelectorAll('.lcm-live-viewport iframe'));
    const scaled = [];
    for (const frame of liveFrames) {
      for (let node = frame.parentElement; node && node !== document.documentElement; node = node.parentElement) {
        const transform = getComputedStyle(node).transform;
        if (!transform || transform === 'none') continue;
        const matrix = new DOMMatrix(transform);
        const scaleX = Math.hypot(matrix.a, matrix.b);
        const scaleY = Math.hypot(matrix.c, matrix.d);
        if (Math.abs(scaleX - 1) > 0.001 || Math.abs(scaleY - 1) > 0.001) {
          scaled.push({className: node.className, transform, scaleX, scaleY});
        }
      }
    }
    const overflows = Array.from(document.querySelectorAll(
      '.slides > section.present, section.present .lcm-slot, section.present .lcm-layout'
    )).map(slide => ({
      id: slide.id || slide.className,
      x: slide.scrollWidth - slide.clientWidth,
      y: slide.scrollHeight - slide.clientHeight,
    })).filter(sample => sample.x > 1 || sample.y > 1);
    return {
      viewport: {width: innerWidth, height: innerHeight},
      stage: {left: box.left, top: box.top, width: box.width, height: box.height},
      scaled,
      overflows,
    };
  })()`);
}

async function assertGeometry(devtools, width, height) {
  await setViewport(devtools, width, height);
  const sample = await inspectGeometry(devtools);
  assert(Math.abs(sample.stage.width / sample.stage.height - 16 / 9) < 0.003,
    `stage ratio drifted at ${width} × ${height}`);
  assert(sample.stage.left >= -1 && sample.stage.top >= -1,
    `stage escaped viewport at ${width} × ${height}`);
  assert(sample.stage.left + sample.stage.width <= width + 1,
    `stage exceeds horizontal viewport at ${width} × ${height}`);
  assert(sample.stage.top + sample.stage.height <= height + 1,
    `stage exceeds vertical viewport at ${width} × ${height}`);
  assert(sample.scaled.length === 0,
    `live viewport has a scaled ancestor: ${JSON.stringify(sample.scaled)}`);
  assert(sample.overflows.length === 0,
    `slide content overflows: ${JSON.stringify(sample.overflows)}`);
  return sample;
}

async function click(devtools, selector) {
  await waitUntil(devtools, `document.getAnimations().every(a => a.playState !== 'running')`,
    'click target is still moving');
  const point = await evaluate(devtools, `(() => {
    const node = document.querySelector(${JSON.stringify(selector)});
    if (!node) throw new Error('Missing click target');
    const r = node.getBoundingClientRect();
    return { x: r.x + r.width / 2, y: r.y + r.height / 2 };
  })()`);
  await devtools.command("Input.dispatchMouseEvent", {
    type: "mousePressed", ...point, button: "left", buttons: 1, clickCount: 1,
  });
  await devtools.command("Input.dispatchMouseEvent", {
    type: "mouseReleased", ...point, button: "left", buttons: 0, clickCount: 1,
  });
}

async function menuAction(devtools, handler) {
  await click(devtools, '.lcm-deck-status button[data-action="menu"]');
  await waitUntil(devtools, 'Reveal.getPlugin("menu").isOpen()', 'menu did not open');
  // The vendor menu animates for 300 ms. Wait for its actual position to settle.
  await waitUntil(devtools, `Math.abs(document.querySelector('.slide-menu.active')
    .getBoundingClientRect().left) < 1`, 'menu did not settle');
  await click(devtools, '.slide-menu-toolbar [data-panel="Custom0"]');
  await click(devtools, '[onclick="RevealMenuToolHandlers.' + handler + '(event)"]');
}

const pages = await fetch(`${debuggingUrl}/json`);
assert(pages.ok, `cannot reach Chrome DevTools at ${debuggingUrl}`);
const targets = await pages.json();
const target = targets.find(candidate => candidate.type === "page");
assert(target, "no Chrome page target is available");

const devtools = new DevTools(target.webSocketDebuggerUrl);
await devtools.connect();
await devtools.command("Page.enable");
await devtools.command("Runtime.enable");
await devtools.command("Network.enable");
await devtools.command("Network.setCacheDisabled", { cacheDisabled: true });
await devtools.command("Emulation.setEmulatedMedia", { media: "" });

try {
  await assertPublishedShell({devtools, baseUrl, evaluate, navigate, waitUntil, setViewport, assert});
  await assertIncrementalLists({devtools, baseUrl, evaluate, navigate, waitUntil, assert});
  await setViewport(devtools, 1920, 1080);
  await navigate(devtools, `${baseUrl}/presentations/`);
  const layoutCards = await evaluate(devtools, `[...document.querySelectorAll('.lc-gallery-card')]
    .filter(card => card.querySelector('a[href*="layouts.html#"]')).map(card => ({
      preview: card.querySelector('a[href*="specimen.html#/"]')?.href,
      source: card.querySelector('a[href*="layouts.html#"]').href,
      layout: card.querySelector('a[href*="layouts.html#"]').hash.slice(1)
    }))`);
  assert(layoutCards.length === 6 && layoutCards.every(card => card.preview && card.source),
    "a reusable layout card lacks a preview or copyable source");
  assert(await evaluate(devtools, `document.querySelector('#complete-decks') &&
    document.querySelector('#reusable-slide-layouts') &&
    getComputedStyle(document.querySelector('main')).display === 'block'`),
    "gallery does not distinguish complete decks from reusable slide layouts");
  const mathGuide = await evaluate(devtools,
    `document.querySelector('a[href$="math-notes.html"]')?.href`);
  assert(mathGuide, 'presentation gallery does not expose the equation authoring guide');
  await navigate(devtools, mathGuide);
  assert(await evaluate(devtools, `document.querySelector('code.sourceCode')?.textContent
    .includes('cssId{impedance-imag}') && document.querySelector('.code-copy-button') &&
    !document.querySelector('iframe, .lcm-math-callout') && !window.LCMMathNotes`),
    'equation guide is not a copyable, inert published page');
  for (const card of layoutCards) {
    await navigate(devtools, card.source);
    assert(await evaluate(devtools, `(() => {
      const section = document.getElementById(${JSON.stringify(card.layout)});
      return section?.querySelector('code.sourceCode')?.textContent.includes(
        ${JSON.stringify("{.lcm-layout-")} + ${JSON.stringify(card.layout)} + '}') &&
        Boolean(section.querySelector('.code-copy-button')) &&
        document.querySelectorAll('iframe').length === 0;
    })()`), `${card.layout} source link does not offer an inert copyable slide`);
    await navigate(devtools, card.preview);
    await waitUntil(devtools, `document.documentElement.dataset.lcmDeckReady === 'true'`,
      `${card.layout} preview did not initialize`);
    assert(await evaluate(devtools, `Reveal.getCurrentSlide().querySelector('[data-lcm-layout]')
      ?.dataset.lcmLayout === ${JSON.stringify(card.layout)}`),
      `${card.layout} preview link opens the wrong layout`);
  }
  await navigate(devtools, deckUrl);
  await waitUntil(devtools,
    `document.querySelector('.reveal.ready') && document.documentElement.dataset.lcmDeckMode === 'audience'`,
    "audience deck did not initialize");
  await evaluate(devtools, `localStorage.setItem('lcm.playground.theme', 'light')`);
  await navigate(devtools, deckUrl);
  await waitUntil(devtools, `document.documentElement.dataset.lcmDeckReady === 'true'`,
    "cached-theme deck did not initialize");
  assert(await evaluate(devtools, `document.documentElement.dataset.lcmResolvedTheme === 'light'`),
    "deck ignored the playground's cached light theme");
  await evaluate(devtools, `document.fonts.ready.then(() => new Promise(resolve =>
    requestAnimationFrame(() => requestAnimationFrame(resolve))))`);
  for (const theme of ["dark", "light"]) {
    await evaluate(devtools, `LineCableModelsTheme.select(${JSON.stringify(theme)})`);
    await assertPublishedText(devtools, '#title-slide h1');
    assert(await evaluate(devtools, `(() => {
      const pointer = document.querySelector('.lcm-pointer');
      return getComputedStyle(pointer, '::before').borderTopColor === 'rgb(255, 0, 0)' &&
        getComputedStyle(pointer, '::after').backgroundColor === 'rgb(255, 0, 0)';
    })()`), `presentation laser is not red in ${theme} theme`);
    assert(await evaluate(devtools, `(() => {
      const home = document.querySelector('.lcm-deck-status a[data-action="home"]');
      const menu = document.querySelector('.lcm-deck-status button[data-action="menu"]');
      return home?.getAttribute('aria-label') === 'Return to playground home' &&
        home.title === home.getAttribute('aria-label') &&
        getComputedStyle(home).color === getComputedStyle(menu).color &&
        getComputedStyle(home.querySelector('svg')).stroke === getComputedStyle(home).color &&
        getComputedStyle(home).backgroundColor === getComputedStyle(menu).backgroundColor;
    })()`), `playground home icon lost its accessible label or ${theme} styling`);
    await click(devtools, '.lcm-deck-status button[data-action="menu"]');
    await waitUntil(devtools, `Reveal.getPlugin('menu').isOpen() &&
      Math.abs(document.querySelector('.slide-menu.active').getBoundingClientRect().left) < 1`,
      'presentation tools menu did not settle');
    await click(devtools, '.slide-menu-toolbar [data-panel="Custom0"]');
    assert(await evaluate(devtools, `(() => {
      const tools = document.querySelector('.slide-menu-wrapper');
      return !tools.querySelector('a[data-action="home"], .lcm-deck-home, [onclick*="toggleScrollView"]') &&
        !tools.textContent.includes('Playground home');
    })()`), `Tools contains playground navigation or unsupported scroll mode in ${theme}`);
    await evaluate(devtools, "Reveal.getPlugin('menu').closeMenu()");
    await waitUntil(devtools, "!Reveal.getPlugin('menu').isOpen()", 'presentation menu did not close');
    // The vendor reports closed before its overlay finishes fading. Wait until
    // that overlay stops intercepting the next real text-selection gesture.
    await waitUntil(devtools, `[...document.querySelectorAll('.slide-menu, .slide-menu-overlay')]
      .every(node => node.getAnimations().every(a => a.playState !== 'running'))`,
      'presentation menu closing animation did not finish');
  }
  await click(devtools, '.lcm-deck-status a[data-action="home"]');
  await waitUntil(devtools, `location.pathname === '/' &&
    document.querySelector('#quarto-document-content h1')?.textContent.includes('LineCableModels playground')`,
    'footer home icon did not navigate to the playground landing page');
  await navigate(devtools, deckUrl);
  await waitUntil(devtools, `document.documentElement.dataset.lcmDeckReady === 'true'`,
    'deck did not initialize after verifying playground navigation');
  assert(await evaluate(devtools,
    `Array.from(document.querySelectorAll('iframe[data-lcm-src]')).every(frame => !frame.hasAttribute('src'))`),
    "off-screen live frames were activated eagerly");

  const geometry = [];
  for (const [width, height] of [[1280, 720], [1920, 1080], [1280, 1024], [1366, 768], [2560, 1080]]) {
    for (const index of [0, 1, 2, 4, 5, 6, 8]) {
      await evaluate(devtools, `Reveal.slide(${index})`);
      geometry.push(await assertGeometry(devtools, width, height));
    }
  }

  await assertMathNotes({ devtools, evaluate, waitUntil, click, setViewport, assert });

  const focusResult = await evaluate(devtools, `(() => {
    const range = document.querySelector('#specimen-range');
    const slide = range.closest('.slides > section');
    const indices = Reveal.getIndices(slide);
    Reveal.slide(indices.h, indices.v || 0);
    range.focus();
    return {before: Reveal.getIndices(), value: Number(range.value)};
  })()`);
  await devtools.command("Input.dispatchKeyEvent", {
    type: "keyDown", key: "ArrowRight", code: "ArrowRight", windowsVirtualKeyCode: 39,
  });
  await devtools.command("Input.dispatchKeyEvent", {
    type: "keyUp", key: "ArrowRight", code: "ArrowRight", windowsVirtualKeyCode: 39,
  });
  const focusAfter = await evaluate(devtools, `({
    index: Reveal.getIndices(),
    value: Number(document.querySelector('#specimen-range').value),
    focused: document.activeElement === document.querySelector('#specimen-range')
  })`);
  assert(focusAfter.index.h === focusResult.before.h && focusAfter.index.v === focusResult.before.v,
    "Reveal stole an arrow key from a focused range control");
  assert(focusAfter.focused, "range control lost focus after its arrow key");
  assert(focusAfter.value > focusResult.value, "range control did not consume its arrow key");

  await evaluate(devtools, `document.activeElement?.blur()`);
  await devtools.command("Input.dispatchKeyEvent", {
    type: "keyDown", key: "l", code: "KeyL", windowsVirtualKeyCode: 76,
  });
  await devtools.command("Input.dispatchKeyEvent", {
    type: "keyUp", key: "l", code: "KeyL", windowsVirtualKeyCode: 76,
  });
  assert(await evaluate(devtools, `document.documentElement.classList.contains('lcm-pointer-enabled')`),
    "presentation pointer did not toggle");
  assert(await evaluate(devtools, `Reveal.getIndices().h === ${focusAfter.index.h}`),
    "L toggled the laser but also advanced the slide");
  assert(await evaluate(devtools, `document.querySelector('.lcm-deck-status button[aria-pressed="true"]')
    ?.textContent === 'Laser on · L'`), "laser state is not visible");

  await devtools.command("Fetch.enable", { patterns: [{ urlPattern: "*/presentations/probe" }] });
  const pendingLiveRequest = devtools.once("Fetch.requestPaused");
  const persistentBefore = await evaluate(devtools, `(() => {
    const frame = document.querySelector('iframe[data-lcm-src="/presentations/probe"]');
    const slide = frame.closest('.slides > section');
    const indices = Reveal.getIndices(slide);
    Reveal.slide(indices.h, indices.v || 0);
    return {
      indices,
      src: frame.getAttribute('src')
    };
  })()`);
  const pausedLive = await pendingLiveRequest;
  assert(await evaluate(devtools, `document.querySelector('.lcm-deck-status output')
      .textContent === 'Loading live view…' &&
    document.documentElement.classList.contains('lcm-pointer-enabled') &&
    !document.querySelector('.lcm-deck-status button[aria-pressed]').disabled`),
    "live loading is not reported, or blocks the already usable laser");
  await devtools.command("Fetch.continueRequest", { requestId: pausedLive.requestId });
  await devtools.command("Fetch.disable");
  await waitUntil(devtools,
    `document.querySelector('iframe[data-lcm-src="/presentations/probe"]')?.contentDocument?.querySelector('.lc-presentation-probe')`,
    "Bonito canvas probe did not mount on first entry");
  await waitUntil(devtools,
    `document.querySelector('.lcm-deck-status output').textContent === 'Ready'`,
    "status did not report readiness after the child mounted");
  await evaluate(devtools, `document.querySelector('iframe[data-lcm-src="/presentations/probe"]')
    .contentWindow.__lcmPersistenceToken = 'retained'`);
  persistentBefore.session = await evaluate(devtools,
    `document.querySelector('iframe[data-lcm-src="/presentations/probe"]')
      .contentDocument.querySelector('.lc-presentation-probe').dataset.sessionId`);
  persistentBefore.src = await evaluate(devtools,
    `document.querySelector('iframe[data-lcm-src="/presentations/probe"]').getAttribute('src')`);
  await evaluate(devtools, `Reveal.next()`);
  await evaluate(devtools, `Reveal.slide(${persistentBefore.indices.h}, ${persistentBefore.indices.v ?? 0})`);
  await waitUntil(devtools,
    `document.querySelector('iframe[data-lcm-src="/presentations/probe"]')?.contentWindow?.__lcmPersistenceToken === 'retained'`,
    "live frame identity was replaced during navigation");
  await waitUntil(devtools,
    `(() => {
      const root = document.querySelector('iframe[data-lcm-src="/presentations/probe"]')?.contentDocument?.documentElement;
      return Number(root?.dataset.lcmDeckSlideEnterCount || 0) >= 1 && root.classList.contains('lcm-deck-visible');
    })()`,
    "child did not settle into its entered lifecycle state");
  const persistentAfter = await evaluate(devtools, `(() => {
    const frame = document.querySelector('iframe[data-lcm-src="/presentations/probe"]');
    const child = frame.contentDocument;
    return {
      session: child.querySelector('.lc-presentation-probe').dataset.sessionId,
      src: frame.getAttribute('src'),
      event: child.documentElement.dataset.lcmDeckEvent,
      enterCount: Number(child.documentElement.dataset.lcmDeckSlideEnterCount || 0),
      visible: child.documentElement.classList.contains('lcm-deck-visible')
    };
  })()`);
  assert(persistentAfter.session === persistentBefore.session, "Bonito session identity changed");
  assert(persistentAfter.src === persistentBefore.src, "live frame source changed during navigation");
  assert(persistentAfter.enterCount >= 1 && persistentAfter.visible,
    `child did not receive slide lifecycle messages: ${JSON.stringify(persistentAfter)}`);

  await menuAction(devtools, 'overview');
  await waitUntil(devtools,
    `document.querySelector('.reveal')?.classList.contains('overview') &&
      document.documentElement.classList.contains('lcm-overview-active')`,
    "overview mode did not activate");
  const overviewResult = await evaluate(devtools, `(() => {
    const slides = Array.from(document.querySelectorAll('.slides > section'));
    const frames = Array.from(document.querySelectorAll('.lcm-live-viewport iframe'));
    const placeholders = Array.from(document.querySelectorAll('.lcm-live-placeholder'));
    return {
      slides: slides.length,
      visibleSlides: slides.filter(slide => getComputedStyle(slide).display !== 'none').length,
      hiddenFrames: frames.filter(frame => getComputedStyle(frame).visibility === 'hidden').length,
      visiblePlaceholders: placeholders.filter(node => getComputedStyle(node).display !== 'none').length,
      warning: Boolean(document.querySelector('.lcm-deck-warning'))
    };
  })()`);
  assert(overviewResult.visibleSlides === overviewResult.slides,
    `overview omitted slides: ${JSON.stringify(overviewResult)}`);
  assert(overviewResult.hiddenFrames === overviewResult.visiblePlaceholders,
    `overview did not replace every live frame: ${JSON.stringify(overviewResult)}`);
  assert(!overviewResult.warning, "overview raised a false transformed-ancestor warning");
  assert(await evaluate(devtools, `(() => {
    const r = Reveal.getCurrentSlide().getBoundingClientRect();
    const s = getComputedStyle(document.querySelector('.reveal .slides'));
    return r.width > 100 && r.height > 50 && r.left >= 0 && r.top >= 0 &&
      r.right <= innerWidth && r.bottom <= innerHeight - 34 && s.overflow === 'visible';
  })()`), "current overview thumbnail is clipped or the neighbors are hidden");
  // Select a visible neighboring thumbnail with a real mouse click, then return
  // to the probe to check that neither operation remounts its session.
  await click(devtools, '.slides > section:nth-child(7)');
  await waitUntil(devtools, '!Reveal.isOverview() && Reveal.getIndices().h === 6',
    'clicking an overview thumbnail did not navigate and close overview');
  await evaluate(devtools, `Reveal.slide(${persistentBefore.indices.h})`);
  await menuAction(devtools, 'overview');

  await evaluate(devtools, `Reveal.toggleOverview()`);
  await waitUntil(devtools,
    `!document.querySelector('.reveal')?.classList.contains('overview') &&
      !document.documentElement.classList.contains('lcm-overview-active') &&
      !document.documentElement.classList.contains('lcm-viewport-settling')`,
    "overview mode did not restore the audience layout");
  const restoredOverview = await evaluate(devtools, `(() => {
    const stage = document.querySelector('.reveal .slides').getBoundingClientRect();
    const slide = Reveal.getCurrentSlide();
    const box = slide.getBoundingClientRect();
    const frame = slide.querySelector('iframe[data-lcm-src="/presentations/probe"]');
    return {
      warning: Boolean(document.querySelector('.lcm-deck-warning')),
      display: getComputedStyle(slide).display,
      frameVisible: frame && getComputedStyle(frame).visibility === 'visible' &&
        getComputedStyle(frame).display !== 'none',
      placeholderHidden: getComputedStyle(slide.querySelector('.lcm-live-placeholder')).display === 'none',
      retained: frame?.contentWindow?.__lcmPersistenceToken,
      stage: { left: stage.left, top: stage.top, width: stage.width, height: stage.height },
      slide: { left: box.left, top: box.top, width: box.width, height: box.height,
        transform: getComputedStyle(slide).transform },
      aligned: Math.abs(box.left - stage.left) <= 2 && Math.abs(box.top - stage.top) <= 2 &&
        Math.abs(box.width - stage.width) <= 2 && Math.abs(box.height - stage.height) <= 2
    };
  })()`);
  assert(!restoredOverview.warning && restoredOverview.display === "grid" &&
    restoredOverview.frameVisible && restoredOverview.placeholderHidden &&
    restoredOverview.retained === "retained" && restoredOverview.aligned,
    `overview exit did not restore the live slide: ${JSON.stringify(restoredOverview)}`);

  const printLifecycle = await evaluate(devtools, `(() => {
    window.dispatchEvent(new Event('beforeprint'));
    const frame = document.querySelector('iframe[data-lcm-src="/presentations/probe"]');
    const during = {
      printing: document.documentElement.classList.contains('lcm-printing'),
      hidden: getComputedStyle(frame).display === 'none',
      placeholder: getComputedStyle(frame.nextElementSibling).display !== 'none'
    };
    window.dispatchEvent(new Event('afterprint'));
    return {
      during,
      restored: !document.documentElement.classList.contains('lcm-printing') &&
        !document.documentElement.classList.contains('lcm-static-live') &&
        frame.contentWindow.__lcmPersistenceToken === 'retained'
    };
  })()`);
  assert(printLifecycle.during.printing && printLifecycle.during.hidden &&
    printLifecycle.during.placeholder && printLifecycle.restored,
    `browser print lifecycle is not reversible: ${JSON.stringify(printLifecycle)}`);

  await waitUntil(devtools,
    `!document.documentElement.classList.contains('lcm-viewport-settling')`,
    "deck did not settle after returning to the canvas slide");

  const clickTarget = await evaluate(devtools, `(() => {
    const frame = document.querySelector('iframe[data-lcm-src="/presentations/probe"]');
    const frameBox = frame.getBoundingClientRect();
    const canvas = frame.contentDocument.querySelector('.lc-presentation-probe canvas');
    const box = canvas.getBoundingClientRect();
    return {
      x: frameBox.left + box.left + box.width * 0.73,
      y: frameBox.top + box.top + box.height * 0.41,
      localX: box.width * 0.73,
      localY: box.height * 0.41
    };
  })()`);
  await devtools.command("Input.dispatchMouseEvent", {
    type: "mousePressed", x: clickTarget.x, y: clickTarget.y,
    button: "left", buttons: 1, clickCount: 1,
  });
  await devtools.command("Input.dispatchMouseEvent", {
    type: "mouseReleased", x: clickTarget.x, y: clickTarget.y,
    button: "left", buttons: 0, clickCount: 1,
  });
  await waitUntil(devtools,
    `Boolean(document.querySelector('iframe[data-lcm-src="/presentations/probe"]')?.contentDocument?.querySelector('.lc-presentation-probe')?.dataset.lastX)`,
    "visible canvas click did not reach the Bonito child");
  const clickResult = await evaluate(devtools, `(() => {
    const probe = document.querySelector('iframe[data-lcm-src="/presentations/probe"]')
      .contentDocument.querySelector('.lc-presentation-probe');
    return {x: Number(probe.dataset.lastX), y: Number(probe.dataset.lastY)};
  })()`);
  assert(Math.abs(clickResult.x - clickTarget.localX) <= 2,
    `canvas x coordinate drifted by ${Math.abs(clickResult.x - clickTarget.localX)} px`);
  assert(Math.abs(clickResult.y - clickTarget.localY) <= 2,
    `canvas y coordinate drifted by ${Math.abs(clickResult.y - clickTarget.localY)} px`);

  await devtools.command("Input.dispatchMouseEvent", { type: "mouseMoved",
    x: clickTarget.x, y: clickTarget.y });
  const pointerResult = await evaluate(devtools, `(() => {
    const pointer = document.querySelector('.lcm-pointer');
    return {x: parseFloat(pointer.style.left), y: parseFloat(pointer.style.top)};
  })()`);
  assert(Math.abs(pointerResult.x - clickTarget.x) <= 2 &&
    Math.abs(pointerResult.y - clickTarget.y) <= 2, "laser did not track inside the live iframe");
  // A focused canvas forwards L to the deck; a focused editor must keep it.
  await devtools.command("Input.dispatchKeyEvent", {
    type: "keyDown", key: "l", code: "KeyL", windowsVirtualKeyCode: 76,
  });
  await devtools.command("Input.dispatchKeyEvent", {
    type: "keyUp", key: "l", code: "KeyL", windowsVirtualKeyCode: 76,
  });
  assert(await evaluate(devtools, `Reveal.getIndices().h === ${persistentBefore.indices.h} &&
    !document.documentElement.classList.contains('lcm-pointer-enabled')`),
    "focused-canvas L failed to toggle, or navigated the deck");
  await evaluate(devtools, `(() => {
    const child = document.querySelector('iframe[data-lcm-src="/presentations/probe"]').contentDocument;
    const input = child.createElement('input'); input.id = 'test-editor';
    child.body.append(input); input.focus();
  })()`);
  await devtools.command("Input.dispatchKeyEvent", {
    type: "keyDown", key: "l", code: "KeyL", windowsVirtualKeyCode: 76, text: "l",
  });
  await devtools.command("Input.dispatchKeyEvent", {
    type: "keyUp", key: "l", code: "KeyL", windowsVirtualKeyCode: 76,
  });
  assert(await evaluate(devtools, `document.querySelector('iframe[data-lcm-src="/presentations/probe"]')
    .contentDocument.querySelector('#test-editor').value === 'l' &&
    Reveal.getIndices().h === ${persistentBefore.indices.h} &&
    !document.documentElement.classList.contains('lcm-pointer-enabled')`),
    "presentation hotkeys intercepted text input");
  await evaluate(devtools, `document.querySelector('iframe[data-lcm-src="/presentations/probe"]')
    .contentDocument.querySelector('#test-editor').remove(); document.body.focus()`);

  for (const theme of ["dark", "light"]) {
    await evaluate(devtools, `(() => {
      const select = document.querySelector('[data-lcm-theme-selector]');
      select.value = ${JSON.stringify(theme)}; select.dispatchEvent(new Event('change'));
    })()`);
    await waitUntil(devtools, `document.querySelector('iframe[data-lcm-src="/presentations/probe"]')
      .contentDocument.documentElement.dataset.lcmResolvedTheme === ${JSON.stringify(theme)}`,
      "theme selector did not synchronize the live child");
    assert(await evaluate(devtools, `(() => {
      const root = getComputedStyle(document.documentElement);
      const child = getComputedStyle(document.querySelector('iframe[data-lcm-src="/presentations/probe"]')
        .contentDocument.documentElement);
      return root.getPropertyValue('--lc-bg') === child.getPropertyValue('--lc-bg') &&
        document.documentElement.dataset.lcmResolvedTheme === ${JSON.stringify(theme)};
    })()`), "deck and live widget use different palettes");
  }

  // Exercise Chrome's real beforeprint/afterprint path, not just synthetic events.
  const directPdf = await devtools.command("Page.printToPDF", {
    printBackground: true, preferCSSPageSize: true,
  });
  if (process.env.LCM_PRESENTATION_ARTIFACTS) {
    await writeFile(process.env.LCM_PRESENTATION_ARTIFACTS + "/direct-print.pdf",
      Buffer.from(directPdf.data, "base64"));
  }
  await waitUntil(devtools, `!document.documentElement.classList.contains('lcm-printing') &&
    document.querySelector('iframe[data-lcm-src="/presentations/probe"]')
      .contentWindow.__lcmPersistenceToken === 'retained'`, "Chrome printing did not restore the live session");

  const audienceReport = {
    geometryChecks: geometry.length,
    viewports: [...new Set(geometry.map(sample =>
      sample.viewport.width + " × " + sample.viewport.height))],
    focus: focusAfter,
    persistence: persistentAfter,
    click: clickResult,
  };

  await menuAction(devtools, 'togglePdfExport');
  await waitUntil(devtools, `document.documentElement?.dataset.lcmDeckMode === 'static'`,
    'the actual PDF menu action did not open static preview');
  assert(await evaluate(devtools, `location.search.includes('lcm-print') &&
    !document.querySelector('.pdf-page')`), 'PDF menu invoked competing Reveal pagination');

  for (const mode of ["view=print", "print-pdf", "lcm-print", "receiver"]) {
    await navigate(devtools, `${deckUrl}?${mode}`);
    await waitUntil(devtools,
      `document.documentElement?.dataset.lcmDeckMode === 'static'`,
      `${mode} did not select static live mode`);
    const staticResult = await evaluate(devtools, `(() => {
      const frames = Array.from(document.querySelectorAll('.lcm-live-viewport iframe'));
      const placeholders = Array.from(document.querySelectorAll('.lcm-live-placeholder'));
      return {
        frames: frames.length,
        loaded: frames.filter(frame => frame.hasAttribute('src')).length,
        placeholders: placeholders.length,
        visible: placeholders.filter(node => getComputedStyle(node).display !== 'none').length,
        linked: placeholders.filter(node => node.querySelector('a[href]')).length,
      };
    })()`);
    assert(staticResult.frames > 0, `${mode} contains no live-boundary specimen`);
    assert(staticResult.loaded === 0, `${mode} loaded a live application`);
    assert(staticResult.visible === staticResult.frames, `${mode} did not show every placeholder`);
    assert(staticResult.linked === staticResult.frames, `${mode} placeholder is missing a link`);
    assert(await evaluate(devtools, `document.documentElement.dataset.lcmMathNotes === 'static' &&
      !document.querySelector('.lcm-math-callout, .lcm-math-term') &&
      [...document.querySelectorAll('.lcm-math-note')].every(note => note.hidden)`),
      `${mode} installed interactive equation explanations`);
    if (mode !== 'receiver') {
      assert(await evaluate(devtools, `(() => {
        const slides = [...document.querySelectorAll('.slides > section')];
        return !document.querySelector('.pdf-page') &&
          document.documentElement.scrollWidth <= innerWidth &&
          slides.length === 9 && slides.every(slide => {
            const r = slide.getBoundingClientRect();
            return r.left >= 0 && r.right <= innerWidth && r.width > 500 &&
              Math.abs(r.width / r.height - 16 / 9) < 0.01;
          });
      })()`), `${mode} clipped pages in the browser PDF preview`);
    }
  }

  await navigate(devtools, `${deckUrl}?lcm-print`);
  await devtools.command("Emulation.setEmulatedMedia", { media: "print" });
  const printStyle = await evaluate(devtools, `(() => {
    const resolveColor = token => {
      const probe = document.createElement('span');
      probe.style.color = 'var(' + token + ')';
      document.body.appendChild(probe);
      const color = getComputedStyle(probe).color;
      probe.remove();
      return color;
    };
    const heading = document.querySelector('.slides > section h2');
    const paragraph = document.querySelector('.slides > section p');
    const link = document.querySelector('.lcm-live-placeholder a');
    const reveal = document.querySelector('.reveal');
    const liveSlide = Array.from(document.querySelectorAll('.slides > section'))
      .find(slide => slide.querySelector('iframe[data-lcm-src="/widgets/control-panel"]'));
    const layout = liveSlide.querySelector('.lcm-layout');
    const slot = layout.querySelector('.lcm-slot');
    const frame = slot.querySelector('.lcm-live-viewport');
    const rect = node => {
      const box = node.getBoundingClientRect();
      return { x: box.x, y: box.y, width: box.width, height: box.height };
    };
    return {
      heading: getComputedStyle(heading).color,
      expectedHeading: resolveColor('--lc-heading'),
      text: getComputedStyle(paragraph).color,
      expectedText: resolveColor('--lc-text'),
      link: getComputedStyle(link).color,
      expectedLink: resolveColor('--lc-link'),
      background: getComputedStyle(reveal).backgroundColor,
      expectedBackground: (() => {
        const probe = document.createElement('span');
        probe.style.backgroundColor = 'var(--lc-bg)';
        document.body.appendChild(probe);
        const color = getComputedStyle(probe).backgroundColor;
        probe.remove();
        return color;
      })(),
      geometry: {
        slide: rect(liveSlide),
        heading: rect(liveSlide.querySelector('h2')),
        layout: rect(layout),
        slot: rect(slot),
        frame: rect(frame),
        sectionRows: getComputedStyle(liveSlide).gridTemplateRows,
        layoutColumns: getComputedStyle(layout).gridTemplateColumns
      }
    };
  })()`);
  assert(printStyle.heading === printStyle.expectedHeading &&
    printStyle.text === printStyle.expectedText &&
    printStyle.link === printStyle.expectedLink &&
    printStyle.background === printStyle.expectedBackground,
    `print media lost the authored palette: ${JSON.stringify(printStyle)}`);
  assert(printStyle.geometry.heading.height < 100 &&
    printStyle.geometry.layout.y -
      (printStyle.geometry.heading.y + printStyle.geometry.heading.height) < 40,
    `print grid stretched the heading row: ${JSON.stringify(printStyle.geometry)}`);
  assert(printStyle.geometry.frame.width <= printStyle.geometry.slot.width + 1 &&
    printStyle.geometry.frame.height <= printStyle.geometry.slot.height + 1,
    `live boundary escaped its print slot: ${JSON.stringify(printStyle.geometry)}`);
  // Clear emulation. Forcing "screen" also forces screen CSS during printToPDF.
  await devtools.command("Emulation.setEmulatedMedia", { media: "" });
  // Browser preview must also fit smaller screens, not only the large audit viewport.
  await setViewport(devtools, 1280, 720);
  assert(await evaluate(devtools, `document.documentElement.scrollWidth <= innerWidth &&
    [...document.querySelectorAll('.slides > section')].every(s =>
      s.getBoundingClientRect().right <= innerWidth && s.getBoundingClientRect().height <= innerHeight - 66)`),
    "PDF preview overflows a 1280px screen");
  if (process.env.LCM_PRESENTATION_ARTIFACTS) {
    const previewPdf = await devtools.command("Page.printToPDF", {
      printBackground: true, preferCSSPageSize: true,
    });
    await writeFile(process.env.LCM_PRESENTATION_ARTIFACTS + "/preview-print.pdf",
      Buffer.from(previewPdf.data, "base64"));
  }

  process.stdout.write(`${JSON.stringify({
    audience: audienceReport,
    overview: overviewResult,
    overviewRestored: restoredOverview,
    printLifecycle,
    printStyle,
    staticModes: ["view=print", "print-pdf", "lcm-print", "receiver"]
  }, null, 2)}\n`);
} finally {
  devtools.close();
}
