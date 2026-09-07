(() => {
  'use strict';
  const root = document.documentElement;
  const query = new URLSearchParams(location.search);
  const printView = query.has('lcm-print');
  const receiver = query.has('receiver');
  const staticMode = printView || receiver;
  const frames = [];
  const readyFrames = new WeakSet();
  const connections = new WeakMap();
  const timeouts = new WeakMap();
  let deck, status, readiness, laserButton, pageNumber, pointer;
  let settling = false, printing = false, laser = false, overview = false;
  let generation = 0, settleTimer = 0, previousPrintClass = false;
  let pointerPosition = { x: innerWidth / 2, y: innerHeight / 2 };

  const slideFrames = slide => slide ? [...slide.querySelectorAll('.lcm-live-viewport iframe')] : [];
  const send = (frame, type, detail = {}) => frame.contentWindow?.postMessage(
    { namespace: 'lcm-deck', type, detail }, location.origin);
  const activeFrames = () => slideFrames(deck?.getCurrentSlide());
  const usable = () => !overview && !printing && !staticMode;

  function updateStatus() {
    if (!status) return;
    const active = activeFrames();
    const failed = active.some(frame => frame.dataset.lcmState === 'unavailable');
    const pending = active.filter(frame => !readyFrames.has(frame)).length;
    readiness.textContent = printView ? 'PDF preview · live views replaced by links' :
      receiver ? 'Presenter preview' : overview ? 'Slide overview' :
      !deck?.isReady() ? 'Loading presentation…' : failed ? 'Live view unavailable' :
      pending ? 'Loading live view…' : settling ? 'Resizing…' : 'Ready';
    const busy = !deck?.isReady() || (usable() && (pending || settling));
    status.dataset.state = failed && usable() ? 'error' : busy ? 'busy' : 'ready';
    laserButton.textContent = 'Laser ' + (laser ? 'on' : 'off') + ' · L';
    laserButton.setAttribute('aria-pressed', String(laser));
    laserButton.disabled = !usable();
    pageNumber.textContent = printView ? (deck?.getTotalSlides() ?? 0) + ' slides' :
      ((deck?.getIndices().h ?? 0) + 1) + ' / ' + (deck?.getTotalSlides() ?? 0);
    root.classList.toggle('lcm-pointer-enabled', laser && usable());
    frames.forEach(frame => {
      try { frame.contentDocument?.documentElement.classList.toggle('lcm-deck-laser', laser && usable()); }
      catch (_) { /* Failed routes may no longer be same-origin. */ }
    });
  }

  function movePointer(x, y) {
    pointerPosition = { x, y };
    if (!pointer) return;
    pointer.style.left = x + 'px';
    pointer.style.top = y + 'px';
  }

  function toggleLaser() {
    if (!usable()) return;
    laser = !laser;
    movePointer(pointerPosition.x, pointerPosition.y);
    updateStatus();
  }

  const editable = target => target?.closest?.(
    'input, select, textarea, button, a[href], [contenteditable="true"], [role="textbox"], [role="slider"], .lcm-math-term, .lcm-math-callout');
  function onKey(event) {
    if (event.ctrlKey || event.metaKey || event.altKey || event.isComposing) return;
    if (editable(event.target)) return;
    const key = event.key.toLowerCase();
    if (key === 'l' || key === 'e' || (printView && key === 'escape')) {
      event.preventDefault();
      event.stopImmediatePropagation();
      if (event.repeat) return;
      if (key === 'l') toggleLaser();
      else togglePdf();
    }
  }

  function syncTheme() {
    const theme = window.LineCableModelsTheme;
    if (!theme) return;
    frames.forEach(frame => send(frame, 'lcm:theme', {
      preference: theme.preference(), resolved: theme.resolved()
    }));
    deck?.getRevealElement().classList.toggle('has-light-background', theme.resolved() === 'light');
    deck?.getRevealElement().classList.toggle('has-dark-background', theme.resolved() === 'dark');
  }

  function connectFrame(frame) {
    const child = frame.contentWindow;
    try {
      if (!child?.document.querySelector('.lc-widget-shell')) return;
      if (connections.get(frame)?.document === child.document) return;
      connections.get(frame)?.abort.abort();
      const abort = new child.AbortController();
      connections.set(frame, { document: child.document, abort });
      child.addEventListener('keydown', onKey, { capture: true, signal: abort.signal });
      child.addEventListener('pointermove', event => {
        const box = frame.getBoundingClientRect();
        movePointer(box.left + event.clientX, box.top + event.clientY);
      }, { passive: true, signal: abort.signal });
      syncTheme();
    } catch (_) { /* Failed routes retain their fallback. */ }
  }

  function acknowledge(frame) {
    readyFrames.add(frame);
    clearTimeout(timeouts.get(frame));
    frame.dataset.lcmState = 'ready';
    frame.parentElement.classList.remove('lcm-live-loading', 'lcm-live-unavailable');
    connectFrame(frame);
    send(frame, activeFrames().includes(frame) && !overview ? 'lcm:slide-enter' : 'lcm:slide-leave');
    settle('live-ready');
    updateStatus();
  }

  function activate(slide) {
    if (staticMode || overview) return;
    slideFrames(slide).forEach(frame => {
      if (frame.hasAttribute('src')) return;
      const route = frame.dataset.lcmSrc;
      if (!route?.startsWith('/') || new URL(route, location.origin).origin !== location.origin) return;
      frame.dataset.lcmState = 'loading';
      frame.parentElement.classList.add('lcm-live-loading');
      const waiting = document.createElement('div');
      waiting.className = 'lcm-live-wait';
      waiting.textContent = 'Loading live view…';
      frame.parentElement.append(waiting);
      frame.addEventListener('load', () => {
        connectFrame(frame);
        send(frame, 'lcm:host-ready');
      });
      frame.src = route;
      timeouts.set(frame, setTimeout(() => {
        if (readyFrames.has(frame)) return;
        frame.dataset.lcmState = 'unavailable';
        frame.parentElement.classList.remove('lcm-live-loading');
        frame.parentElement.classList.add('lcm-live-unavailable');
        updateStatus();
      }, 15000));
    });
    updateStatus();
  }

  function checkGeometry() {
    if (!usable()) return;
    let invalid = false;
    activeFrames().filter(frame => readyFrames.has(frame)).forEach(frame => {
      for (let node = frame.parentElement; node && node !== root; node = node.parentElement) {
        const transform = getComputedStyle(node).transform;
        if (transform === 'none') continue;
        const matrix = new DOMMatrixReadOnly(transform);
        if (Math.abs(matrix.a - 1) > 0.001 || Math.abs(matrix.d - 1) > 0.001 ||
            Math.abs(matrix.b) > 0.001 || Math.abs(matrix.c) > 0.001) invalid = true;
      }
    });
    status.classList.toggle('lcm-layout-error', invalid);
    if (invalid) readiness.textContent = 'Layout error · reload presentation';
  }

  function settle(reason) {
    clearTimeout(settleTimer);
    const revision = ++generation;
    if (staticMode || printing || overview) return;
    settling = true;
    root.classList.add('lcm-viewport-settling');
    activeFrames().forEach(frame => send(frame, 'lcm:viewport-settling', { reason }));
    updateStatus();
    settleTimer = setTimeout(() => requestAnimationFrame(() => requestAnimationFrame(() => {
      if (revision !== generation || !usable()) return;
      settling = false;
      root.classList.remove('lcm-viewport-settling');
      activeFrames().forEach(frame => send(frame, 'lcm:viewport-settled', {
        width: frame.clientWidth, height: frame.clientHeight, devicePixelRatio
      }));
      updateStatus();
      checkGeometry();
    })), 220);
  }

  function setOverview(active) {
    overview = active;
    root.classList.toggle('lcm-overview-active', active);
    ++generation;
    clearTimeout(settleTimer);
    settling = false;
    root.classList.remove('lcm-viewport-settling');
    activeFrames().forEach(frame => send(frame, active ? 'lcm:slide-leave' : 'lcm:slide-enter'));
    if (!active) { activate(deck.getCurrentSlide()); settle('overviewhidden'); }
    updateStatus();
  }

  function enterPrint() {
    if (printing) return;
    if (overview) deck.toggleOverview(false);
    printing = true;
    previousPrintClass = root.classList.contains('print-pdf');
    root.classList.add('lcm-printing', 'print-pdf', 'lcm-static-live');
    frames.forEach(frame => send(frame, 'lcm:print-mode', { active: true }));
    updateStatus();
  }

  function leavePrint() {
    printing = false;
    root.classList.remove('lcm-printing');
    if (!previousPrintClass) root.classList.remove('print-pdf');
    if (!staticMode) root.classList.remove('lcm-static-live');
    frames.forEach(frame => send(frame, 'lcm:print-mode', { active: false }));
    settle('afterprint');
    updateStatus();
  }

  function togglePdf() {
    const url = new URL(location.href);
    url.searchParams.delete('view');
    url.searchParams.delete('print-pdf');
    if (printView) url.searchParams.delete('lcm-print');
    else url.searchParams.set('lcm-print', '');
    location.assign(url);
  }

  function playgroundHomeLink() {
    const link = document.createElement('a');
    link.href = '/';
    link.className = 'lcm-deck-home';
    link.dataset.action = 'home';
    link.title = 'Return to playground home';
    link.setAttribute('aria-label', link.title);
    link.innerHTML = '<svg viewBox="0 0 24 24" aria-hidden="true" focusable="false">' +
      '<path d="m3 10 9-7 9 7M5 9v12h14V9M9 21v-8h6v8"/></svg>';
    return link;
  }

  function installStatus() {
    status = document.createElement('footer');
    status.className = 'lcm-deck-status';
    status.setAttribute('aria-label', 'Presentation controls');
    readiness = document.createElement('output');
    readiness.setAttribute('role', 'status');
    readiness.setAttribute('aria-live', 'polite');
    readiness.textContent = 'Loading presentation…';
    status.append(readiness);
    const button = (text, action) => {
      const node = document.createElement('button');
      node.type = 'button'; node.textContent = text;
      node.addEventListener('click', action);
      status.append(node);
      return node;
    };
    laserButton = button('Laser off · L', toggleLaser);
    status.append(playgroundHomeLink());
    if (!staticMode) {
      button('Menu', () => deck?.getPlugin('menu')?.toggle()).dataset.action = 'menu';
    }
    button(printView ? 'Return to slides · E' : 'PDF preview · E', togglePdf);
    if (printView) button('Print / Save PDF', () => window.print());
    const label = document.createElement('label');
    label.textContent = 'Theme ';
    const selector = document.createElement('select');
    selector.className = 'lc-control-select';
    selector.dataset.lcmThemeSelector = '';
    selector.setAttribute('aria-label', 'Presentation color theme');
    ['system', 'dark', 'light'].forEach(value => selector.add(new Option(
      value[0].toUpperCase() + value.slice(1), value)));
    selector.value = window.LineCableModelsTheme?.preference() ?? 'dark';
    selector.addEventListener('change', () => window.LineCableModelsTheme.select(selector.value));
    label.append(selector); status.append(label);
    pageNumber = document.createElement('span');
    status.append(pageNumber);
    document.body.append(status);
    pointer = document.createElement('div');
    pointer.className = 'lcm-pointer'; pointer.setAttribute('aria-hidden', 'true');
    document.body.append(pointer);
    movePointer(pointerPosition.x, pointerPosition.y);
  }

  function start() {
    frames.push(...document.querySelectorAll('.lcm-live-viewport iframe[data-lcm-src]'));
    installStatus();
    if (staticMode) root.classList.add('lcm-static-live');
    window.addEventListener('keydown', onKey, true);
    document.addEventListener('pointermove', event => movePointer(event.clientX, event.clientY), { passive: true });
    window.addEventListener('message', event => {
      if (event.origin !== location.origin || event.data?.namespace !== 'lcm-deck' ||
          event.data.type !== 'lcm:child-ready') return;
      const frame = frames.find(frame => frame.contentWindow === event.source);
      if (frame) acknowledge(frame);
    });
    window.addEventListener('lcm:theme-changed', syncTheme);
    window.addEventListener('beforeprint', enterPrint);
    window.addEventListener('afterprint', leavePrint);
    ['resize', 'orientationchange', 'fullscreenchange'].forEach(name =>
      window.addEventListener(name, () => settle(name), { passive: true }));
    const initialize = () => {
      if (!window.Reveal?.isReady()) { setTimeout(initialize, 25); return; }
      deck = window.Reveal;
      deck.registerPlugin('lcm-math-notes', window.LCMMathNotes);
      // Controls receive their own keys; filter only at Reveal's dispatch boundary.
      deck.configure({
        keyboardCondition: event => !editable(event.target),
        scrollActivationWidth: 0,
        autoAnimate: false,
        hideInactiveCursor: false
      });
      root.dataset.lcmDeckMode = staticMode ? 'static' : 'audience';
      const pdf = deck.getPlugin('pdf-export');
      if (pdf) pdf.togglePdfExport = togglePdf;
      deck.addKeyBinding({ keyCode: 76, key: 'L', description: 'Toggle laser pointer' }, toggleLaser);
      deck.addKeyBinding({ keyCode: 69, key: 'E', description: 'PDF preview' }, togglePdf);
      // Keep Tools presentation-specific; playground navigation lives in the footer.
      // Scroll mode is unsupported, so remove its slot rather than repurposing it.
      document.querySelector('[onclick*="toggleScrollView"]')?.closest('li')?.remove();
      deck.removeKeyBinding(82);
      deck.on('slidechanged', event => {
        slideFrames(event.previousSlide).forEach(frame => send(frame, 'lcm:slide-leave'));
        activate(event.currentSlide);
        if (!overview) slideFrames(event.currentSlide).forEach(frame => send(frame, 'lcm:slide-enter'));
        settle('slidechanged'); updateStatus();
      });
      deck.on('overviewshown', () => setOverview(true));
      deck.on('overviewhidden', () => setOverview(false));
      if (printView) {
        deck.configure({ keyboard: false, overview: false, controls: false, progress: false });
        window.scrollTo(0, 0);
      }
      syncTheme(); activate(deck.getCurrentSlide()); settle('ready'); updateStatus();
      const attribution = document.querySelector('.reveal > .footer');
      if (attribution) {
        const label = document.createElement('span');
        label.className = 'lcm-deck-attribution';
        label.textContent = attribution.textContent.trim();
        status.insertBefore(label, laserButton);
      }
      root.dataset.lcmDeckReady = 'true';
    };
    initialize();
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start, { once: true });
  else start();
})();
