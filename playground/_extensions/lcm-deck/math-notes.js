/* Equation explanations: MathJax owns typesetting; Quarto's Popper owns placement.
 * Only explicitly authored terms are interactive. No widget/iframe traversal. */
(() => {
  'use strict';
  window.LCMMathNotes = {
    id: 'lcm-math-notes',
    init(deck) {
      const root = document.documentElement;
      const source = [...deck.getRevealElement().querySelectorAll('.lcm-math-note[data-lcm-math-target]')];
      if (!source.length) return;
      const staticMode = new URLSearchParams(location.search).has('lcm-print') ||
        new URLSearchParams(location.search).has('receiver');
      if (staticMode) { root.dataset.lcmMathNotes = 'static'; return; }

      const abort = new AbortController();
      const entries = source.map(note => ({ note, target: note.dataset.lcmMathTarget,
        origin: note.parentNode, next: note.nextSibling, slide: note.closest('.slides > section'), anchor: null }));
      const stage = deck.getSlidesElement();
      let active = null, popper = null, timer = 0, positionFrame = 0, loadingTimer = 0;
      let destroyed = false, hub = null, hook = null;
      const listen = (node, event, callback, options = {}) =>
        node.addEventListener(event, callback, { ...options, signal: abort.signal });

      const status = document.createElement('span');
      status.className = 'lcm-math-notes-status';
      status.setAttribute('role', 'status');
      status.setAttribute('aria-live', 'polite');
      document.querySelector('.lcm-deck-status')?.append(status);
      const report = (state, text) => {
        root.dataset.lcmMathNotes = state;
        status.textContent = text;
        status.title = state === 'ready' ? 'Click an outlined equation term to show its explanation' : text;
      };
      report('loading', 'Equation notes loading…');

      const panel = document.createElement('aside');
      panel.className = 'lcm-math-callout';
      // Avoid collisions with authored IDs, including on pages with other dialogs.
      panel.id = 'lcm-math-callout';
      while (document.getElementById(panel.id)) panel.id += '-popup';
      panel.hidden = true;
      panel.tabIndex = -1;
      panel.setAttribute('role', 'dialog');
      panel.setAttribute('aria-modal', 'false');
      panel.setAttribute('aria-label', 'Equation explanation');
      const closeButton = document.createElement('button');
      closeButton.type = 'button';
      closeButton.className = 'lcm-math-close';
      closeButton.textContent = '×';
      closeButton.setAttribute('aria-label', 'Close equation explanation');
      panel.append(closeButton);
      const leader = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
      leader.classList.add('lcm-math-leader');
      leader.setAttribute('aria-hidden', 'true');
      leader.setAttribute('hidden', '');
      const path = document.createElementNS(leader.namespaceURI, 'path');
      leader.append(path);
      document.body.append(leader, panel);

      function close(restoreFocus = false) {
        if (!active) return;
        const previous = active;
        active = null;
        popper?.destroy(); popper = null;
        cancelAnimationFrame(positionFrame);
        previous.anchor?.setAttribute('aria-expanded', 'false');
        previous.note.hidden = true;
        previous.origin.insertBefore(previous.note,
          previous.next?.parentNode === previous.origin ? previous.next : null);
        panel.hidden = true;
        leader.setAttribute('hidden', '');
        if (restoreFocus && previous.anchor?.isConnected) previous.anchor.focus({ preventScroll: true });
        else if (panel.contains(document.activeElement)) document.activeElement.blur();
      }

      function drawLeader({ state }) {
        if (!active) return;
        const a = active.anchor.getBoundingClientRect(), b = panel.getBoundingClientRect();
        const side = state.placement.split('-')[0];
        const clamp = (n, low, high) => Math.max(low, Math.min(n, high));
        let d;
        if (side === 'top' || side === 'bottom') {
          const x = (a.left + a.right) / 2;
          const y = side === 'top' ? a.top - 5 : a.bottom + 5;
          const by = side === 'top' ? b.bottom : b.top;
          const bx = clamp(x, b.left + 8, b.right - 8), mid = (y + by) / 2;
          const lip = side === 'top' ? 3 : -3;
          d = `M ${a.left} ${y + lip} V ${y} H ${a.right} V ${y + lip} M ${x} ${y} V ${mid} H ${bx} V ${by}`;
        } else {
          const y = (a.top + a.bottom) / 2;
          const x = side === 'left' ? a.left - 5 : a.right + 5;
          const bx = side === 'left' ? b.right : b.left;
          const by = clamp(y, b.top + 8, b.bottom - 8), mid = (x + bx) / 2;
          const lip = side === 'left' ? 3 : -3;
          d = `M ${x + lip} ${a.top} H ${x} V ${a.bottom} H ${x + lip} M ${x} ${y} H ${mid} V ${by} H ${bx}`;
        }
        leader.setAttribute('viewBox', `0 0 ${innerWidth} ${innerHeight}`);
        path.setAttribute('d', d);
        leader.removeAttribute('hidden');
      }

      function sizePanel() {
        const box = stage.getBoundingClientRect();
        panel.style.maxWidth = Math.max(0, box.width - 24) + 'px';
        panel.style.maxHeight = Math.max(0, box.height - 24) + 'px';
      }

      function updatePosition() {
        cancelAnimationFrame(positionFrame);
        positionFrame = requestAnimationFrame(() => {
          if (!active) return;
          if (!active.anchor.isConnected || active.slide !== deck.getCurrentSlide() ||
              getComputedStyle(active.anchor).visibility === 'hidden') { close(); return; }
          sizePanel(); popper?.update();
        });
      }

      function open(entry, keyboard = false) {
        if (deck.isOverview() || root.classList.contains('lcm-printing') || entry.slide !== deck.getCurrentSlide()) return;
        if (active === entry) { close(keyboard); return; }
        close();
        active = entry;
        panel.append(entry.note);
        entry.note.hidden = false;
        panel.hidden = false;
        entry.anchor.setAttribute('aria-expanded', 'true');
        sizePanel();
        popper = window.Popper.createPopper(entry.anchor, panel, {
          placement: 'top', strategy: 'fixed',
          modifiers: [
            { name: 'offset', options: { offset: [0, 20] } },
            { name: 'flip', options: { boundary: stage, padding: 12,
              fallbackPlacements: ['bottom', 'right', 'left'] } },
            { name: 'preventOverflow', options: { boundary: stage, padding: 12, altAxis: true, tether: false } },
            // Overlay coordinates stay in CSS pixels; no transformed slide ancestry.
            { name: 'computeStyles', options: { gpuAcceleration: false, adaptive: false } },
            { name: 'lcmLeader', enabled: true, phase: 'afterWrite', fn: drawLeader }
          ]
        });
        if (keyboard) panel.focus({ preventScroll: true });
      }

      function bind() {
        if (destroyed || !window.Popper?.createPopper) return;
        for (const entry of entries) {
          const matches = document.querySelectorAll('[id="' + entry.target + '"]');
          const anchor = matches?.length === 1 ? matches[0] : null;
          if (!anchor?.closest('.math') || anchor.closest('.slides > section') !== entry.slide || anchor === entry.anchor) continue;
          if (active === entry) close();
          entry.anchor = anchor;
          entry.note.hidden = true;
          entry.note.classList.remove('lcm-math-note-fallback');
          anchor.classList.add('lcm-math-term');
          anchor.setAttribute('role', 'button');
          anchor.setAttribute('tabindex', deck.isOverview() ? '-1' : '0');
          anchor.setAttribute('aria-label', 'Explain ' + anchor.textContent.trim());
          anchor.setAttribute('aria-haspopup', 'dialog');
          anchor.setAttribute('aria-controls', panel.id);
          anchor.setAttribute('aria-expanded', 'false');
          listen(anchor, 'click', event => {
            if (deck.isOverview() || root.classList.contains('lcm-printing')) return;
            if (event.target.closest('.lcm-math-term') !== anchor) return;
            event.stopImmediatePropagation();
            if (!window.getSelection()?.isCollapsed) return;
            event.preventDefault();
            open(entry, event.detail === 0);
          }, { capture: true });
          listen(anchor, 'keydown', event => {
            if (deck.isOverview() || root.classList.contains('lcm-printing')) return;
            if (event.target.closest('.lcm-math-term') !== anchor) return;
            if (event.key !== 'Enter' && event.key !== ' ') return;
            event.preventDefault(); event.stopImmediatePropagation();
            if (!event.repeat) open(entry, true);
          }, { capture: true });
        }
        if (entries.every(entry => entry.anchor?.isConnected)) {
          clearTimeout(timer);
          report('ready', 'Equation notes ready');
        }
        updatePosition();
      }

      listen(closeButton, 'click', () => close(true));
      listen(document, 'pointerdown', event => {
        if (active && !panel.contains(event.target) && !active.anchor.contains(event.target)) close();
      }, { capture: true });
      listen(window, 'blur', () => close()); // Includes moving focus into a live iframe.
      listen(window, 'keydown', event => {
        if (active && event.key === 'Escape') {
          event.preventDefault(); event.stopImmediatePropagation(); close(true);
        }
      }, { capture: true });
      listen(panel, 'keydown', event => {
        // The dialog is non-modal; Tab may leave it. Presentation keys belong to
        // its text/links while focused, not to Reveal's document listener.
        event.stopPropagation();
      });
      listen(panel, 'focusout', () => queueMicrotask(() => {
        if (active && !panel.contains(document.activeElement) && document.activeElement !== active.anchor) close();
      }));
      listen(window, 'beforeprint', () => close());
      listen(window, 'resize', updatePosition, { passive: true });
      listen(window, 'lcm:theme-changed', updatePosition);
      const observer = new ResizeObserver(updatePosition);
      observer.observe(stage); observer.observe(panel);
      const changed = () => { close(); hub?.Queue(bind); };
      const overview = () => {
        close();
        entries.forEach(entry => entry.anchor?.setAttribute('tabindex', deck.isOverview() ? '-1' : '0'));
      };
      const fragment = () => updatePosition();
      deck.on('slidechanged', changed);
      deck.on('overviewshown', overview);
      deck.on('overviewhidden', overview);
      deck.on('fragmenthidden', fragment);

      function connectMath() {
        if (destroyed || hub) return;
        clearTimeout(loadingTimer);
        if (!window.MathJax?.Hub?.Register || !window.Popper?.createPopper) {
          if (root.dataset.lcmMathNotes !== 'unavailable') loadingTimer = setTimeout(connectMath, 100);
          return;
        }
        hub = window.MathJax.Hub;
        // End Math covers initial processing, updates and explicit rerendering.
        hook = hub.Register.MessageHook('End Math', bind);
        hub.Queue(bind);
      }
      listen(document, 'load', event => {
        if (event.target.tagName === 'SCRIPT') connectMath();
      }, { capture: true });
      timer = setTimeout(() => {
        clearTimeout(loadingTimer);
        report('unavailable', 'Equation notes unavailable');
        // Preserve readable explanations if rendering/dependencies fail. Successful
        // late MathJax processing can still bind the terms without a reload.
        entries.filter(entry => !entry.anchor).forEach(entry => {
          entry.note.hidden = false;
          entry.note.classList.add('lcm-math-note-fallback');
        });
      }, 15000);
      connectMath();

      this.destroy = () => {
        close(); destroyed = true;
        clearTimeout(timer); clearTimeout(loadingTimer); cancelAnimationFrame(positionFrame);
        abort.abort(); observer.disconnect();
        if (hook) hub.signal.RemoveHook(hook);
        deck.off('slidechanged', changed); deck.off('overviewshown', overview);
        deck.off('overviewhidden', overview); deck.off('fragmenthidden', fragment);
        for (const entry of entries) {
          entry.anchor?.classList.remove('lcm-math-term');
          for (const name of ['role', 'tabindex', 'aria-label', 'aria-haspopup', 'aria-controls', 'aria-expanded'])
            entry.anchor?.removeAttribute(name);
        }
        panel.remove(); leader.remove(); status.remove();
        delete root.dataset.lcmMathNotes;
      };
      listen(window, 'pagehide', event => event.persisted ? close() : this.destroy());
    }
  };
})();
