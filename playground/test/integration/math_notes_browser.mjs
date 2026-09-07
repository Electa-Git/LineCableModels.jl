/* Real-browser contract, shared by the full presentation gate and focused checks. */
export async function assertMathNotes({ devtools, evaluate, waitUntil, click, setViewport, assert }) {
  const run = expression => evaluate(devtools, expression);
  const wait = (expression, message) => waitUntil(devtools, expression, message, 30000);
  const visible = `!document.querySelector('.lcm-math-callout').hidden`;
  const key = async (key, code, windowsVirtualKeyCode) => {
    for (const type of ['keyDown', 'keyUp']) await devtools.command('Input.dispatchKeyEvent', {
      type, key, code, windowsVirtualKeyCode
    });
  };
  const bounded = async () => {
    await wait(`(() => {
      const panel = document.querySelector('.lcm-math-callout');
      const p = panel.getBoundingClientRect(), s = Reveal.getSlidesElement().getBoundingClientRect();
      return !panel.hidden && p.width > 40 && p.height > 20 &&
        p.left >= s.left && p.top >= s.top && p.right <= s.right && p.bottom <= s.bottom &&
        !document.querySelector('.lcm-math-leader').hasAttribute('hidden');
    })()`, 'equation explanation/leader did not fit the slide');
  };
  await run(`Reveal.slide(Reveal.getIndices(document.querySelector('#balanced-evidence')).h)`);
  await wait(`document.documentElement.dataset.lcmMathNotes === 'ready'`, 'equation notes did not become ready');
  assert(await run(`document.querySelectorAll('.lcm-math-term').length === 2 &&
    !document.querySelector('.MathJax_Error') &&
    [...document.querySelectorAll('.lcm-math-note')].every(note => note.hidden)`),
  'equation markup failed to render or exposed explanation text before activation');
  const before = await run(`JSON.stringify({ index: Reveal.getIndices(),
    rect: document.querySelector('#impedance-imag').getBoundingClientRect().toJSON(),
    math: document.querySelector('#balanced-evidence script[type^="math/tex"]').textContent })`);

  for (const theme of ['dark', 'light']) {
    await run(`LineCableModelsTheme.select('${theme}')`);
    await click(devtools, '#impedance-imag');
    await bounded();
    assert(await run(`(() => {
      const panel = document.querySelector('.lcm-math-callout');
      const sample = document.createElement('span');
      sample.style.color = 'var(--lc-text)'; sample.style.background = 'var(--lc-panel-bg)';
      document.body.append(sample);
      const matches = getComputedStyle(sample).color === getComputedStyle(panel).color &&
        getComputedStyle(sample).backgroundColor === getComputedStyle(panel).backgroundColor;
      sample.remove();
      return matches && panel.textContent.includes('Imaginary part of impedance') &&
        getComputedStyle(panel).caretColor === 'rgba(0, 0, 0, 0)' &&
        getComputedStyle(panel).userSelect === 'text' &&
        document.querySelector('#impedance-imag').getAttribute('aria-expanded') === 'true';
    })()`), `math callout violates shared ${theme} appearance or selectable read-only text`);
    await click(devtools, '#admittance-imag');
    await bounded();
    assert(await run(`document.querySelector('#impedance-imag').getAttribute('aria-expanded') === 'false' &&
      document.querySelector('.lcm-math-callout').textContent.includes('Imaginary part of admittance') &&
      document.querySelectorAll('.lcm-math-callout .lcm-math-note').length === 1`), 'switching terms left a stale note');
    await click(devtools, '#admittance-imag');
    assert(await run(`!(${visible})`), 'clicking the same term did not close its note');
  }
  assert(await run(`JSON.stringify({ index: Reveal.getIndices(),
    rect: document.querySelector('#impedance-imag').getBoundingClientRect().toJSON(),
    math: document.querySelector('#balanced-evidence script[type^="math/tex"]').textContent })`) === before,
    'opening explanations changed slide position, equation source, or equation geometry');

  await click(devtools, '#impedance-imag');
  await click(devtools, '#balanced-evidence h2');
  assert(await run(`!(${visible})`), 'outside click did not dismiss the explanation');
  await run(`document.querySelector('#impedance-imag').focus()`);
  await key('Enter', 'Enter', 13);
  await bounded();
  assert(await run(`document.activeElement.classList.contains('lcm-math-callout')`),
    'keyboard activation did not focus the non-modal explanation');
  await key('l', 'KeyL', 76);
  await key(' ', 'Space', 32);
  assert(await run(`Reveal.getCurrentSlide().id === 'balanced-evidence' && ${visible} &&
    !document.documentElement.classList.contains('lcm-pointer-enabled')`),
    'focused callout keys escaped into presentation controls');
  await key('Escape', 'Escape', 27);
  assert(await run(`!(${visible}) && !Reveal.isOverview() &&
    document.activeElement.id === 'impedance-imag'`), 'Escape did not close and restore focus without opening overview');
  await key(' ', 'Space', 32);
  await bounded();
  await click(devtools, '.lcm-math-close');
  assert(await run(`!(${visible})`), 'close button did not dismiss note');

  assert(await run(`(() => {
    const anchor = document.querySelector('#impedance-imag');
    const range = document.createRange(); range.selectNodeContents(anchor);
    const selection = getSelection(); selection.removeAllRanges(); selection.addRange(range);
    anchor.click();
    const unchanged = document.querySelector('.lcm-math-callout').hidden && !selection.isCollapsed;
    selection.removeAllRanges(); return unchanged;
  })()`), 'selecting equation text activated a note or lost its selection');

  await click(devtools, '#impedance-imag');
  for (const [width, height] of [[1280, 720], [800, 600], [1920, 1080]]) {
    await setViewport(devtools, width, height);
    await bounded();
  }
  // Move only the annotated term to each stage corner. This exercises Popper's
  // edge handling without adding slide-level transforms or touching a live view.
  for (const corner of ['top-left', 'top-right', 'bottom-left', 'bottom-right']) {
    await run(`(() => {
      const anchor = document.querySelector('#impedance-imag');
      anchor.style.cssText = ''; anchor.style.display = 'inline-block';
      const a = anchor.getBoundingClientRect(), s = Reveal.getSlidesElement().getBoundingClientRect();
      const x = '${corner}'.includes('left') ? s.left + 16 : s.right - a.width - 16;
      const y = '${corner}'.includes('top') ? s.top + 16 : s.bottom - a.height - 16;
      anchor.style.transform = 'translate(' + (x - a.left) + 'px,' + (y - a.top) + 'px)';
      window.dispatchEvent(new Event('resize'));
    })()`);
    await bounded();
  }
  await run(`document.querySelector('#impedance-imag').style.cssText = ''; window.dispatchEvent(new Event('resize'))`);
  await bounded();
  await run(`Reveal.toggleOverview(true)`);
  assert(await run(`!(${visible})`), 'overview retained a floating equation note');
  assert(await run(`document.querySelector('#impedance-imag').tabIndex === -1`),
    'overview retained an interactive term tab stop');
  await click(devtools, '#impedance-imag');
  await wait(`!Reveal.isOverview() && document.querySelector('#impedance-imag').tabIndex === 0`,
    'an equation term intercepted thumbnail selection or did not restore keyboard access');
  await click(devtools, '#impedance-imag');
  await run(`window.dispatchEvent(new Event('beforeprint'))`);
  assert(await run(`!(${visible}) && getComputedStyle(document.querySelector('#impedance-imag')).outlineStyle === 'none'`),
    'print retained an explanation or its interactive outline');
  await run(`window.dispatchEvent(new Event('afterprint'))`);
  await click(devtools, '#impedance-imag');
  await run(`Reveal.next()`);
  assert(await run(`!(${visible})`), 'slide navigation left a stale note');
  await run(`Reveal.prev()`);

  await run(`new Promise(resolve => MathJax.Hub.Queue(['Rerender', MathJax.Hub,
    document.querySelector('#balanced-evidence')], resolve))`);
  await wait(`document.querySelector('#impedance-imag')?.classList.contains('lcm-math-term')`,
    'MathJax rerender lost the term binding');
  await click(devtools, '#impedance-imag');
  await bounded();
  await click(devtools, '#impedance-imag');
  assert(await run(`!(${visible})`), 'rerender installed duplicate toggle handlers');

  // Explicit teardown must remove handlers, portals and the MathJax hook. A late
  // dependency recovery must hide fallback prose again and bind exactly once.
  await run(`window.__notesPopper = window.Popper;
    Reveal.getPlugin('lcm-math-notes').destroy(); window.Popper = null;
    Reveal.getPlugin('lcm-math-notes').init(Reveal)`);
  await wait(`document.documentElement.dataset.lcmMathNotes === 'unavailable'`, 'missing Popper stayed loading forever');
  assert(await run(`[...document.querySelectorAll('.lcm-math-note')].every(note => !note.hidden) &&
    document.querySelectorAll('.lcm-math-term').length === 0`), 'dependency failure silently lost explanation content');
  await run(`window.Popper = window.__notesPopper; delete window.__notesPopper;
    const script = document.createElement('script'); document.head.append(script);
    script.dispatchEvent(new Event('load')); script.remove()`);
  await wait(`document.documentElement.dataset.lcmMathNotes === 'ready'`, 'late dependency recovery did not bind notes');
  await click(devtools, '#impedance-imag');
  await bounded();
  await click(devtools, '#impedance-imag');
  assert(await run(`!(${visible}) && document.querySelectorAll('.lcm-math-callout').length === 1 &&
    [...document.querySelectorAll('.lcm-math-note')].every(note => note.hidden)`),
    'reinitialization leaked callouts or fallback text');
  await run(`document.activeElement?.blur()`);
  console.log('Equation note browser contract passed: themes, geometry, input, lifecycle, rerender and failure recovery');
}
