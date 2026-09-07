// Exercise Quarto/Reveal's native fragments. No LCM sequencing controller.
export async function assertIncrementalLists({devtools, baseUrl, evaluate, navigate, waitUntil, assert}) {
  const read = expression => evaluate(devtools, expression);
  async function key(key, code) {
    for (const type of ['keyDown', 'keyUp']) await devtools.command('Input.dispatchKeyEvent', {
      type, key, code: key, windowsVirtualKeyCode: code,
    });
  }
  for (const [file, id, count, lists] of [
    ['starter', 'one-clear-argument', 4, 1],
    ['specimen', 'top-band-plus-two-supporting-regions', 6, 2],
  ]) {
    await navigate(devtools, `${baseUrl}/presentations/${file}.html`);
    await waitUntil(devtools, `document.documentElement.dataset.lcmDeckReady === 'true'`,
      `${file} did not initialize`);
    for (const theme of ['dark', 'light']) {
      const initial = await read(`(async () => {
        const select = document.querySelector('[data-lcm-theme-selector]');
        select.value = ${JSON.stringify(theme)}; select.dispatchEvent(new Event('change'));
        const slide = document.getElementById(${JSON.stringify(id)});
        const indices = Reveal.getIndices(slide);
        Reveal.slide(indices.h, 0, -1);
        document.activeElement?.blur();
        await document.fonts.ready;
        const fragments = [...slide.querySelectorAll('li.fragment')];
        return {
          count: fragments.length,
          lists: slide.querySelectorAll('.lcm-slot > ul, .lcm-slot > ol').length,
          nested: Boolean(slide.querySelector('li.fragment li.fragment')),
          numbered: Boolean(slide.querySelector('ol > li.fragment')),
          hidden: fragments.every(item => !item.classList.contains('visible')),
          boxes: fragments.map(item => {
            const r = item.getBoundingClientRect(); return [r.x, r.y, r.width, r.height];
          })
        };
      })()`);
      assert(initial.count === count && initial.lists === lists && initial.nested &&
        (lists === 1 || initial.numbered), `${file}: missing native nested/numbered list examples`);
      assert(initial.hidden && initial.boxes.every(box => box[2] > 0 && box[3] > 0),
        `${file}/${theme}: hidden fragments did not reserve their normal layout`);
      // Reveal fades fragments on return to a visited slide; visibility changes
      // at the end of that native transition, not at class removal.
      await waitUntil(devtools, `[...Reveal.getCurrentSlide().querySelectorAll('li.fragment')]
        .every(item => getComputedStyle(item).visibility === 'hidden')`,
        `${file}/${theme}: fragments did not finish hiding`);

      async function state(visible) {
        await waitUntil(devtools, `(() => {
          const slide = Reveal.getCurrentSlide();
          return slide.id === ${JSON.stringify(id)} &&
            [...slide.querySelectorAll('li.fragment')].every((item, index) =>
              item.classList.contains('visible') === (index < ${visible}));
        })()`, `${file}/${theme}: fragments did not follow source order at step ${visible}`);
        assert(await read(`(() => {
          const boxes = ${JSON.stringify(initial.boxes)};
          return [...Reveal.getCurrentSlide().querySelectorAll('li.fragment')].every((item, index) => {
            const r = item.getBoundingClientRect();
            return [r.x, r.y, r.width, r.height].every((v, axis) => Math.abs(v - boxes[index][axis]) < 1);
          });
        })()`), `${file}/${theme}: fragment reveal moved allocated content`);
      }
      // The built-in on-screen next button and real keyboard advance share the sequence.
      await read(`document.querySelector('.controls .navigate-right').click(); document.activeElement?.blur()`);
      await state(1);
      for (let n = 2; n <= count; n++) { await key('ArrowRight', 39); await state(n); }
      for (let n = count - 1; n >= 0; n--) { await key('ArrowLeft', 37); await state(n); }
      for (let n = 1; n <= count; n++) { await key('ArrowRight', 39); await state(n); }
      await key('ArrowRight', 39);
      assert(await read(`Reveal.getCurrentSlide().id !== ${JSON.stringify(id)}`),
        `${file}/${theme}: final advance stayed on the completed list`);
    }
    await navigate(devtools, `${baseUrl}/presentations/${file}.html?lcm-print`);
    await waitUntil(devtools, `document.documentElement.dataset.lcmDeckMode === 'static'`,
      `${file} print preview did not initialize`);
    assert(await read(`[...document.querySelectorAll('li.fragment')].every(item => {
      const style = getComputedStyle(item); return style.visibility === 'visible' && style.opacity === '1';
    })`), `${file}: print preview hides incremental content`);
  }
  console.log('Native incremental lists: nested bullets, two-list ordering, reverse, reserved geometry and print passed');
}
