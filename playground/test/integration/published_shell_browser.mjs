// The document shell is shared by every published HTML page, not live apps.
export async function assertPublishedShell({devtools, baseUrl, evaluate, navigate, waitUntil, setViewport, assert}) {
  const read = expression => evaluate(devtools, expression);
  const open = async path => {
    await navigate(devtools, new URL(path, baseUrl).href);
    assert(await read(`Boolean(document.querySelector('#quarto-sidebar'))`),
      `${path}: published page omitted the shared sidebar`);
  };
  await setViewport(devtools, 1440, 800);
  await open('/');
  const paths = await read(`[...new Set([...document.querySelectorAll(
    '#quarto-sidebar a.sidebar-item-text[href]')].map(a => new URL(a.href))
    .filter(url => url.origin === location.origin && !url.hash)
    .map(url => url.pathname))]`);
  paths.push('/presentations/layouts.html', '/presentations/math-notes.html');

  async function toggleNavigation() {
    const point = await read(`(() => {
      const r = document.querySelector('.quarto-btn-toggle').getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    })()`);
    for (const type of ['mousePressed', 'mouseReleased']) await devtools.command('Input.dispatchMouseEvent', {
      type, ...point, button: 'left', buttons: type === 'mousePressed' ? 1 : 0, clickCount: 1,
    });
  }

  async function inspect(compact = false, expanded = true) {
    await read(`document.fonts.ready.then(() => window.scrollTo({top: 0, behavior: 'instant'}))`);
    const sample = await read(`(() => {
      const sidebar = document.querySelector('#quarto-sidebar');
      const menu = sidebar.querySelector('.sidebar-menu-container');
      const main = document.querySelector('main#quarto-document-content');
      const footer = document.querySelector('footer.footer');
      const box = node => { const r = node.getBoundingClientRect();
        return {top: r.top, bottom: r.bottom, left: r.left, right: r.right}; };
      return {
        sidebar: box(sidebar), main: box(main), footer: box(footer),
        position: [sidebar, main, footer].map(node => getComputedStyle(node).position),
        overflow: [sidebar, menu, main].map(node => getComputedStyle(node).overflowY),
        borders: [sidebar, footer].map(node => getComputedStyle(node).borderRightWidth),
        menuFits: menu.scrollHeight <= menu.clientHeight + 1,
        lastLink: box([...sidebar.querySelectorAll('.sidebar-item-text')].at(-1)),
        hidden: getComputedStyle(sidebar).display === 'none',
        width: innerWidth, height: innerHeight,
        documentHeight: document.documentElement.scrollHeight,
        horizontalOverflow: document.documentElement.scrollWidth > innerWidth + 1,
      };
    })()`);
    const where = `${await read('location.pathname')} at ${sample.width}px`;
    assert(sample.position.every(value => value === 'static'), `${where}: shell escaped document flow`);
    assert(!sample.horizontalOverflow, `${where}: horizontal page overflow`);
    assert(sample.footer.top >= sample.main.bottom - 1, `${where}: footer overlays page content`);
    if (expanded) {
      assert(!sample.hidden && sample.menuFits, `${where}: navigation has its own clipped scroll area`);
      assert(sample.overflow.every(value => value === 'visible'), `${where}: nested document scroll owner`);
      assert(sample.sidebar.bottom >= sample.lastLink.bottom, `${where}: navigation ends before its links`);
      assert(Math.abs(sample.sidebar.bottom - sample.footer.top) <= 1 || compact,
        `${where}: sidebar/footer contour is interrupted`);
    } else assert(sample.hidden, `${where}: compact navigation did not close`);
    if (compact) {
      assert(sample.main.left >= 0 && sample.main.right <= sample.width + 1,
        `${where}: hidden sidebar still reserves page width`);
      if (expanded) assert(sample.main.top >= sample.sidebar.bottom,
        `${where}: expanded compact menu overlays content`);
      assert(sample.footer.left === 0 && sample.footer.right >= sample.width - 20,
        `${where}: compact footer is missing or incorrectly sized`);
    } else {
      assert(sample.sidebar.top === 0 && sample.main.left >= sample.sidebar.right,
        `${where}: page overlaps the sidebar`);
      assert(sample.borders.every(value => value === '1px') &&
        Math.abs(sample.sidebar.right - sample.footer.right) <= 1,
        `${where}: sidebar/footer contours differ`);
    }
    // Scroll to the real page end. No fixed footer or nested main should hide it.
    await read(`window.scrollTo({top: document.documentElement.scrollHeight, behavior: 'instant'})`);
    assert(await read(`Math.abs(document.querySelector('footer.footer').getBoundingClientRect()
      .bottom - innerHeight) <= 1`), `${where}: footer cannot be reached by document scrolling`);
    await read(`window.scrollTo({top: 0, behavior: 'instant'})`);
  }

  for (const theme of ['dark', 'light']) {
    await read(`localStorage.setItem('lcm.playground.theme', ${JSON.stringify(theme)})`);
    await setViewport(devtools, 1440, 800);
    for (const path of paths) {
      await open(path);
      assert(await read(`document.documentElement.dataset.lcmResolvedTheme === ${JSON.stringify(theme)}`),
        `${path}: shell did not adopt the ${theme} theme`);
      await inspect();
    }
    await setViewport(devtools, 768, 640);
    for (const path of paths) { await open(path); await inspect(true, false); }
    for (const width of [1024, 901, 900, 768, 390]) {
      await setViewport(devtools, width, 640);
      await open('/');
      if (width > 900) { await inspect(); continue; }
      await inspect(true, false);
      await toggleNavigation();
      await waitUntil(devtools, `document.querySelector('#quarto-sidebar').classList.contains('show')`,
        'compact menu did not open');
      await inspect(true, true);
      await toggleNavigation();
      await waitUntil(devtools, `!document.querySelector('#quarto-sidebar').matches('.show, .collapsing')`,
        'compact menu did not close');
      await inspect(true, false);
    }
  }
  await setViewport(devtools, 1920, 1080);
  console.log(`Published document shell: ${paths.length} routes, both themes, desktop/compact flow passed`);
}
