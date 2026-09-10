// The document shell is shared by every published HTML page, not live apps.
export async function assertPublishedShell({devtools, baseUrl, evaluate, navigate, waitUntil, setViewport, assert}) {
  const response=await fetch(new URL('/runtime/api/capabilities',baseUrl), {signal:AbortSignal.timeout(10000)});
  const capabilities=response.status===404 ? {ui_hosts:false} : await response.json();
  assert(typeof capabilities.ui_hosts==='boolean', 'Harness could not establish the publisher capability contract');
  const expectedEmbedState=capabilities.ui_hosts ? 'unavailable' : 'live';
  // This matrix measures the document shell and the embed rectangle, not the
  // inside of hundreds of repeatedly mounted Julia sessions. Real widget
  // rendering is covered by visual_contract_browser and catalogue_browser.
  // Keep real src/state/lazy attributes while giving matrix frames inert bodies.
  const geometryFixture=capabilities.ui_hosts ? null : await devtools.command('Page.addScriptToEvaluateOnNewDocument',{source:`(() => {
    const src=Object.getOwnPropertyDescriptor(HTMLIFrameElement.prototype,'src');
    Object.defineProperty(HTMLIFrameElement.prototype,'src',{...src,set(value) {
      if (this.hasAttribute('data-lc-published-src')) {
        this.srcdoc='<!doctype html><title>Document geometry fixture</title>';
      }
      src.set.call(this,value);
    }});
  })()`});
  const read = expression => evaluate(devtools, expression);
  const open = async path => {
    await navigate(devtools, new URL(path, baseUrl).href);
    assert(await read(`Boolean(document.querySelector('#quarto-sidebar'))`),
      `${path}: published page omitted the shared sidebar`);
    assert(await read(`Boolean(document.querySelector('#quarto-bootstrap')?.sheet?.cssRules.length)`),
      `${path}: compiled document stylesheet was not served`);
  };
  await setViewport(devtools, 1440, 800);
  await open('/');
  const paths = await read(`[...new Set([...document.querySelectorAll(
    '#quarto-sidebar a.sidebar-item-text[href]')].map(a => new URL(a.href))
    .filter(url => url.origin === location.origin && !url.hash)
    .map(url => url.pathname))]`);
  paths.push('/presentations/layouts.html', '/presentations/math-notes.html');
  // Public navigation is intentionally short. The complete developer sidebar
  // remains part of the same contour/scroll contract and must stay in coverage.
  await open('/dev/');
  const developerPaths = await read(`[...new Set([...document.querySelectorAll(
    '#quarto-sidebar a.sidebar-item-text[href]')].map(a => new URL(a.href))
    .filter(url => url.origin === location.origin && !url.hash)
    .map(url => url.pathname))]`);
  paths.push(...developerPaths.filter(path => !paths.includes(path)));

  async function inspectDocumentGeometry() {
    await waitUntil(devtools, `[...document.querySelectorAll('iframe[data-lc-published-src]')].every(frame => {
      const state=frame.closest('.lc-published-viewport')?.dataset.lcPreviewState;
      return state==='live' || state==='unavailable';
    })`, 'Published embeds did not resolve their shared live/inactive state');
    await read(`Promise.all(document.getAnimations().filter(animation =>
      Number.isFinite(animation.effect.getComputedTiming().endTime))
      .map(animation => animation.finished.catch(() => {})))`);
    const sample = await read(`(() => {
      window.scrollTo({top: 0, behavior: 'instant'});
      const main = document.querySelector('main#quarto-document-content');
      const css = getComputedStyle(main), rect = main.getBoundingClientRect();
      const measure=document.createElement('i');
      measure.style.cssText='position:fixed;visibility:hidden;display:block;width:'+css.maxWidth;
      document.body.append(measure); const maxWidth=measure.getBoundingClientRect().width; measure.remove();
      const box = node => { const r = node.getBoundingClientRect(); return {
        name: node.id || node.className || node.tagName, left:r.left, right:r.right, top:r.top, bottom:r.bottom}; };
      return {path:location.pathname, width:innerWidth, height:innerHeight,
        left:rect.left + parseFloat(css.paddingLeft), right:rect.right - parseFloat(css.paddingRight),
        header:box(document.querySelector('#quarto-header')),
        headerBorder:parseFloat(getComputedStyle(document.querySelector('#quarto-header')).borderBottomWidth),
        headerBorderColor:getComputedStyle(document.querySelector('#quarto-header')).borderBottomColor,
        nav:box(document.querySelector('.quarto-secondary-nav')),
        utilities:box(document.querySelector('.lc-publisher-theme-control')),
        title:box(main.querySelector('.quarto-title .title')),
        children:[...main.children].filter(node => node.getBoundingClientRect().height > 0 &&
          !['SCRIPT','STYLE'].includes(node.tagName)).map(box),
        context:document.querySelector('.quarto-secondary-nav-title').textContent.trim(),
        active:document.querySelector('#quarto-sidebar a.sidebar-link.active .menu-text')?.textContent.trim(),
        sidebar:box(document.querySelector('#quarto-sidebar')),
        main:box(main), footer:box(document.querySelector('footer.footer')),
        lastNavigation:box([...document.querySelectorAll('#quarto-sidebar .sidebar-item')].at(-1)),
        gutter:parseFloat(css.paddingLeft), maxWidth,
        availableWidth:document.querySelector('#quarto-content').getBoundingClientRect().width,
        logo:[...document.querySelectorAll('#quarto-sidebar img.sidebar-logo')].some(img => img.complete && img.naturalWidth > 0),
        headerBackground:getComputedStyle(document.querySelector('.quarto-secondary-nav')).backgroundColor,
        sidebarBackground:getComputedStyle(document.querySelector('#quarto-sidebar')).backgroundColor,
        scrollWidth:document.documentElement.scrollWidth, clientWidth:document.documentElement.clientWidth,
        scrollbarGutter:getComputedStyle(document.documentElement).scrollbarGutter,
        scrollHeight:document.documentElement.scrollHeight, clientHeight:document.documentElement.clientHeight,
        contentRight:document.querySelector('#quarto-content').getBoundingClientRect().right,
        escaped:[...document.querySelectorAll('body *')].filter(node => {
          const r=node.getBoundingClientRect(), s=getComputedStyle(node);
          return r.width>0 && r.height>0 && r.right>document.documentElement.clientWidth+1 && s.visibility==='visible';
        }).slice(0,5).map(box),
      };
    })()`);
    const where = `${sample.path} ${sample.width}×${sample.height}`;
    const available = sample.availableWidth - (sample.width > 900 ? sample.sidebar.right : 0);
    const expectedLeft = (sample.width > 900 ? sample.sidebar.right : 0) +
      Math.max(0, available - sample.maxWidth) / 2 + sample.gutter;
    assert(Math.abs(sample.left - expectedLeft) <= 1,
      `${where}: document was relocated by a descendant's Quarto page-full tracks: ${JSON.stringify({left:sample.left,expectedLeft,maxWidth:sample.maxWidth,available})}`);
    assert(sample.logo, `${where}: publisher omitted its application identity image`);
    assert(sample.headerBackground === sample.sidebarBackground,
      `${where}: compact/desktop header escaped the shared theme`);
    assert(sample.headerBorder === 1 && sample.headerBorderColor !== sample.headerBackground &&
      sample.nav.bottom <= sample.header.bottom - sample.headerBorder + 0.1 &&
      sample.utilities.bottom <= sample.nav.bottom + 0.1,
      `${where}: a header child paints over the bottom contour: ${JSON.stringify({header:sample.header,nav:sample.nav,utilities:sample.utilities})}`);
    assert(sample.scrollbarGutter === 'auto' && Math.abs(sample.header.right-sample.contentRight)<=1,
      `${where}: document/header reserves an extra empty scrollbar strip`);
    if (sample.scrollHeight <= sample.clientHeight) assert(Math.abs(sample.header.right-sample.width)<=1,
      `${where}: a page without vertical overflow still leaves a scrollbar gap`);
    const footerHeight=sample.footer.bottom-sample.footer.top;
    if (sample.width>900 && Math.max(sample.main.bottom,sample.lastNavigation.bottom+16)+footerHeight<=sample.height) {
      assert(sample.scrollHeight<=sample.clientHeight,
        `${where}: fitting contents still overflow because the shell guessed the footer height`);
    }
    if (sample.path==='/' && sample.width>=1920 && sample.height>=1047) {
      assert(sample.scrollHeight<=sample.clientHeight,
        `${where}: desktop landing page scrolls only to expose trailing spacing/footer`);
    }
    assert(sample.scrollWidth <= sample.clientWidth + 1, `${where}: content expands the document horizontally: ${JSON.stringify({scroll:sample.scrollWidth,client:sample.clientWidth,escaped:sample.escaped})}`);
    assert(Math.abs(sample.title.left - sample.left) <= 1, `${where}: h1 escaped the document column`);
    assert(sample.title.top >= sample.header.bottom + 8,
      `${where}: fixed header covers the title (short-screen header reservation)`);
    for (const child of sample.children) {
      assert(Math.abs(child.left - sample.left) <= 1 && Math.abs(child.right - sample.right) <= 1,
        `${where}: ${child.name} independently centres/shrinks/expands: ${JSON.stringify({column:[sample.left,sample.right],child})}`);
    }
    const canvas = await read(`(() => {
      const root=document.querySelector('.lc-cs-shell'); if (!root) return null;
      const figure=root.querySelector('.lc-cs-figure'), drawing=figure.querySelector('svg');
      return {figureHeight:figure.getBoundingClientRect().height, drawingHeight:drawing.getBoundingClientRect().height,
        figureDisplay:getComputedStyle(figure).display, overflow:getComputedStyle(root.querySelector('.lc-cs-canvas')).overflowY,
        compact:root.getBoundingClientRect().width <= 34 * parseFloat(getComputedStyle(document.documentElement).fontSize),
        state:root.dataset.sidebarState, workspaceWidth:root.querySelector('.lc-cs-workspace').getBoundingClientRect().width};
    })()`);
    if (canvas) {
      assert(canvas.figureDisplay === 'grid' && canvas.figureHeight >= 288 && canvas.drawingHeight >= 224 && canvas.overflow === 'auto',
        `${where}: Quarto displaced the owned canvas/figure layout: ${JSON.stringify(canvas)}`);
      if (canvas.compact) assert(canvas.state === 'collapsed' && canvas.workspaceWidth >= 200,
        `${where}: compact navigation crushes the work area`);
    }
    if (sample.active) assert(sample.context === sample.active,
      `${where}: breadcrumb describes a layout class rather than the declared navigation identity`);
    await inspectPublishedEmbeds(where);
  }

  async function inspectPublishedEmbeds(where, expectedState=expectedEmbedState) {
    const embeds=await read(`(() => {
      const box=n=>{const r=n.getBoundingClientRect();return {left:r.left,right:r.right,top:r.top,bottom:r.bottom,width:r.width,height:r.height};};
      return [...document.querySelectorAll('.lc-published-viewport')].map(root=>{
        const frame=root.querySelector('iframe'), note=root.querySelector('.lc-published-placeholder');
        const p=note.querySelector('p'), c=getComputedStyle(note), f=getComputedStyle(frame);
        return {title:frame.title,state:root.dataset.lcPreviewState,root:box(root),frame:box(frame),note:box(note),text:box(p),
          inset:['paddingTop','paddingRight','paddingBottom','paddingLeft'].map(k=>parseFloat(c[k])),
          framePadding:f.padding, frameHidden:frame.hidden, noteHidden:note.hidden, src:frame.getAttribute('src'),
          color:c.color,bg:getComputedStyle(root).backgroundColor,noteDisplay:c.display,
          height:parseFloat(getComputedStyle(root).height),textMargin:getComputedStyle(p).margin};
      });
    })()`);
    for (const embed of embeds) {
      assert(embed.state===expectedState,
        `${where}: ${embed.title} is ${embed.state}, but this publisher requires ${expectedState} previews`);
      assert(embed.root.width>0 && embed.root.height>0 && Math.abs(embed.root.height-embed.height)<=1,
        `${where}: ${embed.title} lost its reserved viewport`);
      if (embed.state==='unavailable') {
        assert(embed.frameHidden && !embed.src && !embed.noteHidden && embed.noteDisplay==='grid',
          `${where}: inactive ${embed.title} loaded a frame or lost its fallback`);
        assert(embed.inset.every(n=>n>=8 && n===embed.inset[0]) && embed.textMargin==='0px' &&
          embed.text.left>=embed.root.left+embed.inset[0]-1 && embed.text.right<=embed.root.right-embed.inset[0]+1 &&
          embed.text.top>=embed.root.top+embed.inset[0]-1 && embed.color!==embed.bg,
          `${where}: inactive ${embed.title} has unowned/clipped/unreadable inner content: ${JSON.stringify(embed)}`);
      } else {
        assert(!embed.frameHidden && embed.noteHidden && embed.framePadding==='0px' &&
          Math.abs(embed.frame.width-embed.root.width)<=1 && Math.abs(embed.frame.height-embed.root.height)<=1,
          `${where}: running ${embed.title} acquired extra wrapper spacing: ${JSON.stringify(embed)}`);
      }
    }
  }

  async function inspectSourceDialog() {
    await open('/templates/full-canvas.html');
    await read(`document.querySelector('#quarto-code-tools-source').click()`);
    await waitUntil(devtools, `!!document.querySelector('.modal.show')`, 'Source dialog did not open');
    // Bootstrap intentionally ignores dismissal while its opening transition
    // is in flight. Wait for the rendered dialog, not merely the show class.
    await read(`Promise.all(document.querySelector('.modal.show').getAnimations({subtree:true})
      .map(animation => animation.finished.catch(() => {})))`);
    const sample = await read(`(() => {
      const modal=document.querySelector('.modal.show'), panel=modal.querySelector('.modal-content');
      const probe=document.createElement('i'); panel.append(probe);
      const token=name=>{probe.style.color='var('+name+')';return getComputedStyle(probe).color;};
      const result={background:getComputedStyle(panel).backgroundColor, expected:token('--lc-panel-bg'),
        close:getComputedStyle(panel.querySelector('.btn-close')).color, muted:token('--lc-muted')};
      probe.remove(); return result;
    })()`);
    assert(sample.background === sample.expected && sample.close === sample.muted,
      `Source dialog retained compiled Darkly styles: ${JSON.stringify(sample)}`);
    await read(`document.querySelector('.modal.show .btn-close').click()`);
    await waitUntil(devtools, `!document.querySelector('.modal.show')`, 'Source dialog did not close');
  }

  async function toggleNavigation() {
    const point = await read(`(() => {
      const r = document.querySelector('.quarto-btn-toggle').getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    })()`);
    for (const type of ['mousePressed', 'mouseReleased']) await devtools.command('Input.dispatchMouseEvent', {
      type, ...point, button: 'left', buttons: type === 'mousePressed' ? 1 : 0, clickCount: 1,
    });
  }

  async function inspectNavigationStates() {
    const leaves = '.sidebar-item:not(.sidebar-item-section) > .sidebar-item-container > a.sidebar-link';
    const count = await read(`document.querySelectorAll('#quarto-sidebar ${leaves}').length`);
    // Exercise actual pointer and keyboard states, including root-level links.
    for (let index = 0; index < count; ++index) {
      const target = `document.querySelectorAll('#quarto-sidebar ${leaves}')[${index}]`;
      const snapshot = async () => read(`(() => {
        const link = ${target}, css = getComputedStyle(link);
        const token = name => {
          const probe = document.createElement('i'); probe.style.color = 'var(' + name + ')';
          link.append(probe); const value = getComputedStyle(probe).color; probe.remove(); return value;
        };
        return {label: link.textContent.trim(), color: css.color, background: css.backgroundColor,
          active: link.classList.contains('active'), border: css.borderLeftColor,
          hovered: link.matches(':hover'), focused: link.matches(':focus-visible'),
          text: token('--lc-text'), strong: token('--lc-strong-text'),
          hoverBg: token('--lc-hover-bg'), activeBg: token('--lc-active-bg'), accent: token('--lc-link')};
      })()`);
      const check = (s, interacting) => {
        assert(s.color === (s.active || interacting ? s.strong : s.text),
          `${s.label}: incorrect navigation text state ${JSON.stringify(s)}`);
        assert(s.background === (s.active ? s.activeBg : interacting ? s.hoverBg : 'rgba(0, 0, 0, 0)'),
          `${s.label}: incorrect navigation background state`);
        assert(s.border === (s.active ? s.accent : 'rgba(0, 0, 0, 0)'), `${s.label}: incorrect selection indicator`);
      };
      await read(`document.activeElement?.blur(); ${target}.scrollIntoView({block: 'center', behavior: 'instant'});
        new Promise(resolve => requestAnimationFrame(resolve))`);
      await devtools.command('Input.dispatchMouseEvent', {type: 'mouseMoved', x: 1430, y: 10});
      check(await snapshot(), false);
      const point = await read(`(() => { const r = ${target}.getBoundingClientRect();
        return {x: r.x + r.width / 2, y: r.y + r.height / 2}; })()`);
      await devtools.command('Input.dispatchMouseEvent', {type: 'mouseMoved', ...point});
      let state = await snapshot();
      assert(state.hovered, `${state.label}: pointer missed the navigation row`); check(state, true);
      await devtools.command('Input.dispatchMouseEvent', {type: 'mouseMoved', x: 1430, y: 10});
      check(await snapshot(), false);
      // Establish keyboard modality, then target each link without activating it.
      await devtools.command('Input.dispatchKeyEvent', {type: 'keyDown', key: 'Tab', code: 'Tab', windowsVirtualKeyCode: 9});
      await devtools.command('Input.dispatchKeyEvent', {type: 'keyUp', key: 'Tab', code: 'Tab', windowsVirtualKeyCode: 9});
      await read(`${target}.focus({preventScroll: true})`);
      state = await snapshot();
      assert(state.focused, `${state.label}: keyboard focus is not visible`); check(state, true);
      await read(`${target}.blur()`); check(await snapshot(), false);
    }
    await read(`window.scrollTo({top: 0, behavior: 'instant'})`);
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
        links: [...sidebar.querySelectorAll('.sidebar-item:not(.sidebar-item-section) > .sidebar-item-container > a.sidebar-link')]
          .map(link => {
            const css = getComputedStyle(link), r = link.getBoundingClientRect();
            return {text: link.textContent.trim(), display: css.display,
              padding: [css.paddingTop, css.paddingRight, css.paddingBottom, css.paddingLeft],
              fontSize: css.fontSize, lineHeight: css.lineHeight, border: css.borderLeftWidth,
              left: r.left, width: r.width, parentWidth: link.parentElement.getBoundingClientRect().width};
          }),
      };
    })()`);
    const where = `${await read('location.pathname')} at ${sample.width}px`;
    const text = await read(`(() => {
      const main = document.querySelector('main#quarto-document-content');
      const probe = document.createElement('i'); probe.style.color = 'var(--lc-heading)'; main.append(probe);
      const heading = getComputedStyle(probe).color; probe.remove();
      return {heading, headings: [...main.querySelectorAll('.quarto-title-block .title, section.level1 > h1, section.level2 > h2, section.level3 > h3, section.level4 > h4, section.level5 > h5, section.level6 > h6')]
        .map(node => ({text: node.textContent.trim(), color: getComputedStyle(node).color}))};
    })()`);
    for (const heading of text.headings) assert(heading.color === text.heading,
      `${where}: heading ${JSON.stringify(heading.text)} escaped the shared theme`);
    assert(sample.position.every(value => value === 'static'), `${where}: shell escaped document flow`);
    assert(!sample.horizontalOverflow, `${where}: horizontal page overflow`);
    assert(sample.footer.top >= sample.main.bottom - 1, `${where}: footer overlays page content`);
    if (expanded) {
      assert(sample.links.length > 0, `${where}: no navigation leaves checked`);
      const reference = sample.links[0];
      for (const link of sample.links) {
        const row = `${where}: navigation row ${JSON.stringify(link.text)}`;
        assert(link.display === 'block' && Math.abs(link.width - link.parentWidth) <= 1,
          `${row} is not a full-width navigation target`);
        assert(link.padding.every(value => parseFloat(value) > 0) && link.border === '3px',
          `${row} lost its spacing or active-indicator reservation`);
        assert(JSON.stringify(link.padding) === JSON.stringify(reference.padding) &&
          link.fontSize === reference.fontSize && link.lineHeight === reference.lineHeight &&
          Math.abs(link.left - reference.left) <= 1,
          `${row} differs from its peers in sizing or alignment`);
      }
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
    for (const path of ['/', '/dev/']) { await open(path); await inspectNavigationStates(); }
    for (const path of paths) {
      await open(path);
      assert(await read(`document.documentElement.dataset.lcmResolvedTheme === ${JSON.stringify(theme)}`),
        `${path}: shell did not adopt the ${theme} theme`);
      await inspect();
    }
    await inspectSourceDialog();
    // Include wide desktops (where conflicting max-widths first diverge),
    // short projection screens and small phones. Compare actual document
    // edges, not only colours or the outer shell's lack of overflow.
    for (const [width, height] of [[2560,1440], [1920,1080], [1920,1047], [1440,800], [1024,600], [901,640], [900,640], [768,640], [390,844]]) {
      await setViewport(devtools, width, height);
      for (const path of paths) { await open(path); await inspectDocumentGeometry(); }
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
  // Error-state layout must not depend on whether this harness is a standalone
  // publisher or a runtime gateway. Only the capability response is substituted;
  // the page, shortcode, script, theme and CSS are the real published assets.
  for (const failure of ['offline','503','invalid-json']) {
    const {identifier}=await devtools.command('Page.addScriptToEvaluateOnNewDocument',{source:`(() => {
      const original=window.fetch;
      window.fetch=(url, options)=>new URL(url,location.href).pathname==='/runtime/api/capabilities' ?
        ${failure==='offline' ? "Promise.reject(new TypeError('offline fixture'))" :
          `Promise.resolve(new Response('${failure==='invalid-json' ? '<html>unavailable</html>' : '{}'}', {status:${failure==='503' ? 503 : 200}}))`} : original(url,options);
    })()`});
    try {
      for (const theme of ['dark','light']) {
        await read(`localStorage.setItem('lcm.playground.theme', ${JSON.stringify(theme)})`);
        await setViewport(devtools, 390, 844);
        await open('/templates/bonito-widget.html');
        await waitUntil(devtools, `document.querySelector('.lc-published-viewport')?.dataset.lcPreviewState==='unavailable'`,
          `${failure}: preview escaped the shared unavailable state`);
        assert(await read(`document.querySelector('.lc-published-placeholder p').textContent.includes('service is unavailable')`),
          `${failure}: missing actionable service-unavailable message`);
        await inspectPublishedEmbeds(`${failure} / ${theme}`, 'unavailable');
      }
    } finally { await devtools.command('Page.removeScriptToEvaluateOnNewDocument',{identifier}); }
  }
  await setViewport(devtools, 1920, 1080);
  if (geometryFixture) await devtools.command('Page.removeScriptToEvaluateOnNewDocument',geometryFixture);
  console.log(`Published document shell: ${paths.length} routes, both themes, 9 viewport sizes; fitting-page height, alignment, painted header contour, scrollbar edge, live/inactive embed insets and desktop/compact flow passed`);
}
