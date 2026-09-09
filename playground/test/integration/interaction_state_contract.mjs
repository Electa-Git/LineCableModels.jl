import assert from 'node:assert/strict';

// Exercise actual pointer/keyboard transitions, not pre-painted state classes.
// Runs against real template/widget/workbench routes on the broker-free fixture.
export async function checkInteractionStates({base, command, evaluate, wait, selectTheme}) {
  const report = [];
  const q = JSON.stringify;
  const settle = () => evaluate(`(async () => {
    const frame = () => new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)));
    await frame();
    // Compare settled states, not a color halfway through a legitimate fade.
    // Never wait for indefinite busy/spinner animations.
    await Promise.all(document.getAnimations()
      .filter(a => Number.isFinite(a.effect.getComputedTiming().endTime))
      .map(a => a.finished.catch(() => {})));
    await frame();
  })()`);
  const away = async () => {
    await command('Input.dispatchMouseEvent', {type: 'mouseMoved', x: 1599, y: 999});
    await settle();
  };
  const pointAt = async selector => {
    const point = await evaluate(`(() => {
      const el = document.querySelector(${q(selector)});
      el.scrollIntoView({block: 'nearest', inline: 'nearest'});
      const r = el.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    })()`);
    await command('Input.dispatchMouseEvent', {type: 'mouseMoved', ...point});
    // CDP acknowledges input before the next rendered hover frame. Do not
    // compare a genuine hover against the preceding unhovered frame.
    await settle();
    return point;
  };
  const click = async selector => {
    const point = await pointAt(selector);
    await command('Input.dispatchMouseEvent', {type: 'mousePressed', button: 'left', buttons: 1, clickCount: 1, ...point});
    await command('Input.dispatchMouseEvent', {type: 'mouseReleased', button: 'left', buttons: 0, clickCount: 1, ...point});
    await settle();
  };
  const navigate = async (route, selector) => {
    await command('Page.navigate', {url: base + route});
    await wait(`location.pathname === ${q(route)} && document.readyState === 'complete' &&
      Boolean(document.querySelector(${q(selector)}))`, 'route did not mount: ' + route);
    await settle();
    // The template fixture also serves the visual/X-ray audit. Picking mode
    // deliberately intercepts clicks, so turn it off through its real launcher.
    await evaluate(`(() => {
      const toggle = document.querySelector('.lc-xray-host')?.shadowRoot?.querySelector('.xray-toggle');
      if (toggle?.getAttribute('aria-pressed') === 'true') toggle.click();
    })()`);
    await away();
  };
  const theme = async value => {
    await selectTheme(value);
    await wait(`document.documentElement.dataset.lcmResolvedTheme === ${q(value)}`, 'theme did not update');
    await settle();
  };
  const style = selector => evaluate(`(() => {
    const el = document.querySelector(${q(selector)}), s = getComputedStyle(el);
    return {color:s.color, bg:s.backgroundColor, weight:s.fontWeight,
      left:s.borderLeftColor, bottom:s.borderBottomColor, shadow:s.boxShadow};
  })()`);
  const selected = (selector, marker) => `document.querySelector(${q(selector)}).matches(${q(marker)})`;
  const expectSelected = async (selector, marker) => wait(selected(selector, marker), 'selection did not update: ' + selector);
  const key = async (key, code, virtual) => {
    await command('Input.dispatchKeyEvent', {type:'keyDown', key, code, windowsVirtualKeyCode:virtual,
      ...(key === 'Enter' ? {text:'\r', unmodifiedText:'\r'} : {})});
    await command('Input.dispatchKeyEvent', {type:'keyUp', key, code, windowsVirtualKeyCode:virtual});
    await settle();
  };

  for (const [route, family, shell, toggle] of [
    ['/templates/collapsible-sidebar.html', '.lc-cs-nav-item', '.lc-cs-shell', '[data-sidebar-toggle]'],
    ['/workbenches/template', '.lc-wb-nav-item', '.lc-wb-shell', '[data-lc-wb-sidebar-toggle]'],
  ]) {
    console.log('Checking navigation state transitions: ' + route);
    await navigate(route, family);
    const first = family + '[aria-label="Cable geometry"]';
    const second = family + '[aria-label="System overview"]';
    const current = '[aria-current="page"]';
    // Wait for Bonito handlers as well as markup; repeating the same selection
    // is idempotent and never dispatches a numerical job in this fixture.
    await wait(`(() => {document.querySelector(${q(second)}).click(); return ${selected(second, current)};})()`,
      'navigation behavior did not initialize');
    for (const value of ['light', 'dark', 'light']) {
      await theme(value);
      for (const collapsed of [false, true]) {
        const state = collapsed ? 'collapsed' : 'expanded';
        if (!await evaluate(`document.querySelector(${q(shell)}).dataset.sidebarState === ${q(state)}`)) {
          await click(toggle);
          await wait(`document.querySelector(${q(shell)}).dataset.sidebarState === ${q(state)}`, 'collapse did not update');
          await settle();
        }
        await away();
        const normal = await style(first);
        assert.equal(normal.bg, 'rgba(0, 0, 0, 0)', `${route}: unselected navigation is painted`);
        await pointAt(first);
        const hoverState = await evaluate(`(() => {
          const el = document.querySelector(${q(first)}), r = el.getBoundingClientRect();
          return {hovered: el.matches(':hover'), rect: r.toJSON(),
            hit: document.elementFromPoint(r.x + r.width / 2, r.y + r.height / 2)?.outerHTML.slice(0, 220),
            background: getComputedStyle(el).backgroundColor};
        })()`);
        assert.notEqual(hoverState.background, normal.bg,
          `${route}: genuine hover is missing: ${JSON.stringify(hoverState)}`);
        if (collapsed) {
          await wait(`getComputedStyle(document.querySelector(${q(first)}), '::after').opacity === '1'`, 'rail tooltip did not appear');
          assert.notEqual(await evaluate(`getComputedStyle(document.querySelector(${q(first)}), '::after').display`), 'none');
        }
        await click(first);
        await expectSelected(first, current);
        await click(second);
        await expectSelected(second, current);
        await away();
        assert.deepEqual(await style(first), normal, `${route} ${value} ${state}: stale highlight after selection changed`);
        assert.equal(await evaluate(`document.querySelectorAll(${q(family + current)}).length`), 1, 'multiple active items');
        // The marker and tooltip cannot continue to claim the old selection.
        if (family === '.lc-cs-nav-item') {
          if (!collapsed) {
            assert.equal(await evaluate(`getComputedStyle(document.querySelector(${q(first)}), '::after').content`), 'none');
            assert.equal(await evaluate(`getComputedStyle(document.querySelector(${q(second)}), '::after').content`), '""');
          } else {
            await wait(`getComputedStyle(document.querySelector(${q(first)}), '::after').opacity === '0'`, 'old rail tooltip persisted');
          }
          assert(!await evaluate(`Boolean(document.querySelector(${q(shell)}).querySelector('[data-tooltip*="active"], [data-tooltip*="hovered"]'))`));
        }
        const disabled = family + ':disabled';
        const disabledBefore = await style(disabled);
        await click(disabled);
        assert.deepEqual(await style(disabled), disabledBefore, 'disabled navigation responds visually to hover/click');
        await expectSelected(second, current);
        await away();
        // System overview precedes Cable geometry in both shells. A real Tab
        // must give geometry an outline without silently selecting it.
        await click(second);
        await key('Tab', 'Tab', 9);
        await away();
        assert(await evaluate(`document.querySelector(${q(first)}).matches(':focus-visible')`), 'keyboard focus cue missing');
        assert.notEqual(await evaluate(`getComputedStyle(document.querySelector(${q(first)})).outlineStyle`), 'none');
        assert.deepEqual(await style(first), normal, 'keyboard focus was painted as persistent selection');
        await key('Enter', 'Enter', 13);
        await expectSelected(first, current);
        await click(second);
        await expectSelected(second, current);
        await away();
        report.push({route, theme:value, state});
      }
    }
  }

  // Selection semantics differ legitimately (tabs, checked radios, data rows).
  // Each must restore the previous item's exact normal styling, not just ARIA.
  const groups = [
    ['/widgets/ribbon', '.lc-ribbon-tab', '[aria-selected="true"]'],
    ['/fixture/standalone', '.lc-ribbon-tab', '[aria-selected="true"]'],
    ['/fixture/workbench', '.lc-ribbon-tab', '[aria-selected="true"]'],
    ['/widgets/tabs', '.bw-tab', '.bw-active'],
    ['/workbenches/template', '.lc-wb-dock-tab', '.is-active'],
    ['/widgets/form-toolkit', '.lc-segmented-control-choice', ':has(input:checked)'],
    ['/widgets/data-view-toolkit', '.lc-data-row', '[aria-selected="true"]'],
  ];
  for (const [route, family, marker] of groups) {
    console.log('Checking selection round trip: ' + route + ' ' + family);
    await navigate(route, family);
    await evaluate(`window.stateItems = [...document.querySelectorAll(${q(family)})].slice(0,2);
      stateItems.forEach((el,i) => el.dataset.stateAudit = String(i));`);
    const first = '[data-state-audit="1"]', second = '[data-state-audit="0"]';
    await wait(`(() => {document.querySelector(${q(second)}).click(); return ${selected(second,marker)};})()`, 'selection behavior did not initialize');
    for (const value of ['light', 'dark']) {
      await theme(value);
      await away();
      const normal = await style(first);
      await click(first);
      await expectSelected(first, marker);
      await click(second);
      await expectSelected(second, marker);
      await away();
      assert.deepEqual(await style(first), normal, `${route} ${value}: old selection still highlighted`);
      assert.equal(await evaluate(`document.querySelectorAll(${q(family + marker)}).length`), 1, `${route}: multiple selected items`);
      report.push({route, theme:value, family});
    }
  }
  // Toolbar actions are not selections: pointer focus after clicking an action
  // must not leave its hover paint behind. Configured active buttons stay active.
  for (const route of ['/widgets/toolbar', '/fixture/standalone', '/fixture/workbench']) {
    await navigate(route, '.lc-toolbar-button');
    const button = '.lc-toolbar-button:not(.lc-toolbar-button-active):not(:disabled)';
    for (const value of ['light', 'dark']) {
      await theme(value);
      await away();
      const normal = await style(button);
      await click(button);
      await away();
      assert.deepEqual(await style(button), normal, `${route} ${value}: action retained hover after click`);
      report.push({route, theme:value, family:'toolbar action'});
    }
  }
  return report;
}
