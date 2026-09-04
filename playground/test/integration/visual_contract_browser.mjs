#!/usr/bin/env node

const baseUrl = process.argv[2] ?? "http://127.0.0.1:8080";
const debuggingUrl = process.argv[3] ?? "http://127.0.0.1:9222";

const routes = [
  ["/widgets/slider", ".lc-widget-shell"],
  ["/widgets/toggle", ".lc-widget-shell"],
  ["/widgets/dropdown", ".lc-native-select"],
  ["/widgets/text-input", ".lc-native-field"],
  ["/widgets/number-spinner", ".lc-native-number"],
  ["/widgets/actions", ".lc-widget-actions"],
  ["/widgets/progress", ".lc-progress"],
  ["/widgets/console", ".lc-console"],
  ["/widgets/tabs", ".lc-tabs-demo"],
  ["/widgets/toolbar", ".lc-toolbar-select"],
  ["/widgets/ribbon", ".lc-ribbon"],
  ["/widgets/control-panel", ".lc-native-select"],
  ["/widgets/form-toolkit", ".lc-form"],
  ["/widgets/overlay-toolkit", ".lc-dialog"],
  ["/widgets/data-view-toolkit", ".lc-viewport-frame"],
  ["/widgets/repeater", ".lc-repeater"],
  ["/widgets/job-panel", ".lc-job-panel"],
  ["/widgets/file-upload", ".lc-upload-field"],
  ["/widgets/geographic-map", ".lc-map-component"],
  ["/widgets/power-system-canvas", ".lc-power-system-canvas"],
];

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
    const result = new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject });
    });
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
    throw new Error(result.exceptionDetails.text ?? "browser evaluation failed");
  }
  return result.result.value;
}

async function navigate(devtools, url, selector) {
  const loaded = devtools.once("Page.loadEventFired");
  await devtools.command("Page.navigate", { url });
  await Promise.race([
    loaded,
    new Promise((_, reject) => setTimeout(
      () => reject(new Error(`page load timed out: ${url}`)),
      15_000,
    )),
  ]);
  const deadline = Date.now() + 15_000;
  while (Date.now() < deadline) {
    const ready = await evaluate(
      devtools,
      `location.href === ${JSON.stringify(url)} &&
        document.readyState === "complete" &&
        Boolean(document.querySelector(${JSON.stringify(selector)}))`,
    );
    if (ready) return;
    await new Promise(resolve => setTimeout(resolve, 100));
  }
  throw new Error(`selector did not mount: ${selector} at ${url}`);
}

async function waitUntil(devtools, expression, message, timeout = 5_000) {
  const deadline = Date.now() + timeout;
  while (Date.now() < deadline) {
    if (await evaluate(devtools, expression)) return;
    await new Promise(resolve => setTimeout(resolve, 50));
  }
  throw new Error(message);
}

async function clickButton(devtools, label, scope = "document") {
  const clicked = await evaluate(devtools, `(() => {
    const root = ${scope};
    const button = Array.from(root.querySelectorAll('button')).find(
      candidate => candidate.textContent.trim() === ${JSON.stringify(label)}
    );
    if (!button) return false;
    button.click();
    return true;
  })()`);
  assert(clicked, `button not found: ${label}`);
}

async function hoverSelector(devtools, selector) {
  const point = await evaluate(devtools, `(() => {
    const element = document.querySelector(${JSON.stringify(selector)});
    if (!element) return null;
    element.scrollIntoView({block: "center", inline: "center"});
    const box = element.getBoundingClientRect();
    return {x: box.left + box.width / 2, y: box.top + box.height / 2};
  })()`);
  assert(point, `hover target not found: ${selector}`);
  await devtools.command("Input.dispatchMouseEvent", {
    type: "mouseMoved",
    x: point.x,
    y: point.y,
  });
  await new Promise(resolve => setTimeout(resolve, 75));
}

async function setTheme(devtools, theme) {
  await navigate(devtools, `${baseUrl}/`, "#quarto-sidebar");
  await evaluate(
    devtools,
    `localStorage.setItem("lcm.playground.theme", ${JSON.stringify(theme)})`,
  );
}

async function inspectPalette(devtools) {
  return evaluate(devtools, `(() => {
    const style = getComputedStyle(document.documentElement);
    const page = document.querySelector(".lc-wb-page");
    return {
      theme: page?.dataset.resolvedTheme ??
        document.documentElement.dataset.lcmResolvedTheme,
      text: style.getPropertyValue("--lc-text").trim(),
      heading: style.getPropertyValue("--lc-heading").trim(),
      option: style.getPropertyValue("--lc-option-bg").trim(),
      selected: style.getPropertyValue("--lc-option-selected-bg").trim(),
      horizontalOverflow: document.documentElement.scrollWidth > window.innerWidth + 1,
    };
  })()`);
}

async function waitForPalette(devtools, theme) {
  const deadline = Date.now() + 5_000;
  let palette = await inspectPalette(devtools);
  while (palette.theme !== theme && Date.now() < deadline) {
    await new Promise(resolve => setTimeout(resolve, 50));
    palette = await inspectPalette(devtools);
  }
  return palette;
}

async function inspectSelect(devtools, selector) {
  return evaluate(devtools, `(() => {
    const select = document.querySelector(${JSON.stringify(selector)});
    if (!select) return null;
    const options = Array.from(select.options);
    const selected = options.find(option => option.selected) ?? options[0];
    const ordinary = options.find(option => !option.selected) ?? selected;
    const selectStyle = getComputedStyle(select);
    const selectedStyle = getComputedStyle(selected);
    const ordinaryStyle = getComputedStyle(ordinary);
    return {
      contract: select.classList.contains("lc-control-select"),
      scheme: selectStyle.colorScheme,
      selectedColor: selectedStyle.color,
      selectedBackground: selectedStyle.backgroundColor,
      ordinaryColor: ordinaryStyle.color,
      ordinaryBackground: ordinaryStyle.backgroundColor,
    };
  })()`);
}

async function waitForSelectTheme(devtools, selector, theme) {
  const deadline = Date.now() + 5_000;
  let sample = await inspectSelect(devtools, selector);
  while (sample?.scheme !== theme && Date.now() < deadline) {
    await new Promise(resolve => setTimeout(resolve, 50));
    sample = await inspectSelect(devtools, selector);
  }
  return sample;
}

async function inspectShell(devtools) {
  return evaluate(devtools, `(() => {
    const sidebar = document.querySelector("#quarto-sidebar");
    const footer = document.querySelector("footer.footer");
    const sidebarBox = sidebar.getBoundingClientRect();
    const footerBox = footer.getBoundingClientRect();
    const sidebarStyle = getComputedStyle(sidebar);
    const footerStyle = getComputedStyle(footer);
    return {
      location: location.href,
      bodyClass: document.body.className,
      sidebarClass: sidebar.className,
      sidebarPosition: sidebarStyle.position,
      sidebarDisplay: sidebarStyle.display,
      viewportHeight: innerHeight,
      sidebarTop: sidebarBox.top,
      sidebarBottom: sidebarBox.bottom,
      sidebarRight: sidebarBox.right,
      footerBottom: footerBox.bottom,
      footerRight: footerBox.right,
      sidebarBorderRight: sidebarStyle.borderRightWidth,
      footerBorderRight: footerStyle.borderRightWidth,
      footerBorderTop: footerStyle.borderTopWidth,
    };
  })()`);
}

async function inspectCodeContract(devtools) {
  return evaluate(devtools, `(() => {
    const code = document.querySelector('pre.sourceCode code.sourceCode:has(span.kw)');
    const surface = code?.closest('div.sourceCode') ?? code?.parentElement;
    if (!code || !surface) return null;
    const parse = value => (value.match(/[\\d.]+/g) ?? []).slice(0, 3).map(Number);
    const luminance = value => {
      const channels = parse(value).map(channel => {
        const normalized = channel / 255;
        return normalized <= 0.04045
          ? normalized / 12.92
          : ((normalized + 0.055) / 1.055) ** 2.4;
      });
      return 0.2126 * channels[0] + 0.7152 * channels[1] + 0.0722 * channels[2];
    };
    const contrast = (foreground, background) => {
      const light = Math.max(luminance(foreground), luminance(background));
      const dark = Math.min(luminance(foreground), luminance(background));
      return (light + 0.05) / (dark + 0.05);
    };
    const background = getComputedStyle(surface).backgroundColor;
    const roles = Array.from(code.querySelectorAll('span[class]'))
      .filter(element => element.textContent.trim().length > 0)
      .map(element => ({
        role: element.className,
        color: getComputedStyle(element).color,
      }))
      .filter((sample, index, all) =>
        all.findIndex(candidate => candidate.role === sample.role) === index)
      .map(sample => ({
        ...sample,
        contrast: contrast(sample.color, background),
      }));
    const base = getComputedStyle(code).color;
    return {background, base, baseContrast: contrast(base, background), roles};
  })()`);
}

async function inspectSidebarSpecimen(devtools) {
  return evaluate(devtools, `(() => {
    const resolveColor = variable => {
      const probe = document.createElement('i');
      probe.style.color = 'var(' + variable + ')';
      document.body.append(probe);
      const value = getComputedStyle(probe).color;
      probe.remove();
      return value;
    };
    const resolveBackground = variable => {
      const probe = document.createElement('i');
      probe.style.backgroundColor = 'var(' + variable + ')';
      document.body.append(probe);
      const value = getComputedStyle(probe).backgroundColor;
      probe.remove();
      return value;
    };
    const hovered = getComputedStyle(document.querySelector('.lc-cs-nav-item.is-hover-preview'));
    const active = getComputedStyle(document.querySelector('.lc-cs-nav-item.is-active'));
    const scene = getComputedStyle(document.querySelector('.lc-cs-figure > svg'));
    return {
      hovered: {color: hovered.color, background: hovered.backgroundColor},
      active: {color: active.color, background: active.backgroundColor},
      scene: {background: scene.backgroundColor},
      tokens: {
        strong: resolveColor('--lc-strong-text'),
        hover: resolveBackground('--lc-hover-bg'),
        active: resolveBackground('--lc-active-bg'),
        scene: resolveBackground('--lc-scene-bg'),
      },
    };
  })()`);
}

async function inspectSourceButton(devtools) {
  return evaluate(devtools, `(() => {
    const resolveColor = variable => {
      const probe = document.createElement('i');
      probe.style.color = 'var(' + variable + ')';
      document.body.append(probe);
      const value = getComputedStyle(probe).color;
      probe.remove();
      return value;
    };
    const resolveBackground = variable => {
      const probe = document.createElement('i');
      probe.style.backgroundColor = 'var(' + variable + ')';
      document.body.append(probe);
      const value = getComputedStyle(probe).backgroundColor;
      probe.remove();
      return value;
    };
    const button = document.querySelector('#quarto-code-tools-source');
    if (!button) return null;
    const style = getComputedStyle(button);
    return {
      hovered: button.matches(':hover'),
      color: style.color,
      background: style.backgroundColor,
      tokens: {
        strong: resolveColor('--lc-strong-text'),
        hover: resolveBackground('--lc-hover-bg'),
      },
    };
  })()`);
}

async function inspectWorkbenchAffordances(devtools) {
  return evaluate(devtools, `(() => {
    const split = document.querySelector('.lc-wb-split-horizontal');
    const handle = split?.querySelector(':scope > .lc-wb-splitter .lc-wb-splitter-handle');
    const first = split?.querySelector(':scope > .lc-wb-split-first');
    const inspector = document.querySelector('.lc-wb-inspector');
    const xrayHost = document.querySelector('.lc-xray-host');
    const xrayToggle = xrayHost?.shadowRoot?.querySelector('.xray-toggle');
    const splitBox = split?.getBoundingClientRect();
    const handleBox = handle?.getBoundingClientRect();
    const inspectorBox = inspector?.getBoundingClientRect();
    const xrayBox = xrayToggle?.getBoundingClientRect();
    return {
      split: splitBox ? {left: splitBox.left, top: splitBox.top} : null,
      handle: handleBox ? {
        left: handleBox.left,
        top: handleBox.top,
        width: handleBox.width,
        height: handleBox.height,
      } : null,
      firstWidth: first?.getBoundingClientRect().width ?? 0,
      legacyRule: split
        ? getComputedStyle(split.querySelector(':scope > .lc-wb-splitter'), '::after').content
        : null,
      inspectorInteractiveChildren: inspector
        ? inspector.querySelectorAll('.lc-wb-inspector-header button, .lc-wb-inspector-header [role="button"]').length
        : 0,
      inspectorLeft: inspectorBox?.left ?? null,
      xray: xrayBox ? {
        top: xrayBox.top,
        right: xrayBox.right,
      } : null,
    };
  })()`);
}

async function inspectXRayGrid(devtools) {
  return evaluate(devtools, `(() => {
    const shadow = document.querySelector('.lc-xray-host')?.shadowRoot;
    const grids = Array.from(shadow?.querySelectorAll('.xray-grid') ?? []);
    const samples = grids.map(grid => {
      const columns = grid.classList.contains('xray-grid-4') ? 4 : 2;
      const cells = Array.from(grid.querySelectorAll(':scope > .xray-cell'));
      const rows = [];
      for (let index = 0; index < cells.length; index += columns) {
        const boxes = cells.slice(index, index + columns).map(cell => cell.getBoundingClientRect());
        rows.push({
          topSpread: Math.max(...boxes.map(box => box.top)) - Math.min(...boxes.map(box => box.top)),
          bottomSpread: Math.max(...boxes.map(box => box.bottom)) - Math.min(...boxes.map(box => box.bottom)),
        });
      }
      return {
        columns,
        firstDivider: cells[0]?.getBoundingClientRect().right ?? null,
        rows,
      };
    });
    return {
      samples,
      maximumRowDrift: Math.max(0, ...samples.flatMap(sample =>
        sample.rows.flatMap(row => [row.topSpread, row.bottomSpread]))),
      firstDividerDrift: Math.max(0, ...samples.map(sample => sample.firstDivider)) -
        Math.min(...samples.map(sample => sample.firstDivider)),
    };
  })()`);
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
await devtools.command("Emulation.setDeviceMetricsOverride", {
  width: 1600,
  height: 1000,
  deviceScaleFactor: 1,
  mobile: false,
});

try {
  const report = {};
  for (const theme of ["dark", "light"]) {
    await setTheme(devtools, theme);
    await navigate(devtools, `${baseUrl}/`, "#quarto-sidebar");
    const baseline = await waitForPalette(devtools, theme);
    assert(baseline.theme === theme, `home did not resolve ${theme} theme`);
    const shell = await inspectShell(devtools);
    assert(Math.abs(shell.sidebarTop) <= 1, `${theme} sidebar does not start at viewport top`);
    assert(
      Math.abs(shell.sidebarBottom - shell.viewportHeight) <= 1,
      `${theme} sidebar does not reach viewport bottom`,
    );
    assert(
      Math.abs(shell.footerBottom - shell.viewportHeight) <= 1,
      `${theme} footer does not meet viewport bottom`,
    );
    assert(
      Math.abs(shell.sidebarRight - shell.footerRight) <= 1,
      `${theme} sidebar/footer right contours are misaligned`,
    );
    assert(shell.sidebarBorderRight === "1px", `${theme} sidebar contour is missing`);
    assert(shell.footerBorderRight === "1px", `${theme} footer contour is missing`);
    assert(shell.footerBorderTop === "1px", `${theme} footer top contour is missing`);

    const specimens = {};
    for (const [route, selector] of routes) {
      await navigate(devtools, `${baseUrl}${route}`, selector);
      const palette = await waitForPalette(devtools, theme);
      assert(palette.theme === theme, `${route} did not resolve ${theme} theme`);
      assert(palette.text === baseline.text, `${route} text token drifted in ${theme}`);
      assert(palette.heading === baseline.heading, `${route} heading token drifted in ${theme}`);
      assert(palette.option === baseline.option, `${route} option token drifted in ${theme}`);
      assert(!palette.horizontalOverflow, `${route} overflows horizontally in ${theme}`);
      specimens[route] = palette;
    }

    const selectors = [
      ["/", ".lc-publisher-theme-select"],
      ["/widgets/dropdown", ".lc-native-select"],
      ["/widgets/toolbar", ".lc-toolbar-select"],
      ["/widgets/ribbon", ".lc-toolbar-select"],
      ["/widgets/control-panel", ".lc-native-select"],
      ["/widgets/form-toolkit", "#earth_model"],
      ["/widgets/form-toolkit", "#outputs"],
      ["/widgets/repeater", ".lc-native-select"],
      ["/widgets/geographic-map", ".lc-map-field select"],
      ["/workbenches/template", ".lc-wb-choice-select"],
      ["/workbenches/template", ".lc-wb-theme-select"],
    ];
    const selectReport = [];
    for (const [route, selector] of selectors) {
      await navigate(devtools, `${baseUrl}${route}`, selector);
      const sample = await waitForSelectTheme(devtools, selector, theme);
      assert(sample?.contract, `${route} ${selector} is outside the select contract`);
      assert(sample.scheme === theme, `${route} ${selector} uses the wrong color scheme`);
      const values = [
        sample.selectedColor,
        sample.selectedBackground,
        sample.ordinaryColor,
        sample.ordinaryBackground,
      ];
      assert(
        values.every(value => value && value !== "rgba(0, 0, 0, 0)"),
        `${route} ${selector} has an unresolved foreground/background in ${theme}`,
      );
      selectReport.push({ route, selector, ...sample });
    }
    const reference = selectReport[0];
    for (const sample of selectReport.slice(1)) {
      for (const property of [
        "scheme",
        "selectedColor",
        "selectedBackground",
        "ordinaryColor",
        "ordinaryBackground",
      ]) {
        assert(
          sample[property] === reference[property],
          `${sample.route} ${sample.selector} drifted from the shared ${property} contract in ${theme}`,
        );
      }
    }

    await navigate(devtools, `${baseUrl}/workbenches/template`, ".lc-wb-splitter-handle");
    const workbenchBefore = await inspectWorkbenchAffordances(devtools);
    assert(workbenchBefore.handle?.width >= 18 && workbenchBefore.handle?.height >= 18,
      `split handle is not a usable target in ${theme}`);
    assert(workbenchBefore.handle.top - workbenchBefore.split.top <= 16,
      `split handle is not anchored near the leading edge in ${theme}`);
    assert(["none", "normal"].includes(workbenchBefore.legacyRule),
      `legacy full-height split rule is still rendered in ${theme}`);
    assert(workbenchBefore.inspectorInteractiveChildren === 0,
      `inspector header exposes an inert control in ${theme}`);
    if (workbenchBefore.xray) {
      assert(workbenchBefore.xray.top <= 8,
        `X-RAY launcher is not in the top chrome in ${theme}`);
      assert(workbenchBefore.xray.right <= workbenchBefore.inspectorLeft,
        `X-RAY launcher overlaps the inspector in ${theme}`);
    }
    await evaluate(devtools, `globalThis.lcmXRay?.enable()`);
    await hoverSelector(devtools, ".lc-wb-menubar");
    await waitUntil(devtools,
      `Boolean(document.querySelector('.lc-xray-host')?.shadowRoot?.querySelector('.xray-panel:not([hidden]) .xray-grid-4'))`,
      `X-RAY metadata grid did not render in ${theme}`);
    const xrayGrid = await inspectXRayGrid(devtools);
    assert(xrayGrid.maximumRowDrift <= 0.1,
      `X-RAY cell borders drift by ${xrayGrid.maximumRowDrift}px in ${theme}`);
    assert(xrayGrid.firstDividerDrift <= 0.5,
      `X-RAY key column drifts by ${xrayGrid.firstDividerDrift}px between sections in ${theme}`);
    await evaluate(devtools, `globalThis.lcmXRay?.disable()`);
    const handlePoint = {
      x: workbenchBefore.handle.left + workbenchBefore.handle.width / 2,
      y: workbenchBefore.handle.top + workbenchBefore.handle.height / 2,
    };
    await devtools.command("Input.dispatchMouseEvent", {
      type: "mousePressed",
      x: handlePoint.x,
      y: handlePoint.y,
      button: "left",
      buttons: 1,
      clickCount: 1,
    });
    await devtools.command("Input.dispatchMouseEvent", {
      type: "mouseMoved",
      x: handlePoint.x + 48,
      y: handlePoint.y,
      button: "left",
      buttons: 1,
    });
    await devtools.command("Input.dispatchMouseEvent", {
      type: "mouseReleased",
      x: handlePoint.x + 48,
      y: handlePoint.y,
      button: "left",
      buttons: 0,
      clickCount: 1,
    });
    const workbenchAfter = await inspectWorkbenchAffordances(devtools);
    assert(workbenchAfter.firstWidth - workbenchBefore.firstWidth >= 30,
      `split handle did not resize its pane in ${theme}`);

    await navigate(devtools, `${baseUrl}/workbenches/`,
      "pre.sourceCode code.sourceCode span.kw");
    const codeContract = await inspectCodeContract(devtools);
    assert(codeContract, `code contract did not mount in ${theme}`);
    assert(codeContract.baseContrast >= 4.5,
      `base code text fails contrast in ${theme}: ${codeContract.baseContrast}`);
    assert(codeContract.roles.length >= 4,
      `syntax-role coverage is unexpectedly narrow in ${theme}`);
    for (const role of codeContract.roles) {
      assert(role.contrast >= 4.5,
        `syntax role ${role.role} fails contrast in ${theme}: ${role.contrast}`);
    }

    await navigate(devtools, `${baseUrl}/templates/collapsible-sidebar.html`,
      ".lc-cs-shell");
    const sidebarSpecimen = await inspectSidebarSpecimen(devtools);
    assert(sidebarSpecimen.hovered.color === sidebarSpecimen.tokens.strong,
      `sidebar hover foreground escaped the ${theme} palette`);
    assert(sidebarSpecimen.hovered.background === sidebarSpecimen.tokens.hover,
      `sidebar hover background escaped the ${theme} palette`);
    assert(sidebarSpecimen.active.color === sidebarSpecimen.tokens.strong,
      `sidebar active foreground escaped the ${theme} palette`);
    assert(sidebarSpecimen.active.background === sidebarSpecimen.tokens.active,
      `sidebar active background escaped the ${theme} palette`);
    assert(sidebarSpecimen.scene.background === sidebarSpecimen.tokens.scene,
      `sidebar scene escaped the shared engineering-scene palette in ${theme}`);

    await devtools.command("Input.dispatchMouseEvent", {
      type: "mouseMoved", x: 2, y: 998,
    });
    await hoverSelector(devtools, "#quarto-code-tools-source");
    const sourceButton = await inspectSourceButton(devtools);
    assert(sourceButton?.hovered, `Source button did not enter hover state in ${theme}`);
    assert(sourceButton.color === sourceButton.tokens.strong,
      `Source hover foreground escaped the ${theme} palette`);
    assert(sourceButton.background === sourceButton.tokens.hover,
      `Source hover background escaped the ${theme} palette`);

    await navigate(devtools, `${baseUrl}/widgets/form-toolkit`, ".lc-form");
    const secretBoundary = await evaluate(devtools, `(() => {
      const secret = document.querySelector('#credential');
      return {
        type: secret.type,
        serializedValue: secret.getAttribute('value'),
        reactiveAttribute: secret.hasAttribute('oninput') || secret.hasAttribute('onchange'),
      };
    })()`);
    assert(secretBoundary.type === "password", "SecretInput is not masked");
    assert(secretBoundary.serializedValue === null, "SecretInput serialized a value");
    assert(!secretBoundary.reactiveAttribute, "SecretInput exposes an inline reactive binding");
    await evaluate(devtools, `(() => {
      const field = document.querySelector('#scenario');
      field.value = 'Changed';
      field.dispatchEvent(new Event('input', {bubbles: true}));
      return true;
    })()`);
    await clickButton(devtools, "Restore defaults");
    await waitUntil(devtools,
      `document.querySelector('#scenario').value === 'North corridor'`,
      "form reset did not restore its declared initial value");

    await navigate(devtools, `${baseUrl}/widgets/overlay-toolkit`, ".lc-dialog");
    await new Promise(resolve => setTimeout(resolve, 750));
    await clickButton(devtools, "Message dialog");
    await waitUntil(devtools,
      `document.querySelector('#message_specimen').open`,
      "message dialog did not open");
    await evaluate(devtools,
      `document.querySelector('#message_specimen .lc-dialog-close').click()`);
    await waitUntil(devtools,
      `!document.querySelector('#message_specimen').open`,
      "message dialog did not close");
    await clickButton(devtools, "Confirmation dialog");
    await waitUntil(devtools,
      `document.querySelector('#confirmation_specimen').open`,
      "confirmation dialog did not open");
    await evaluate(devtools,
      `document.querySelector('#confirmation_specimen .lc-dialog-action.is-danger').click()`);
    await waitUntil(devtools,
      `document.querySelector('.lc-widget-value').textContent.includes('confirmation_specimen.confirm')`,
      "confirmation action did not reach the typed event boundary");
    await clickButton(devtools, "Push toast");
    await waitUntil(devtools,
      `document.querySelectorAll('.lc-toast').length === 1`,
      "toast did not mount");

    await navigate(devtools, `${baseUrl}/widgets/data-view-toolkit`, ".lc-viewport-frame");
    await evaluate(devtools,
      `window.__lcmPersistentTable = document.querySelector('.lc-data-table')`);
    await clickButton(devtools, "Loading");
    await waitUntil(devtools,
      `document.querySelector('.lc-viewport-frame').classList.contains('is-loading')`,
      "viewport loading state did not render");
    const persistent = await evaluate(devtools,
      `window.__lcmPersistentTable === document.querySelector('.lc-data-table')`);
    assert(persistent, "viewport state transition remounted persistent content");
    await evaluate(devtools, `(() => {
      const filter = document.querySelector('.lc-data-filter');
      filter.value = 'Resistance';
      filter.dispatchEvent(new Event('input', {bubbles: true}));
      return true;
    })()`);
    try {
      await waitUntil(devtools,
        `document.querySelectorAll('.lc-data-table tbody tr').length === 1`,
        "data-table filter did not update rendered rows");
    } catch (error) {
      const diagnostics = await evaluate(devtools, `(() => ({
        value: document.querySelector('.lc-data-filter')?.value,
        rows: Array.from(document.querySelectorAll('.lc-data-table tbody tr')).map(row => row.innerText),
        body: document.querySelector('.lc-data-table tbody')?.innerHTML,
      }))()`);
      throw new Error(`${error.message}: ${JSON.stringify(diagnostics)}`);
    }

    await navigate(devtools, `${baseUrl}/widgets/`,
      `iframe[title="Typed field and form toolkit"]`);
    const galleryFrames = [
      ["Typed field and form toolkit", ".lc-form"],
      ["Overlay and feedback toolkit", ".lc-dialog"],
      ["Data-view toolkit", ".lc-viewport-frame"],
    ];
    const gallerySizing = [];
    for (const [title, selector] of galleryFrames) {
      await evaluate(devtools, `document.querySelector(
        'iframe[title=${JSON.stringify(title)}]'
      ).scrollIntoView({block: 'center'})`);
      await waitUntil(devtools, `(() => {
        const frame = document.querySelector('iframe[title=${JSON.stringify(title)}]');
        return Boolean(frame?.contentDocument?.querySelector(${JSON.stringify(selector)}));
      })()`, `gallery iframe did not load: ${title}`, 10_000);
      const sizing = await evaluate(devtools, `(() => {
        const frame = document.querySelector('iframe[title=${JSON.stringify(title)}]');
        const root = frame.contentDocument.documentElement;
        return {
          title: ${JSON.stringify(title)},
          clientHeight: frame.clientHeight,
          scrollHeight: root.scrollHeight,
          horizontalOverflow: root.scrollWidth > frame.clientWidth + 1,
          verticalOverflow: root.scrollHeight > frame.clientHeight + 1,
        };
      })()`);
      assert(!sizing.horizontalOverflow, `${title} gallery viewport overflows horizontally`);
      assert(!sizing.verticalOverflow,
        `${title} gallery viewport has an internal scrollbar (${sizing.scrollHeight}px content in ${sizing.clientHeight}px)`);
      gallerySizing.push(sizing);
    }

    report[theme] = {
      baseline,
      shell,
      specimens,
      selects: selectReport,
      codeContract,
      sidebarSpecimen,
      sourceButton,
      workbenchAffordances: {before: workbenchBefore, after: workbenchAfter},
      xrayGrid,
      toolkit: { secretBoundary, persistentViewport: persistent },
      gallerySizing,
    };
  }
  console.log(JSON.stringify(report, null, 2));
} finally {
  devtools.close();
}
