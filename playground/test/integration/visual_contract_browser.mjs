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
          horizontalOverflow: root.scrollWidth > frame.clientWidth + 1,
          verticalOverflow: root.scrollHeight > frame.clientHeight + 1,
        };
      })()`);
      assert(!sizing.horizontalOverflow, `${title} gallery viewport overflows horizontally`);
      assert(!sizing.verticalOverflow, `${title} gallery viewport has an internal scrollbar`);
      gallerySizing.push(sizing);
    }

    report[theme] = {
      baseline,
      shell,
      specimens,
      selects: selectReport,
      toolkit: { secretBoundary, persistentViewport: persistent },
      gallerySizing,
    };
  }
  console.log(JSON.stringify(report, null, 2));
} finally {
  devtools.close();
}
