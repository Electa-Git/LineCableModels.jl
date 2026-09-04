(() => {
const controllers = globalThis.__lcmXRayControllers ||= new Set();
const pending = globalThis.__lcmXRayPending ||= [];

function formatValue(value) {
  if (value === null) return "nothing";
  if (value === undefined) return "undefined";
  if (typeof value === "string") return value;
  if (["number", "boolean", "bigint"].includes(typeof value)) return String(value);
  if (Array.isArray(value)) {
    const preview = value.slice(0, 6).map(formatValue).join(", ");
    return `[${preview}${value.length > 6 ? ", …" : ""}]`;
  }
  try {
    const serialized = JSON.stringify(value);
    if (serialized && serialized.length <= 120) return serialized;
  } catch (_) {
    // Fall through to the object's compact string representation.
  }
  return String(value);
}

function dispatch(entry) {
  let claimed = false;
  for (const controller of controllers) {
    if (!controller.owns(entry.node)) continue;
    controller.accept(entry);
    claimed = true;
  }
  if (!claimed) pending.push(entry);
}

globalThis.__lcmXRayDispatch = dispatch;

function element(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

function appendRows(parent, rows, columns) {
  const grid = element("div", `xray-grid xray-grid-${columns}`);
  for (const row of rows) {
    for (const value of row) grid.append(element("span", "xray-cell", value));
  }
  parent.append(grid);
}

function visitRules(rules, scopes, result, seen) {
  for (const rule of rules) {
    if (rule.selectorText && rule.style) {
      const owned = scopes.some(scope => rule.selectorText.includes(scope));
      if (!owned) continue;
      const declarations = [];
      for (const property of rule.style) {
        const priority = rule.style.getPropertyPriority(property);
        declarations.push({
          property,
          value: rule.style.getPropertyValue(property).trim(),
          priority,
        });
      }
      const key = `${rule.selectorText}\n${JSON.stringify(declarations)}`;
      if (!seen.has(key)) {
        seen.add(key);
        result.push({ selector: rule.selectorText, declarations });
      }
    } else if (rule.cssRules) {
      visitRules(rule.cssRules, scopes, result, seen);
    }
  }
}

function ownedCss(node, scopes) {
  const result = [];
  const seen = new Set();
  for (const sheet of document.styleSheets) {
    try {
      visitRules(sheet.cssRules || [], scopes, result, seen);
    } catch (_) {
      // Cross-origin stylesheets are intentionally outside the owned contract.
    }
  }
  if (node.getAttribute("style")) {
    result.unshift({
      selector: "element.style",
      declarations: [{ property: "style", value: node.getAttribute("style"), priority: "" }],
    });
  }
  return result;
}

function isTruthyFlag(value) {
  return ["1", "true", "yes", "on"].includes(String(value || "").toLowerCase());
}

function installXRay(root, options = {}, styles = "") {
  for (const existing of controllers) {
    if (existing.root === root) existing.destroy();
  }

  const page = root.closest(".lc-wb-page") || document.body;
  const host = document.createElement("div");
  host.className = "lc-xray-host";
  host.setAttribute("aria-live", "polite");
  page.append(host);
  const shadow = host.attachShadow({ mode: "open" });
  shadow.append(element("style", "", styles));

  const outline = element("div", "xray-outline");
  const badge = element("div", "xray-badge");
  outline.append(badge);

  const toggle = element("button", "xray-toggle", "X-RAY");
  toggle.type = "button";
  toggle.title = `${options.shortcut || "Ctrl+Shift+X"} · component diagnostics`;
  toggle.setAttribute("aria-pressed", "false");

  const panel = element("aside", "xray-panel");
  panel.hidden = true;
  panel.setAttribute("aria-label", "Owned component diagnostics");
  const panelHeader = element("header", "xray-panel-header");
  panelHeader.title = "Drag diagnostics window";
  const heading = element("div", "xray-heading");
  const kicker = element("span", "xray-kicker", "OWNED COMPONENT");
  const title = element("strong", "xray-title");
  heading.append(kicker, title);
  const close = element("button", "xray-close", "×");
  close.type = "button";
  close.setAttribute("aria-label", "Close component diagnostics");
  panelHeader.append(heading, close);
  const breadcrumbs = element("nav", "xray-breadcrumbs");
  breadcrumbs.setAttribute("aria-label", "Owned component ancestors");
  const body = element("div", "xray-body");
  const resizeHandle = element("div", "xray-resize-handle");
  resizeHandle.title = "Resize diagnostics window";
  resizeHandle.setAttribute("aria-hidden", "true");
  panel.append(panelHeader, breadcrumbs, body, resizeHandle);
  shadow.append(outline, toggle, panel);

  const registry = new WeakMap();
  let active = Boolean(options.enabled);
  let hovered = null;
  let selected = null;
  let pinned = false;
  let hoverTimer = 0;
  let frame = 0;
  let resizeObserver = null;
  let panelPositioned = false;
  let panelInteraction = null;
  const abort = new AbortController();
  const signal = abort.signal;

  const queryFlag = new URLSearchParams(window.location.search).get("xray");
  if (queryFlag !== null) active = isTruthyFlag(queryFlag);

  function owns(node) {
    return node instanceof Element && (node === root || root.contains(node));
  }

  function accept(entry) {
    if (entry.kind === "register") {
      registry.set(entry.node, entry.metadata);
      entry.node.setAttribute("data-lcm-inspection-id", entry.metadata.id);
      entry.node.setAttribute("data-lcm-component", entry.metadata.name);
      return;
    }
    if (entry.kind === "binding") {
      const metadata = registry.get(entry.node);
      if (!metadata) return;
      const binding = metadata.bindings.find(candidate => candidate.name === entry.name);
      if (binding) binding.value = formatValue(entry.value);
      if (entry.node === selected) render(entry.node);
    }
  }

  function inspectableFromEvent(event) {
    for (const candidate of event.composedPath()) {
      if (!(candidate instanceof Element)) continue;
      if (!owns(candidate)) continue;
      if (registry.has(candidate)) return candidate;
    }
    return null;
  }

  function ancestors(node) {
    const result = [];
    let candidate = node;
    while (candidate && candidate instanceof Element) {
      if (registry.has(candidate)) result.push(candidate);
      if (candidate === root) break;
      candidate = candidate.parentElement;
    }
    return result.reverse();
  }

  function positionOutline() {
    frame = 0;
    if (!active || !selected || !selected.isConnected) {
      outline.classList.remove("is-visible");
      return;
    }
    const bounds = selected.getBoundingClientRect();
    outline.style.left = `${bounds.left}px`;
    outline.style.top = `${bounds.top}px`;
    outline.style.width = `${bounds.width}px`;
    outline.style.height = `${bounds.height}px`;
    outline.classList.add("is-visible");
    const metadata = registry.get(selected);
    badge.textContent = `${metadata.name}  ${Math.round(bounds.width)} × ${Math.round(bounds.height)}`;
    badge.classList.toggle("is-below", bounds.top < 30);
  }

  function scheduleOutline() {
    if (!frame) frame = requestAnimationFrame(positionOutline);
  }

  function clamp(value, minimum, maximum) {
    return Math.min(Math.max(value, minimum), Math.max(minimum, maximum));
  }

  function materializePanelPosition() {
    const bounds = panel.getBoundingClientRect();
    panel.style.left = `${bounds.left}px`;
    panel.style.top = `${bounds.top}px`;
    panel.style.right = "auto";
    panel.style.width = `${bounds.width}px`;
    panelPositioned = true;
    return bounds;
  }

  function constrainPanel() {
    if (panel.hidden || !panelPositioned) return;
    const margin = 8;
    const bounds = panel.getBoundingClientRect();
    const maximumWidth = Math.max(240, window.innerWidth - (2 * margin));
    const maximumHeight = Math.max(160, window.innerHeight - (2 * margin));
    const width = Math.min(bounds.width, maximumWidth);
    const height = Math.min(bounds.height, maximumHeight);
    const left = clamp(bounds.left, margin, window.innerWidth - width - margin);
    const top = clamp(bounds.top, margin, window.innerHeight - height - margin);
    panel.style.left = `${left}px`;
    panel.style.top = `${top}px`;
    panel.style.width = `${width}px`;
    if (panel.style.height || bounds.height > maximumHeight) {
      panel.style.height = `${height}px`;
    }
  }

  function beginPanelInteraction(event, kind) {
    if (event.button !== 0) return;
    if (kind === "move" && event.target.closest("button")) return;
    event.preventDefault();
    event.stopPropagation();
    const bounds = materializePanelPosition();
    panelInteraction = {
      kind,
      pointerId: event.pointerId,
      startX: event.clientX,
      startY: event.clientY,
      left: bounds.left,
      top: bounds.top,
      width: bounds.width,
      height: bounds.height,
    };
    const handle = kind === "move" ? panelHeader : resizeHandle;
    handle.setPointerCapture(event.pointerId);
    panel.classList.add(kind === "move" ? "is-moving" : "is-resizing");
  }

  function updatePanelInteraction(event) {
    const interaction = panelInteraction;
    if (!interaction || interaction.pointerId !== event.pointerId) return;
    event.preventDefault();
    const margin = 8;
    const deltaX = event.clientX - interaction.startX;
    const deltaY = event.clientY - interaction.startY;
    if (interaction.kind === "move") {
      const left = clamp(
        interaction.left + deltaX,
        margin,
        window.innerWidth - interaction.width - margin,
      );
      const top = clamp(
        interaction.top + deltaY,
        margin,
        window.innerHeight - interaction.height - margin,
      );
      panel.style.left = `${left}px`;
      panel.style.top = `${top}px`;
      return;
    }
    const minimumWidth = Math.min(320, window.innerWidth - interaction.left - margin);
    const minimumHeight = Math.min(180, window.innerHeight - interaction.top - margin);
    const width = clamp(
      interaction.width + deltaX,
      minimumWidth,
      window.innerWidth - interaction.left - margin,
    );
    const height = clamp(
      interaction.height + deltaY,
      minimumHeight,
      window.innerHeight - interaction.top - margin,
    );
    panel.style.width = `${width}px`;
    panel.style.height = `${height}px`;
  }

  function endPanelInteraction(event) {
    const interaction = panelInteraction;
    if (!interaction || interaction.pointerId !== event.pointerId) return;
    const handle = interaction.kind === "move" ? panelHeader : resizeHandle;
    if (handle.hasPointerCapture(event.pointerId)) {
      handle.releasePointerCapture(event.pointerId);
    }
    panel.classList.remove("is-moving", "is-resizing");
    panelInteraction = null;
    constrainPanel();
  }

  function section(label) {
    const sectionNode = element("section", "xray-section");
    sectionNode.append(element("h3", "xray-section-title", label));
    body.append(sectionNode);
    return sectionNode;
  }

  function renderBreadcrumbs(node) {
    breadcrumbs.replaceChildren();
    const chain = ancestors(node);
    chain.forEach((candidate, index) => {
      if (index) breadcrumbs.append(element("span", "xray-separator", "›"));
      const metadata = registry.get(candidate);
      const button = element("button", "xray-crumb", metadata.name);
      button.type = "button";
      button.classList.toggle("is-current", candidate === node);
      button.addEventListener("click", event => {
        event.stopPropagation();
        pinned = true;
        render(candidate);
      });
      breadcrumbs.append(button);
    });
  }

  function render(node) {
    const metadata = registry.get(node);
    if (!metadata) return;
    selected = node;
    title.textContent = metadata.name;
    renderBreadcrumbs(node);
    body.replaceChildren();

    const bounds = node.getBoundingClientRect();
    appendRows(body, [
      ["Julia type", metadata.julia_type],
      ["Rendered size", `${Math.round(bounds.width)} × ${Math.round(bounds.height)} px`],
      ["Source", `${metadata.source.file}:${metadata.source.line}`],
      ["Module", metadata.source.module],
    ], 2);

    if (metadata.parameters.length) {
      const target = section("Parameters");
      appendRows(target, metadata.parameters.map(parameter => [
        parameter.name,
        parameter.value,
        parameter.julia_type,
        parameter.origin,
      ]), 4);
    }

    if (metadata.bindings.length) {
      const target = section("Bindings");
      appendRows(target, metadata.bindings.map(binding => [
        binding.name,
        binding.value,
        binding.julia_type,
        binding.notes || "live",
      ]), 4);
    }

    if (metadata.actions.length) {
      const target = section("Callbacks and actions");
      appendRows(target, metadata.actions.map(action => [
        action.event,
        action.name,
        action.callback_type,
        action.disabled ? "disabled" : action.owner,
      ]), 4);
    }

    const css = ownedCss(node, metadata.css_scopes);
    if (css.length) {
      const target = section("Component-owned CSS");
      for (const rule of css) {
        const block = element("div", "xray-css-rule");
        block.append(element("code", "xray-selector", rule.selector));
        for (const declaration of rule.declarations) {
          const suffix = declaration.priority ? ` !${declaration.priority}` : "";
          block.append(element(
            "code",
            "xray-declaration",
            `${declaration.property}: ${declaration.value}${suffix};`,
          ));
        }
        target.append(block);
      }
    } else if (metadata.css_scopes.length) {
      const target = section("Component-owned CSS");
      target.append(element(
        "p",
        "xray-empty",
        `No authored rules found for ${metadata.css_scopes.join(", ")}`,
      ));
    }

    if (metadata.notes.length) {
      const target = section("Notes");
      for (const note of metadata.notes) target.append(element("p", "xray-note", note));
    }

    panel.hidden = false;
    constrainPanel();
    if (resizeObserver) resizeObserver.disconnect();
    resizeObserver = new ResizeObserver(scheduleOutline);
    resizeObserver.observe(node);
    scheduleOutline();
  }

  function clearSelection() {
    selected = null;
    panel.hidden = true;
    outline.classList.remove("is-visible");
    if (resizeObserver) resizeObserver.disconnect();
  }

  function synchronize() {
    root.dataset.lcmXrayActive = String(active);
    toggle.classList.toggle("is-active", active);
    toggle.textContent = active ? "X-RAY ON" : "X-RAY";
    toggle.setAttribute("aria-pressed", String(active));
    if (!active) {
      pinned = false;
      hovered = null;
      clearSelection();
    }
  }

  function setActive(value) {
    active = Boolean(value);
    synchronize();
  }

  function pointerMove(event) {
    if (!active || event.composedPath().includes(host)) return;
    const candidate = inspectableFromEvent(event);
    if (candidate === hovered) return;
    hovered = candidate;
    window.clearTimeout(hoverTimer);
    if (!candidate) {
      if (!pinned) clearSelection();
      return;
    }
    if (pinned) return;
    hoverTimer = window.setTimeout(() => render(candidate), 110);
  }

  function captureSelection(event) {
    if (!active || event.composedPath().includes(host)) return;
    const candidate = inspectableFromEvent(event);
    if (!candidate) return;
    event.preventDefault();
    event.stopImmediatePropagation();
    pinned = true;
    render(candidate);
  }

  function keyboard(event) {
    if (event.ctrlKey && event.shiftKey && event.key.toLowerCase() === "x") {
      event.preventDefault();
      setActive(!active);
      return;
    }
    if (event.key === "Escape" && active) {
      event.preventDefault();
      pinned = false;
      clearSelection();
    }
  }

  toggle.addEventListener("click", event => {
    event.stopPropagation();
    setActive(!active);
  });
  close.addEventListener("click", event => {
    event.stopPropagation();
    pinned = false;
    clearSelection();
  });
  panelHeader.addEventListener("pointerdown", event => beginPanelInteraction(event, "move"));
  panelHeader.addEventListener("pointermove", updatePanelInteraction);
  panelHeader.addEventListener("pointerup", endPanelInteraction);
  panelHeader.addEventListener("pointercancel", endPanelInteraction);
  resizeHandle.addEventListener("pointerdown", event => beginPanelInteraction(event, "resize"));
  resizeHandle.addEventListener("pointermove", updatePanelInteraction);
  resizeHandle.addEventListener("pointerup", endPanelInteraction);
  resizeHandle.addEventListener("pointercancel", endPanelInteraction);
  document.addEventListener("pointermove", pointerMove, { capture: true, signal });
  document.addEventListener("click", captureSelection, { capture: true, signal });
  document.addEventListener("keydown", keyboard, { capture: true, signal });
  document.addEventListener("scroll", scheduleOutline, { capture: true, passive: true, signal });
  window.addEventListener("resize", () => {
    constrainPanel();
    scheduleOutline();
  }, { passive: true, signal });

  const controller = {
    root,
    owns,
    accept,
    enable: () => setActive(true),
    disable: () => setActive(false),
    toggle: () => setActive(!active),
    destroy() {
      abort.abort();
      window.clearTimeout(hoverTimer);
      if (frame) cancelAnimationFrame(frame);
      if (resizeObserver) resizeObserver.disconnect();
      host.remove();
      controllers.delete(controller);
    },
  };
  controllers.add(controller);

  for (let index = pending.length - 1; index >= 0; index -= 1) {
    if (!owns(pending[index].node)) continue;
    accept(pending[index]);
    pending.splice(index, 1);
  }

  globalThis.lcmXRay = controller;
  synchronize();
  return controller;
}

globalThis.LineCableModelsComponentXRay = { installXRay };
})();
