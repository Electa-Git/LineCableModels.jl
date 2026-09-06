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
  const previewBar = element('div', 'xray-preview-bar');
  const interaction = element('button', 'xray-action', 'Pick components');
  interaction.type = 'button';
  interaction.title = 'Pick: clicks inspect. Interact: clicks operate the application.';
  interaction.setAttribute('aria-pressed', 'false');
  const previewLabel = element('label', 'xray-preview-switch');
  const previewToggle = element('input');
  previewToggle.type = 'checkbox'; previewToggle.checked = true;
  previewToggle.setAttribute('aria-label', 'Enable temporary CSS preview');
  previewLabel.append(previewToggle, 'Preview');
  const childrenLabel = element('label', 'xray-preview-switch');
  const applyChildren = element('input'); applyChildren.type = 'checkbox';
  applyChildren.setAttribute('aria-label', 'Apply to children');
  childrenLabel.title = 'Reset component also clears previews and drafts in all descendant components.';
  childrenLabel.append(applyChildren, 'Apply to children');
  const resetComponent = element('button', 'xray-action', 'Reset component');
  const resetAll = element('button', 'xray-action', 'Reset all');
  resetAll.title = 'Restore authored CSS everywhere in this X-ray host, including all children. Application data is unchanged.';
  const copyChanges = element('button', 'xray-action', 'Copy changes');
  for (const button of [resetComponent, resetAll, copyChanges]) button.type = 'button';
  const previewStatus = element('output', 'xray-preview-status', 'Temporary CSS only · no source changes');
  previewStatus.setAttribute('aria-live', 'polite');
  previewBar.append(interaction);
  if (options.css_preview !== false) previewBar.append(previewLabel, childrenLabel, resetComponent, resetAll, copyChanges, previewStatus);
  const resizeHandle = element("div", "xray-resize-handle");
  resizeHandle.title = "Resize diagnostics window";
  resizeHandle.setAttribute("aria-hidden", "true");
  panel.append(panelHeader, breadcrumbs, previewBar, body, resizeHandle);
  shadow.append(outline, toggle, panel);

  const registry = new WeakMap();
  let active = Boolean(options.enabled);
  let hovered = null;
  let selected = null;
  let interacting = false;
  const bindingCells = new Map();
  const cssFieldRefresh = new Set();
  let sizeCell = null;
  let frame = 0;
  let resizeObserver = null;
  let layoutObserver = null;
  let panelPositioned = false;
  let panelInteraction = null;
  const abort = new AbortController();
  const signal = abort.signal;
  const preview = new globalThis.LineCableModelsCSSPreview.Preview(root, options.css_editors || {}, () => {
    if (selected && !selected.isConnected) {
      selected = null; panel.hidden = true; bindingCells.clear(); cssFieldRefresh.clear();
      if (resizeObserver) resizeObserver.disconnect();
    }
    updatePreviewStatus();
    body.querySelector('.xray-export')?.remove();
    for (const refresh of cssFieldRefresh) refresh();
    scheduleOutline();
  });
  resetComponent.disabled = resetAll.disabled = copyChanges.disabled = true;

  function updatePreviewStatus() {
    const changes = preview.changes();
    const invalid = changes.filter(change => change.error).length;
    const scoped = selected ? preview.changes(selected, {recursive: applyChildren.checked}).length : 0;
    const descendants = selected ? preview.changes(selected, {recursive: true}).length - preview.changes(selected).length : 0;
    previewStatus.textContent = `Temporary CSS · ${scoped} in selection${applyChildren.checked ? ' + children' : ''}` +
      (!applyChildren.checked && descendants ? ` · ${descendants} in children` : '') +
      ` · ${changes.length} total${invalid ? ` · ${invalid} invalid` : ''} · no source writes`;
    resetComponent.title = applyChildren.checked
      ? 'Restore authored CSS for this component and all descendants, including hidden components.'
      : 'Restore only this component’s own previews and drafts; descendants are unchanged.';
    resetComponent.disabled = scoped === 0;
    resetAll.disabled = changes.length === 0;
    copyChanges.disabled = changes.length === 0 || invalid > 0;
    copyChanges.title = invalid ? 'Correct or reset invalid drafts before copying changes.' : 'Copy proposed CSS for manual source review.';
  }

  const queryFlag = new URLSearchParams(window.location.search).get("xray");
  if (queryFlag !== null) active = isTruthyFlag(queryFlag);

  function positionToggle() {
    const inspector = root.querySelector(".lc-wb-inspector");
    const inspectorBounds = inspector?.getBoundingClientRect();
    const inspectorVisible = inspectorBounds && inspectorBounds.width > 1;
    const right = inspectorVisible
      ? Math.max(8, window.innerWidth - inspectorBounds.left + 6)
      : 8;
    host.style.setProperty("--xray-toggle-right", `${right}px`);
  }
  positionToggle();
  layoutObserver = new ResizeObserver(positionToggle);
  layoutObserver.observe(root);

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
      if (entry.node === selected && bindingCells.has(entry.name))
        bindingCells.get(entry.name).textContent = binding.value;
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
    const target = hovered?.isConnected ? hovered : selected;
    if (!active || !target || !target.isConnected) {
      outline.classList.remove("is-visible");
      return;
    }
    const bounds = target.getBoundingClientRect();
    outline.style.left = `${bounds.left}px`;
    outline.style.top = `${bounds.top}px`;
    outline.style.width = `${bounds.width}px`;
    outline.style.height = `${bounds.height}px`;
    outline.classList.add("is-visible");
    const metadata = registry.get(target);
    badge.textContent = `${metadata.name}  ${Math.round(bounds.width)} × ${Math.round(bounds.height)}`;
    badge.classList.toggle("is-below", bounds.top < 30);
    if (sizeCell && selected?.isConnected) {
      const size = selected.getBoundingClientRect();
      sizeCell.textContent = `${Math.round(size.width)} × ${Math.round(size.height)} px`;
    }
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
        render(candidate);
      });
      breadcrumbs.append(button);
    });
  }

  function cssField(node, metadata, rule, declaration) {
    const {property, value, priority} = declaration;
    const spec = options.css_preview === false ? {kind: 'readonly', reason: 'CSS preview is disabled by host policy.'} :
      preview.editor(rule, declaration, metadata.css_editors || {});
    const row = element('div', 'xray-css-field');
    row.dataset.cssProperty = property; row.dataset.cssRule = rule.key;
    const authored = element('code', 'xray-authored', `${property}: ${value}${priority ? ' !important' : ''}`);
    row.append(authored);
    if (spec.kind === 'readonly') {
      row.classList.add('is-readonly');
      row.append(element('small', 'xray-readonly', spec.reason || 'Read-only by component policy.'));
      return row;
    }
    const controls = element('div', 'xray-css-controls');
    const enabled = element('input'); enabled.type = 'checkbox';
    enabled.setAttribute('aria-label', `Enable ${property} override`);
    const editor = element('div', 'xray-value-editor');
    const reset = element('button', 'xray-action', 'Reset'); reset.type = 'button';
    reset.setAttribute('aria-label', `Reset ${property}`);
    const feedback = element('small', 'xray-field-feedback');
    feedback.setAttribute('aria-live', 'polite');
    const difference = element('div', 'xray-css-diff');
    const originalValue = element('code', 'xray-original-value', value);
    const proposedValue = element('code', 'xray-proposed-value');
    const proposalLabel = element('span');
    difference.append(element('span', '', 'Original'), originalValue, element('span', '', '→'), proposalLabel, proposedValue);
    const initial = preview.state(node, rule, property);
    let draft = initial?.draft ?? initial?.value ?? value;
    enabled.checked = initial?.draftEnabled ?? initial?.enabled ?? true;
    const refresh = () => {
      const edit = preview.state(node, rule, property);
      row.classList.toggle('is-invalid', Boolean(edit?.error));
      row.classList.toggle('is-changed', Boolean(edit));
      row.classList.toggle('is-preview-disabled', Boolean(edit && (!edit.enabled || !preview.enabled)));
      difference.hidden = !edit;
      proposalLabel.textContent = edit?.error ? 'Draft (invalid)' : !edit?.enabled ? 'Override (disabled)' : 'Override';
      proposedValue.textContent = edit?.draft ?? edit?.value ?? value;
      if (proposedValue.textContent === '') proposedValue.textContent = '(empty)';
      feedback.textContent = edit?.error
        ? `${edit.error} Not applied; ${edit.enabled && preview.enabled ? `last valid override: ${edit.value}` : 'authored CSS retained'}.`
        : edit ? !edit.enabled ? 'Override disabled · authored CSS retained'
          : !preview.enabled ? 'Preview paused · authored CSS retained'
            : 'Instance override · original conditions retained' : '';
      reset.disabled = !edit;
      for (const input of editor.querySelectorAll('input, select')) input.setAttribute('aria-invalid', String(Boolean(edit?.error)));
    };
    const apply = () => { preview.set(node, rule, declaration, draft.trim(), spec, enabled.checked); refresh(); };
    function mountEditor(mode) {
      editor.replaceChildren();
      const number = draft.match(globalThis.LineCableModelsCSSPreview.simpleNumber);
      if (spec.kind === 'length' || spec.kind === 'number') {
        const numeric = mode ? mode === 'number' : Boolean(number);
        const modeSelect = element('select');
        modeSelect.className = 'lc-control-select';
        modeSelect.setAttribute('aria-label', `${property} value mode`);
        modeSelect.add(new Option('Number', 'number')); modeSelect.add(new Option('Expression', 'expression'));
        modeSelect.value = numeric ? 'number' : 'expression';
        modeSelect.addEventListener('change', () => {
          if (modeSelect.value === 'number' && !number) draft = spec.kind === 'length' ? '0px' : '0';
          mountEditor(modeSelect.value);
          apply();
        });
        editor.append(modeSelect);
        if (numeric) {
          const input = element('input'); input.type = 'number'; input.step = spec.step;
          if (spec.minimum !== null) input.min = spec.minimum;
          if (spec.maximum !== null) input.max = spec.maximum;
          input.value = number?.[1] ?? '0';
          input.setAttribute('aria-label', `${property} preview`);
          const unit = element('select'); unit.setAttribute('aria-label', `${property} unit`);
          unit.className = 'lc-control-select';
          for (const suffix of new Set(['', ...(spec.units || []), number?.[2] || ''])) unit.add(new Option(suffix || 'unitless', suffix));
          unit.value = number?.[2] || '';
          const update = () => { draft = input.value ? input.value + unit.value : ''; apply(); };
          input.addEventListener('input', update); unit.addEventListener('change', update);
          editor.append(input);
          if (spec.kind === 'length') editor.append(unit);
          return;
        }
      }
      if (spec.kind === 'choice' || spec.kind === 'color') {
        const select = element('select'); select.setAttribute('aria-label', `${property} preview`);
        select.className = 'lc-control-select';
        const choices = spec.kind === 'choice' ? spec.choices : ['transparent', 'currentColor'];
        if (spec.kind === 'color') {
          const style = getComputedStyle(node);
          for (const key of [...style].filter(name => name.startsWith('--lc-')).sort()) {
            if (CSS.supports('color', style.getPropertyValue(key).trim())) choices.push(`var(${key})`);
          }
        }
        for (const choice of new Set([value, draft, ...choices])) select.add(new Option(choice, choice));
        select.value = draft;
        select.addEventListener('change', () => { draft = select.value; apply(); });
        editor.append(select);
        if (spec.kind === 'color') {
          const swatch = element('span', 'xray-swatch'); swatch.style.backgroundColor = draft;
          select.addEventListener('change', () => { swatch.style.backgroundColor = draft; });
          editor.append(swatch);
        }
        return;
      }
      const input = element('input'); input.type = 'text'; input.value = draft;
      input.spellcheck = false; input.setAttribute('aria-label', `${property} preview`);
      input.addEventListener('input', () => { draft = input.value; apply(); });
      editor.append(input);
    }
    mountEditor();
    enabled.addEventListener('change', apply);
    reset.addEventListener('click', () => {
      preview.resetProperty(node, rule, property); draft = value; enabled.checked = true;
      mountEditor(); refresh();
    });
    cssFieldRefresh.add(refresh); refresh();
    controls.append(enabled, editor, reset); row.append(controls, difference, feedback);
    return row;
  }

  function render(node) {
    const metadata = registry.get(node);
    if (!metadata) return;
    selected = node;
    title.textContent = metadata.name;
    renderBreadcrumbs(node);
    body.replaceChildren();
    bindingCells.clear();
    cssFieldRefresh.clear();

    const bounds = node.getBoundingClientRect();
    appendRows(body, [
      ["Julia type", metadata.julia_type],
      ["Rendered size", `${Math.round(bounds.width)} × ${Math.round(bounds.height)} px`],
      ["Source", `${metadata.source.file}:${metadata.source.line}`],
      ["Module", metadata.source.module],
    ], 2);
    sizeCell = body.querySelector('.xray-grid').children[3];

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
      metadata.bindings.forEach((binding, index) => bindingCells.set(binding.name,
        target.querySelector('.xray-grid').children[index * 4 + 1]));
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

    const css = globalThis.LineCableModelsCSSPreview.collect(node, metadata.css_scopes);
    if (css.length) {
      const target = section("Component-owned CSS");
      for (const rule of css) {
        const block = element("div", "xray-css-rule");
        block.append(element("code", "xray-selector", rule.selector));
        block.append(element('small', 'xray-css-source', rule.source.name));
        if (rule.conditions.length) block.append(element('small', 'xray-css-condition', rule.conditions.join(' → ')));
        const unresolved = rule.declarations.filter(declaration => !declaration.value);
        if (unresolved.length) block.append(element('small', 'xray-readonly xray-css-unresolved',
          `${unresolved.length} shorthand-derived values are not exposed by the browser; left read-only. See the owning stylesheet.`));
        for (const declaration of rule.declarations.filter(declaration => declaration.value)) {
          block.append(cssField(node, metadata, rule, declaration));
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
    updatePreviewStatus();
    scheduleOutline();
  }

  function clearSelection() {
    preview.reset();
    selected = null;
    cssFieldRefresh.clear();
    interacting = false;
    interaction.textContent = 'Pick components';
    interaction.setAttribute('aria-pressed', 'false');
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
      hovered = null;
      clearSelection();
    }
  }

  function setActive(value) {
    active = Boolean(value);
    synchronize();
  }

  function pointerMove(event) {
    if (!active) return;
    const candidate = event.composedPath().includes(host) ? null : inspectableFromEvent(event);
    if (candidate === hovered) return;
    hovered = candidate;
    scheduleOutline();
  }

  function captureSelection(event) {
    if (!active || interacting || event.composedPath().includes(host)) return;
    const candidate = inspectableFromEvent(event);
    if (!candidate) return;
    event.preventDefault();
    event.stopImmediatePropagation();
    if (candidate !== selected || panel.hidden) render(candidate);
  }

  function keyboard(event) {
    if (event.ctrlKey && event.shiftKey && event.key.toLowerCase() === "x") {
      event.preventDefault();
      setActive(!active);
      return;
    }
    if (event.key === "Escape" && active) {
      event.preventDefault();
      clearSelection();
    }
  }

  toggle.addEventListener("click", event => {
    event.stopPropagation();
    setActive(!active);
  });
  close.addEventListener("click", event => {
    event.stopPropagation();
    clearSelection();
  });
  interaction.addEventListener('click', () => {
    interacting = !interacting;
    interaction.textContent = interacting ? 'Interact with application' : 'Pick components';
    interaction.setAttribute('aria-pressed', String(interacting));
  });
  previewToggle.addEventListener('change', () => preview.toggle(previewToggle.checked));
  applyChildren.addEventListener('change', updatePreviewStatus);
  resetComponent.addEventListener('click', () => {
    if (!selected) return;
    preview.reset(selected, {recursive: applyChildren.checked}); render(selected);
  });
  resetAll.addEventListener('click', () => { preview.reset(); if (selected) render(selected); });
  copyChanges.addEventListener('click', async () => {
    const text = preview.export();
    let output = body.querySelector('.xray-export');
    if (!output) {
      output = element('textarea', 'xray-export'); output.readOnly = true;
      output.setAttribute('aria-label', 'CSS changes for review (read-only)'); body.prepend(output);
    }
    output.value = text;
    try { await navigator.clipboard.writeText(text); previewStatus.textContent = 'Changes copied · review in the owning stylesheet'; }
    catch (_) { output.focus(); output.select(); previewStatus.textContent = 'Copy the selected changes · clipboard permission is unavailable'; }
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
  document.addEventListener('pointerdown', event => {
    if (active && !interacting && !event.composedPath().includes(host) && inspectableFromEvent(event)) {
      event.preventDefault(); event.stopImmediatePropagation();
    }
  }, {capture: true, signal});
  document.addEventListener("keydown", keyboard, { capture: true, signal });
  document.addEventListener("scroll", scheduleOutline, { capture: true, passive: true, signal });
  window.addEventListener("resize", () => {
    positionToggle();
    constrainPanel();
    scheduleOutline();
  }, { passive: true, signal });

  const controller = {
    root,
    owns,
    accept,
    preview,
    enable: () => setActive(true),
    disable: () => setActive(false),
    toggle: () => setActive(!active),
    destroy() {
      abort.abort();
      preview.destroy();
      if (frame) cancelAnimationFrame(frame);
      if (resizeObserver) resizeObserver.disconnect();
      if (layoutObserver) layoutObserver.disconnect();
      host.remove();
      controllers.delete(controller);
    },
  };
  controllers.add(controller);
  const detachObserver = new MutationObserver(() => {
    if (!root.isConnected) { detachObserver.disconnect(); controller.destroy(); }
    else if (selected && !root.contains(selected)) {
      selected = null; panel.hidden = true; bindingCells.clear(); cssFieldRefresh.clear();
      if (resizeObserver) resizeObserver.disconnect();
      scheduleOutline();
    }
  });
  detachObserver.observe(document.documentElement, {childList: true, subtree: true});
  signal.addEventListener('abort', () => detachObserver.disconnect(), {once: true});

  for (let index = pending.length - 1; index >= 0; index -= 1) {
    if (!owns(pending[index].node)) continue;
    accept(pending[index]);
    pending.splice(index, 1);
  }

  globalThis.lcmXRay = controller;
  window.addEventListener('pagehide', () => controller.destroy(), {once: true, signal});
  synchronize();
  return controller;
}

globalThis.LineCableModelsComponentXRay = { installXRay };
})();
