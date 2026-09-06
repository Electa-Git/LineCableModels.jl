(() => {
'use strict';
// Multiple diagnostic hosts in one document must share preview rule identities.
if (globalThis.LineCableModelsCSSPreview) return;
// CSSOM siblings preserve the authored cascade and enclosing conditions. Source
// declarations and files are never rewritten. Only our scoped rules are removed.
const previewRules = new WeakSet();
const ruleIds = new WeakMap();
let ruleSequence = 0;
const units = new Set(['', 'px', 'rem', 'em', '%', 'vw', 'vh', 'dvh', 'ch']);
const simpleNumber = /^(-?(?:\d+\.?\d*|\.\d+))([a-z%]*)$/i;

function selectors(text) {
  const result = [];
  let start = 0, depth = 0, quote = '';
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (ch === '\\') { i++; continue; }
    if (quote) { if (ch === quote) quote = ''; continue; }
    if (ch === '"' || ch === "'") quote = ch;
    else if (ch === '(' || ch === '[') depth++;
    else if (ch === ')' || ch === ']') depth--;
    else if (ch === ',' && !depth) { result.push(text.slice(start, i).trim()); start = i + 1; }
  }
  result.push(text.slice(start).trim());
  return result;
}

function ownedSelectors(selector, scopes) {
  const classes = new Set(scopes.filter(scope => /^\.[\w-]+$/.test(scope)).map(scope => scope.slice(1)));
  return selectors(selector).filter(part => {
    // Class tokens inside attribute strings are not selector ownership.
    const structural = part.replace(/\[(?:[^\]"']|"[^"]*"|'[^']*')*\]/g, '');
    return [...structural.matchAll(/\.([\w-]+)/g)].some(match => classes.has(match[1]));
  });
}

function selectorReason(parts) {
  if (!parts.length) return 'Ownership needs an explicit class scope.';
  if (parts.some(part => /[\\&]|:visited\b|:host\b|::(?!before\b|after\b)/.test(part)))
    return 'This selector requires unsupported nesting, escaping or pseudo-element handling.';
  return '';
}

function scopedSelector(parts, id) {
  const marker = `[data-lcm-preview-scope~="${id}"]`;
  return parts.map(part => {
    const pseudo = part.match(/::(?:before|after)\s*$/)?.[0] || '';
    const base = pseudo ? part.slice(0, -pseudo.length) : part;
    return `:is(${base}):where(${marker}, ${marker} *)${pseudo}`;
  }).join(', ');
}

function declarations(style) {
  // CSSStyleDeclaration iteration expands shorthand into generated longhands.
  // Its serialized declarations retain shorthand and never include inherited CSS.
  const result = [];
  let start = 0, depth = 0, quote = '';
  const text = style.cssText + ';';
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (ch === '\\') { i++; continue; }
    if (quote) { if (ch === quote) quote = ''; continue; }
    if (ch === '"' || ch === "'") quote = ch;
    else if (ch === '(' || ch === '[') depth++;
    else if (ch === ')' || ch === ']') depth--;
    else if (ch === ';' && !depth) {
      const entry = text.slice(start, i), colon = entry.indexOf(':');
      if (colon > 0) {
        const property = entry.slice(0, colon).trim();
        result.push({property, value: style.getPropertyValue(property).trim(),
          priority: style.getPropertyPriority(property)});
      }
      start = i + 1;
    }
  }
  return result;
}

function collect(node, scopes) {
  const records = [];
  const duplicates = new Map();
  const visit = (parent, source, conditions = [], blocked = '') => {
    for (const rule of [...parent.cssRules]) {
      if (previewRules.has(rule)) continue;
      if (rule.selectorText && rule.style) {
        const parts = ownedSelectors(rule.selectorText, scopes);
        if (!parts.length) continue;
        if (!ruleIds.has(rule)) ruleIds.set(rule, ++ruleSequence);
        const signature = JSON.stringify([source.name, conditions, rule.selectorText, rule.style.cssText]);
        if (duplicates.has(signature)) {
          duplicates.get(signature).copies.push({rule, parent}); continue;
        }
        const record = {key: ruleIds.get(rule), rule, parent, source, conditions, parts, copies: [{rule, parent}],
          selector: rule.selectorText,
          reason: blocked || selectorReason(parts) || (!source.owned ? 'Unidentified stylesheet; preview is read-only.' : ''),
          declarations: declarations(rule.style)};
        duplicates.set(signature, record); records.push(record);
      } else if (rule.cssRules && !/^@(?:-\w+-)?keyframes\b/.test(rule.cssText)) {
        const condition = rule.cssText.slice(0, rule.cssText.indexOf('{')).trim();
        visit(rule, source, [...conditions, condition], blocked ||
          (/^@(media|supports|container|layer)\b/.test(condition) ? '' : 'Unsupported grouping rule; read-only.'));
      }
    }
  };
  for (const [index, sheet] of [...node.ownerDocument.styleSheets].entries()) {
    if (sheet.ownerNode?.hasAttribute('data-lcm-xray-preview')) continue;
    const tagged = sheet.ownerNode?.getAttribute('data-lcm-css-source');
    const url = sheet.href ? new URL(sheet.href, location.href) : null;
    const name = tagged || (url?.origin === location.origin ? url.pathname : `Unidentified inline stylesheet ${index + 1}`);
    try { visit(sheet, {name, owned: Boolean(tagged || url?.origin === location.origin)}); }
    catch (_) { /* Cross-origin or unavailable sheets cannot be edited. */ }
  }
  if (node.style.length) records.unshift({key: 'inline', selector: 'element.style',
    source: {name: 'Application-managed inline declarations', owned: false}, conditions: [],
    reason: 'Inline styles may be driven by application state; read-only.',
    declarations: declarations(node.style)});
  return records;
}

function validate(property, value, spec) {
  if (!spec || spec.kind === 'readonly') return 'Property is not enabled for CSS preview.';
  // CSS.supports validates syntax, not safety. This is deliberately narrower.
  if (!value || value.length > 300 || /[;{}\\]|\/\*|url\s*\(|image-set\s*\(|attr\s*\(/i.test(value))
    return 'Use a CSS value without URLs, escapes, comments or extra declarations.';
  if (!CSS.supports(property, value)) return 'The browser does not accept this CSS value.';
  const number = value.match(simpleNumber);
  if (number) {
    const numeric = Number(number[1]), unit = number[2];
    if (!Number.isFinite(numeric)) return 'Enter a finite number.';
    if (!units.has(unit) || (spec.kind === 'number' && unit) ||
        (spec.kind === 'length' && unit && !spec.units.includes(unit))) return 'Unsupported unit.';
    if (spec.minimum !== null && numeric < spec.minimum) return `Minimum is ${spec.minimum}.`;
    if (spec.maximum !== null && numeric > spec.maximum) return `Maximum is ${spec.maximum}.`;
  }
  if (spec.kind === 'choice' && !spec.choices.includes(value) && !/^(var\(--lc-[\w-]+\)|inherit|initial|unset|revert)$/.test(value))
    return 'Choose one of the declared CSS values.';
  if (spec.kind === 'color' && !CSS.supports('color', value)) return 'Choose a color or semantic color token.';
  return '';
}

class Preview {
  constructor(root, catalogue, changed = () => {}) {
    this.root = root; this.catalogue = catalogue; this.changed = changed;
    this.groups = new Map(); this.nodes = new Map(); this.enabled = true;
    this.observer = new MutationObserver(() => {
      for (const node of this.nodes.keys()) if (!node.isConnected) this.reset(node);
    });
    this.observer.observe(root, {childList: true, subtree: true});
  }
  editor(record, declaration, overrides = {}) {
    if (record.reason) return {kind: 'readonly', reason: record.reason};
    if (!declaration.value) return {kind: 'readonly', reason: 'The browser does not expose an authored value for this shorthand-derived declaration.'};
    const base = this.catalogue[declaration.property];
    if (!base) return {kind: 'readonly', reason: 'Not in the CSS preview allowlist (including application variables and transforms).'};
    return overrides[declaration.property] || base;
  }
  state(node, record, property) {
    return this.groups.get(node)?.get(record.key)?.edits.get(property);
  }
  mark(node) {
    if (!this.nodes.has(node)) {
      const id = 'xp-' + crypto.randomUUID();
      this.nodes.set(node, {id, original: node.getAttribute('data-lcm-preview-scope')});
      node.setAttribute('data-lcm-preview-scope', [node.getAttribute('data-lcm-preview-scope'), id].filter(Boolean).join(' '));
    }
    return this.nodes.get(node).id;
  }
  removeRule(group) {
    if (!group.inserted) return;
    for (const {rule, parent} of group.inserted) {
      const index = [...parent.cssRules].indexOf(rule);
      if (index >= 0) parent.deleteRule(index);
    }
    group.inserted = null;
  }
  apply(group) {
    const {record, node, edits} = group;
    const values = [...edits].filter(([, edit]) => edit.enabled);
    if (!this.enabled || !values.length) { this.removeRule(group); return; }
    const selector = scopedSelector(record.parts, this.mark(node));
    const declarations = values.map(([property, edit]) =>
      `${property}: ${edit.value}${edit.priority ? ' !important' : ''};`).join('\n');
    // Insert before removing the previous override so failure leaves it intact.
    const inserted = [];
    try {
      for (const {rule, parent} of record.copies) {
        const index = [...parent.cssRules].indexOf(rule);
        if (index < 0) throw Error('The source rule was replaced; reset and select the component again.');
        const at = parent.insertRule(`${selector} { ${declarations} }`, index + 1);
        const copy = parent.cssRules[at];
        previewRules.add(copy); inserted.push({rule: copy, parent});
      }
    } catch (error) { this.removeRule({inserted}); throw error; }
    this.removeRule(group);
    group.inserted = inserted;
  }
  set(node, record, declaration, value, spec, enabled = true) {
    if (!this.root.contains(node) && node !== this.root) return 'Component is outside this X-ray host.';
    if (record.reason) return record.reason;
    if (!declaration.value) return 'An authored default is unavailable; this declaration is read-only.';
    if (value === declaration.value) { this.resetProperty(node, record, declaration.property); return ''; }
    if (!this.groups.has(node)) this.groups.set(node, new Map());
    this.mark(node);
    const groups = this.groups.get(node);
    if (!groups.has(record.key)) groups.set(record.key, {node, record, edits: new Map(), inserted: null});
    const group = groups.get(record.key), previous = group.edits.get(declaration.property);
    const reject = error => {
      // Keep the last valid preview, but retain the rejected draft so changing
      // selection cannot hide it and every reset includes it.
      group.edits.set(declaration.property, {...(previous || {
        value: declaration.value, enabled: false, priority: declaration.priority, original: declaration.value,
      }), draft: value, draftEnabled: enabled, error});
      this.changed(); return error;
    };
    const error = validate(declaration.property, value, spec);
    if (error) return reject(error);
    group.edits.set(declaration.property, {value, enabled, priority: declaration.priority, original: declaration.value});
    try { this.apply(group); }
    catch (failure) { return reject(failure.message); }
    this.changed(); return '';
  }
  resetProperty(node, record, property) {
    const group = this.groups.get(node)?.get(record.key);
    if (!group) return;
    group.edits.delete(property); this.apply(group);
    if (!group.edits.size) this.groups.get(node).delete(record.key);
    if (!this.groups.get(node).size) this.reset(node);
    this.changed();
  }
  reset(node, {recursive = false} = {}) {
    // Snapshot targets before removing anything; notify once for a whole group.
    const targets = [...this.nodes.keys()].filter(candidate =>
      !node || candidate === node || (recursive && node.contains(candidate)));
    for (const target of targets) {
      for (const group of this.groups.get(target)?.values() || []) this.removeRule(group);
      this.groups.delete(target);
      const mark = this.nodes.get(target);
      if (mark) {
        const remaining = (target.getAttribute('data-lcm-preview-scope') || '').split(' ').filter(id => id !== mark.id).join(' ');
        if (remaining) target.setAttribute('data-lcm-preview-scope', remaining); else target.removeAttribute('data-lcm-preview-scope');
      }
      this.nodes.delete(target);
    }
    this.changed();
  }
  toggle(enabled) {
    this.enabled = enabled;
    for (const groups of this.groups.values()) for (const group of groups.values()) this.apply(group);
    this.changed();
  }
  changes(node, {recursive = false} = {}) {
    const changes = [];
    for (const [candidate, groups] of this.groups) {
      if (node && candidate !== node && !(recursive && node.contains(candidate))) continue;
      for (const group of groups.values()) {
        for (const [property, edit] of group.edits) changes.push({component: group.node.dataset.lcmComponent,
          source: group.record.source.name, selector: group.record.parts.join(', '),
          conditions: group.record.conditions, property, ...edit});
      }
    }
    return changes;
  }
  export() {
    // Invalid drafts are not CSS proposals. The UI also disables copying until
    // they are corrected or reset, so a last-valid value cannot masquerade as a draft.
    return this.changes().filter(change => !change.error).map(change => {
      const clean = text => String(text).replaceAll('*/', '* /').replaceAll('\n', ' ');
      const note = `/* ${clean(change.source)} | ${clean(change.component)} | was: ${clean(change.original)}${change.enabled ? '' : ' | preview disabled'} */`;
      let text = `${change.selector} {\n  ${change.property}: ${change.value}${change.priority ? ' !important' : ''};\n}`;
      for (const condition of [...change.conditions].reverse()) text = `${condition} {\n${text.split('\n').map(line => '  ' + line).join('\n')}\n}`;
      return note + '\n' + text;
    }).join('\n\n');
  }
  destroy() { this.observer.disconnect(); this.reset(); }
}
globalThis.LineCableModelsCSSPreview = {collect, Preview, validate, simpleNumber};
})();
