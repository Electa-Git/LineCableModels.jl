# Owned-component X-ray and CSS preview

Start with `lcm playground start --xray`, open a workbench, and enable X-ray
with its top-bar button or `Ctrl+Shift+X`. Diagnostics are opt-in; no metadata
is emitted without host permission.

## Use it

1. **Hover** outlines an inspectable component. It does not populate or change
   the diagnostic window.
2. **Click** selects that instance and loads its metadata and owned CSS.
   Moving the mouse elsewhere leaves selection and drafts alone.
3. CSS numbers have spinners and unit selectors; keywords use dropdowns;
   colors offer semantic brand tokens. Expression mode keeps `var(...)`,
   `calc(...)`, percentages and compound values explicit instead of replacing
   them with computed pixels. Checkboxes enable/disable individual overrides;
   CSS keyword states use dropdowns, not invented Boolean properties.
4. The mode button shows **Pick components** or **Interact with application**.
   Click it to switch. Pick intercepts component clicks so inspection cannot
   accidentally invoke a command; Interact operates the application without
   changing the selected component. Code metadata and bindings remain read-only.
5. **Preview** compares originals with drafts. Edited properties have an accent
   border and an explicit **Original → Override** comparison. Disabled overrides
   remain highlighted; invalid drafts show the rejected value and last valid
   override, and remain visible when you return to that component.
6. **Reset component** restores only the selected instance's own CSS overrides
   and drafts. Check **Apply to children** to include its entire descendant tree,
   including hidden and deeply nested components. Breadcrumbs select a group or
   the whole Workbench; the status line shows selected, child, and total counts.
7. **Reset all** always clears every X-ray override and draft in this host,
   regardless of selection, **Apply to children**, or whether preview is paused.
   Per-property **Reset** restores that declaration alone. Reset restores the
   authored stylesheet, not browser defaults, and preserves the selected theme.
8. **Copy changes** exports source
   file, original selector, enclosing conditions, old value and proposed CSS.
   Review changes in the owning stylesheet; exporting never saves source files.
   Correct or reset invalid drafts before copying; rejected text is never
   exported as executable CSS. Editing again clears an outdated export preview.

The window remains movable and resizable. Breadcrumbs select ancestors.
Closing it, Escape, disabling X-ray, removing the component, navigating away,
or reloading discards temporary previews. No preview is persisted to storage
or sent to Julia, NATS, or a worker. Normal application actions in Interact mode
retain their existing callbacks and may perform their normal work.
X-ray resets only its own CSS previews and drafts. It does not reset entered
application data, native splitter/dock positions, or changes made outside X-ray.

## Reuse and specialize

The shared implementation lives in `src/diagnostics/ComponentXRay.jl`,
`CssEditors.jl`, `css_preview.js`, and `component_xray.js`. Gallery and workbench
hosts use that same engine; there is no per-page editor or parallel theme.
Controls consume `assets/brand.css` and `assets/control-contract.css` even
inside the inspector's isolated shadow root.

Keep the existing `inspection(::MyComponent)` registration with explicit
simple-class `css_scopes`. Most components need no new Julia code: the central
property catalogue selects an editor. Optional specialization uses dispatch:

```julia
import LineCableModelsPlayground.ComponentXRay: CssEditor, css_editors

css_editors(::MyComponent) = Dict(
    "gap" => CssEditor(:length; units=["px", "rem"], minimum=0, maximum=48, step=1),
    "display" => CssEditor(:choice; choices=["grid", "flex"]),
    "width" => CssEditor(:readonly),
)
```

`ComponentInspection` consumes these hints automatically. Its `css_overrides`
keyword can explicitly supply hints when needed. Hints describe editors, not
current CSS values; the owning stylesheet remains the source of declarations.
Bounds constrain simple numeric values, not the result of arbitrary expressions.
Component overrides cannot enable properties outside the shared allowlist.
`XRayPolicy(; permitted=true, enabled=true, css_preview=false)` retains the
read-only inspector for hosts that should not expose CSS preview.

Embedded first-party styles carry `data-lcm-css-source="path/to/owner.css"`.
Same-origin linked styles use their URL. Unidentified inline stylesheet blocks,
cross-origin sheets, and application-managed `element.style` are not editable.
CSSOM serialization may normalize author syntax; the inspector shows declared
rules and shorthand, not a computed/inherited style dump or claimed source lines.
If the browser exposes blank values for shorthand-derived declarations, the
inspector reports them as unavailable rather than offering invented defaults.

## Cascade and safety limits

Each draft inserts an instance-scoped CSSOM sibling immediately after its
original rule, within the same grouping. Original declarations are untouched.
The scope adds zero specificity; selector state, source order, `!important`,
media/support/container conditions and layers remain in force. Repeated copies
of a portable component's stylesheet are handled together. Consequently a
declaration can still be inactive or overridden by a later, stronger rule.
The window shows selectors and conditions; it does not force hover/focus states.

Reset removes only rules and scope tokens owned by that preview. Previewing a
parent's descendant selector affects matching descendants of that instance,
not other instances elsewhere. Colors change a declaration's token reference;
the editor never redefines global theme variables.
Component-only reset preserves descendant-owned drafts; removing a parent's
declaration can still change the appearance inherited by its descendants.

V1 deliberately leaves application variables, transforms, unsupported
properties, ambiguous/nested selectors, and unsupported grouping rules read-only.
This protects application-controlled geometry and the presentation contract.
Values must pass both the property allowlist and browser CSS validation; URLs,
escapes, comments and injected extra declarations are rejected. This is a
developer diagnostic tool, not an authorization boundary or an arbitrary CSS
editor for untrusted users. No source write-back, bulk/global edits, animation
editor, forced pseudo-states, or undo history is implemented in this version.

## Regression gate

```sh
julia --startup-file=no --project=. test/runtests.jl
bash test/integration/run-xray.sh
```

The browser gate uses an isolated profile, local Julia mock fixture, and the
real template workbench. It checks click-only selection, instance isolation,
repeated stylesheets, conditional rules, theme changes, read-only metadata,
live bindings without input replacement, validation, export, reset, and teardown.
A nested-tree fixture checks component-only, recursive, hidden-descendant, and
global reset isolation plus changed/disabled/invalid comparisons in both themes.
