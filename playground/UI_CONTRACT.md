# UI ownership and verification

The gallery and live applications are consumers of the same components. A
successful calculation alone is not acceptance of an application UI.

Published documents have one column owner: `main#quarto-document-content` in
`theme.scss`. It owns the maximum content width, gutters and fixed-header
clearance. Titles, prose, cards, catalogues, sections and code blocks flow inside
that same column; page classes must not independently centre or constrain them.
The utility header uses those same column tokens. Layout classes describe
composition, not navigation identity. Bounded workbenches, widgets and slides
retain their own viewport contracts rather than inheriting document geometry.

| Concern | Owner | Consumer responsibility |
| --- | --- | --- |
| Palette, including Bootstrap RGB channels | `assets/brand.css` | Use semantic tokens; do not redeclare colours. |
| Theme selection, storage and system preference | `assets/theme-init.html` | Use `theme_script()` in Bonito; forward deck theme messages through its API. |
| Published document shell | `assets/theme.scss` | Supply Quarto navigation/content, independent of menu depth. |
| Published live/inactive embeds | `bonito.lua`, `assets/published-frames.css` / `.js` | The shortcode owns one reserved rectangle. Its inactive message uses the shared inset; live frames have zero outer padding. Never insert bare notices into a live canvas host. |
| Page heading, density and content area | `Toolkit.WorkspacePage`, `assets/workspace.css` | Supply title, content and optional heading slots; use `fill=true` for a bounded viewport. |
| Resizable regions | `WorkbenchUI.SplitPane`, `assets/split-pane.css` | Compose existing views; do not copy its CSS or drag handlers. |
| Dock and inspector insets | `workbench.css`, `--lc-content-inset` in `brand.css` | Dock tab contents supply content, not outer padding; all tabs share the shell inset. |
| Indeterminate activity | `assets/control-contract.css` | Status lines use `lc-activity-status` and explicit `data-busy="true"`; viewport and toolbar indicators share its glyph, theme tokens and reduced-motion policy. Clear busy when evidence becomes stale or work ends. |
| Semantic status | `Toolkit.StatusIndicator`, `control-contract.css`, `brand.css` | Use a visible label plus bold text and `data-tone="neutral\|info\|success\|warning\|danger"`. Success uses the green `--lc-success` token, not an application accent. Unknown/stale is not online. Colours alone never encode the state. |
| Action feedback | `Toolkit.ActionButton`, `control-contract.css` | The action owner supplies busy and disabled state. Browser controls use the same `.lc-button[data-busy]` and `aria-busy` contract with an explicit busy label. Do not invent progress percentages or automatically replay uncertain actions. |
| Fields, buttons, tables, disclosures, frames | `Toolkit`, `forms.css`, `data-views.css` | Reuse the components. Browser-only runtime controls use the same structural classes and styles. |
| Scientific drawings | `ScientificViews`, `scientific-views.css` | Own scientific geometry/series only, not field/button/page styling. |

`ViewportFrame(...; sizing=:content)` fits ordinary forms; the default
`:viewport` reserves canvas height. `sizing=:fill` fits a drawing to its host's
available height; widening a split must not make the drawing taller. An optional
`max_height="32rem"` caps the whole frame (header and footer included). Canvas
frames retain their 12rem usable minimum, so smaller hosts scroll rather than
clip content. For example:

```julia
WorkspacePage("Construction",
    SplitPane(ViewportFrame("Cross-section", drawing; sizing=:fill),
        ViewportFrame("Inputs", fields; sizing=:content); scroll=:parent);
    fill=true)
```

Use `SplitPane(...; scroll=:parent)` for canvas/form compositions: the workbench
view (or standalone document) owns overflow at its right edge, including long
forms, expanded disclosures and the narrow stacked layout. Its wrappers must
not introduce inner pane scrollbars or hide overflow. The default `scroll=:panes`
remains available for deliberately independent editors/logs. The diagnostics
dock is a separate region and keeps its own scrolling. None of these policies
depends on colour theme. Component X-ray metadata belongs to the
component defining each style, rather than duplicating descendants' metadata.
Reactive Bonito `data-*` and `aria-*` state must synchronize HTML attributes
explicitly (`onjs` / `setAttribute`); a DOM-property update alone does not update
CSS selectors or accessibility state. The real-host feedback test covers both
entering and leaving the busy state.

## Document and return-navigation ownership

Published pages allocate the footer's natural height through the shared flex
shell, never a fixed subtraction from the screen height. Components must not
append trailing margins to the document's own bottom inset. A fitting document
has no scroll range; longer content and expanded navigation remain scrollable.

Workbench return navigation uses `Toolkit.NavigationButton`, with button styles
owned by `forms.css`. The sidebar only controls placement and label visibility
in an icon rail; the accessible name and icon remain available. For example:

```julia
CableStudy.app(client; return_button=NavigationButton("Dashboard";
    href="/workbenches/", icon=WorkbenchUI.icon(:workbench)))
```

The default is Home at `/`. Root-relative destinations keep navigation in the
owned site; native links work without a Julia callback or worker connection.

## Runtime status and diagnostic ownership

Terminal loading describes connection/start/restart/stop, not individual input,
output or keepalive requests. The terminal badge reserves its activity glyph's
space and only changes input/cursor options when their values change. Ordinary
typing must not resize the terminal or shift its controls; delayed-acknowledgement
browser tests sample geometry throughout typing in both themes.

Broker connection, worker health, application-run lifetime, assignment, and
executor preparation are independent facts. The shared browser client reads the
owned run as well as inventory: an online worker cannot make a stopped/failed
run accept assignments. Unknown run evidence disables admission. An ended run
offers its existing status page, where a new run can be explicitly started; no
widget silently replaces a run or bypasses the server's authority checks.

All runtime controls in one document and run share one client and a bounded
256-entry action history. Explicit refresh records start and completion/failure;
assignment, preparation and terminal controls record action/phase outcomes.
Background polls record connection/run transitions, not one log per poll.
Entries contain bounded status labels and request identities, never scientific
inputs, credentials, terminal keystrokes or output. This page-local history is
not durable and is not a replacement for operator logs.

The separate **Control events · server** disclosure retains the coordinator's
bounded event stream and labels gaps honestly. In an owned run it displays that
run's events plus shared worker/connection events; inventory-only controls show
all events authorized for the account. Worker inventory does not claim
per-worker preparation: readiness belongs to a particular assigned executor.
CableStudy owns one `WorkerDiagnostics` in its dock. Its `StudyRuntime` content
uses `diagnostics=false`; standalone/presentation consumers retain their own
diagnostic disclosure. New applications must choose one owner per diagnostic
scope rather than mounting the same inventory/log twice.

## Acceptance checks

Run source and rendering tests with `julia --project=playground playground/test/runtests.jl`.
The dedicated app comparison uses real registered UI drivers with a temporary
database, no broker and no scientific worker allocations:

```bash
LCM_RUNTIME_BROWSER_SUITE=parity bash playground/runtime/test/run-bonito.sh
LCM_RUNTIME_BROWSER_SUITE=scientific bash playground/runtime/test/run-bonito.sh
```

The parity test compares computed typography, padding, borders and colours
between the gallery, reference workbench, CableStudy and scientific frames. It
also checks resizing, narrow hosts and preservation of mounted plots. Keep
functional input, failure/recovery, presentation and X-ray tests alongside it.
The construction regression samples every drag step, frame-height caps, dock
expansion, short windows and narrow stacking in both themes. It rejects inner
pane scrolling, clipped or distorted drawings, and unreachable form content.
It also checks dock insets for every reference tab and the actual CableStudy
diagnostics, plus runtime startup/stop/recovery/deadline states using intercepted
HTTP evidence (no worker allocation). Startup elapsed time explicitly measures
time on the page, not server progress. Only reported worker progress is shown as
a percentage. Status updates remain accessible without continuous timer announcements.
The published-shell browser test compares actual left/right content edges and
header clearance on every navigation route in both themes at eight viewport
sizes (390–2560px, including short projection screens). Review rendered
screenshots as well: token equality and a lack of overflow alone do not prove
that a composition is visually correct.
Header children must fit inside its border, not paint over it. The browser alone
allocates the document scrollbar: do not reserve a second/empty root gutter.
Tests enforce these painted boundaries and the actual text insets of every
inactive embed. Run the shell check against a runtime gateway as well as the
standalone publisher: one requires an explicit toolkit run, the other serves
live widget routes directly. Both modes are part of acceptance.
The repeated document-geometry matrix uses inert iframe bodies in standalone
mode; it does not claim to exercise Julia widget internals. Actual widget
rendering is checked separately by `visual_contract_browser.mjs` and the real
deck/workbench/gallery launch tests in `runtime/test/catalogue_browser.mjs`.
Capability-check failures retain that same inactive frame instead of attempting
a different route namespace. Live frames establish their visible rectangle
before assigning `src`, preserving native lazy loading for the gallery.
The widget browser test checks all 22 registered specimens at 390, 768, 1024
and 1600px in both themes. Set `LCM_VISUAL_ARTIFACTS` to an output directory to
retain desktop and narrow screenshots from that test. Stateful hover/selection,
Source dialogs and the sidebar's owned canvas are checked explicitly; a theme
switch must also update vendor-generated publisher surfaces.

## Publishing updates

`lcm playground build` and `lcm presentation build` commit only successful
renders. Both publishers retain bounded byte snapshots and adopt the next
committed revision without restarting UI hosts. Failed renders retain the last
valid publication; the previous revision's non-HTML assets remain available for
in-flight requests. Run builds serially, and use these commands rather than
editing `_site` or calling Quarto directly against a running publication.

Publication limits are 64 MiB per file and 128 MiB per revision, with at most the
current and previous revisions held. Julia application-code or embedded-style
changes still require a new application run; publishing does not mutate an
existing user's session.
