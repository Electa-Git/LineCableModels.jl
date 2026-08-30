# Bonito showcase continuation handoff

This report records the state and design of the Bonito showcase on
`feat/codespaces-showcase` as of 2026-08-30. It is intended to let a new working
thread continue without reconstructing the decisions that led to the current
application.

## Working branch and repository state

- All showcase work belongs on `feat/codespaces-showcase`.
- The repository-root `AGENTS.md` requires an agent to verify and, when safe,
  switch to that branch before changing files or Git state.
- The working tree also contains a large, older set of user-owned changes under
  `src/engine`, `src/earthprops`, `src/transforms`, and the main documentation.
  Those changes existed before the showcase branch was selected and are not
  part of the Bonito work.
- The untracked `docs/src/transmission-line-parameters.md` and
  `test/quality/formulation_catalogue.jl` are also user-owned engine work and
  must remain outside showcase commits.
- The showcase-specific files are `.devcontainer/`, `showcase/bonito/`,
  `showcase/legacy-pluto-showcase-intent.md`, the root branch guard in
  `AGENTS.md`, and the `launch.sh` exception in `.gitignore`.

The branch already contained these showcase commits before the deck migration:

- `a7d030b7 feat(showcase): add Bonito page authoring boundary` introduced the
  Documenter-style Bonito application, its authored pages, the local
  PowerImpedance diagram extraction, the live Julia worker, and focused tests.
- `84040f8a docs(showcase): add legacy Pluto reconstruction brief` captured the
  content and intent of the discarded historical Pluto notebooks without
  reviving their implementation.
- `8e0c08cd build(codespaces): add Julia showcase container` added the Julia
  1.12 Codespaces environment and port forwarding.
- `a5d5fb1b style(showcase): remove trailing blank line` corrected the legacy
  brief's final formatting.

The commit containing this report completes the migration from one Julia file
per rendered page to one Julia file per stateful deck.

## Product direction that is now locked in

The target is a live technical manual that looks and feels like the project's
Documenter site. It is not intended to be a conventional slide deck and it is
not a live Documenter build.

Reveal.js was removed. Its viewport transforms, overlays, navigation ownership,
and resize behavior conflicted with Bonito/WGLMakie and the Documenter-like
shell. The current application deliberately uses ordinary browser layout,
native links, hash navigation within a deck, browser fullscreen, browser print,
and a small amount of shell-owned JavaScript delivered through Bonito.

The visual appearance in `assets/theme.css` is an accepted baseline and should
not be casually restyled. In particular:

- Keep the dark Documenter-inspired sidebar, top bar, typography, spacing, and
  light Makie panels.
- Keep all application styling scoped under the showcase shell except for the
  necessary viewport reset.
- Do not add selectors that reach into Bonito or WGLMakie implementation
  internals.
- Do not restore Reveal.js or another page-owning JavaScript presentation
  framework.

The cable cross-section is intentionally a showcase-local toy renderer. It
uses native Makie primitives and copies the current project colors, but it does
not call, reproduce, or extend `LineCableModels.PlotBuilder`. Cairo parity is
not a product requirement; the Cairo export is only a smoke path.

## Runtime architecture

`showcase/bonito/serve.jl` is the server entry point. It loads
LineCableModels before WGLMakie, activates WGLMakie, mounts the routes, starts
the asynchronous application-case preparation, and shuts down cleanly on an
interrupt. With `--open`, it uses `xdg-open` to launch the default Linux
browser.

`showcase/bonito/app.jl` owns:

- discovery and validation of deck descriptors;
- route generation and the Documenter-like application shell;
- sidebar grouping and search;
- previous/next navigation across every enabled deck page;
- hash navigation among pages in the same deck;
- fullscreen, print, keyboard, and touch controls;
- browser-tab persistence of presentation mode;
- viewport reconciliation after monitor, resolution, maximize, and visual
  viewport changes;
- the connection indicator and numerical-preparation status;
- asynchronous startup preparation for resource-heavy decks.

`showcase/bonito/PageAuthoring.jl` is the authoring boundary. Deck files should
normally use its `slide`, `article`, `prose`, `webgrid`, `webpart`, `control`,
`value_list`, `diagnostic`, `status_line`, `color_key`, `local_image`,
`deck_page`, and `deck_descriptor` functions instead of manually constructing
the shell DOM.

`showcase/bonito/assets/theme.css` owns the complete visual system and responsive
behavior. Font size uses bounded viewport scaling, while browser zoom remains
the deliberate projector-room override. WGLMakie uses its supported
`resize_to = :parent` configuration; the shell never patches canvas internals.

`showcase/bonito/export.jl` activates CairoMakie and saves the default cable toy
figure to `showcase/bonito/build/cable-design.svg`.

## Deck and page model

Every `showcase/bonito/pages/*.jl` file is one self-contained Julia module and
one discoverable deck. The file's final expression returns a `DECK` descriptor.
A deck has:

- a globally unique lowercase-hyphenated `id`;
- a sidebar `group` shared by any number of decks;
- a display `title` and numeric `order`;
- a deck-level `render` flag;
- an optional `setup(session)` function; and
- one or more `deck_page` entries.

`setup(session)` runs once when that deck route creates a browser session. Its
return value is given to every page builder in the file, so widgets,
observables, figures, callbacks, caches, and numerical models can be shared
across all views of the same case. A `deck_page` is a rendering delimiter, not
a new Julia module, Bonito session, or server route.

Pages inside the deck have an ID unique within that deck, a title, a `render`
flag, an optional CSS class, and a `build(session, state)` function. Intra-deck
addresses use URL fragments and change pages without reloading the route. A
cross-deck navigation loads the other deck's route and creates that deck's own
session.

The `group` is only a sidebar category. Multiple independent decks can use
`group = "Application cases"`; each will appear under that heading and retain
its own route, setup function, and state.

Both deck-level and page-level `render = false` flags are read during startup.
Disabled decks receive no route and are never built. Disabled pages remain in
source but do not enter navigation and are not built. Restart the server after
adding files or changing visibility.

The README contains a complete minimal two-page deck example. New files should
follow the numeric naming convention, for example
`pages/050_another_application_case.jl`, and should keep all case-specific
helpers private inside their module.

## Current decks

### Overview

`pages/010_overview.jl` is the landing deck under the `Showcase` group. It
explains the live-manual boundary and shows the real numerical-preparation
progress. It has one page, `overview`.

### Core and insulation

`pages/020_core_and_insulation.jl` is the `Cable design` deck. Its native Bonito
sliders rebuild a small `CableDesign`, derive all displayed values from that
object, and update two persistent Makie circle plots in place. Data-model
construction and tree inspection remain locally contained because the main
data model is still being refactored. It has one page,
`core-and-insulation`.

### Parametric OHL/UGC transition

`pages/030_ohl_ugc_transition.jl` is the first `Application cases` deck. Its
three pages share one setup state:

1. `corridor` shows the network diagram and the OHL/UGC composition control.
2. `b5-impedance` shows the live driving-point impedance at B5.
3. `linearization` explains the frozen operating point and provides the
   explicit recache control.

Startup preparation occurs in an isolated Julia process and reports real
network, power-flow, linearization, and response phases to the landing page. A
prepared reference network/model pair is kept process-wide. Each browser deck
session receives a deep copy, then creates one network diagram and one
impedance figure reused by all three internal pages.

Moving the underground-cable slider mutates only the OHL and UGC corridor
lengths in that session and recomputes the frequency response. It does not
rebuild the network, rerun power flow, or re-linearize. Completed positions are
cached, and stale intermediate slider requests are discarded. The B4-B5 line
in the single-line diagram changes continuously from red for OHL to blue for
UGC. `NetworkState` is intentionally not observable.

The **Re-cache power flow + linearization** button forcefully solves and
linearizes the current session at the selected split, clears that session's
response cache, and refreshes the curve. It does not alter other open sessions.

PowerImpedance is pinned from the public GitHub repository at commit
`457ba27831a3841c97e14b9c832840390df946f4`. The package's not-yet-public
GraphMakie functionality from MR 36 commit
`71183fe1251cd178bc6c8704594d197d4d988414` was copied narrowly into
`PowerImpedanceDiagramExt.jl` and its `projection.jl`/`render.jl` helpers. This
local extraction contains the topology projection and direct Makie renderer,
not the upstream PlotBuilder integration, and it modifies no package source.
Delete the local extraction once equivalent public PowerImpedance functionality
is available.

### Live Julia workspace

`pages/040_live_julia.jl` is a `Developer tools` deck with one page. It combines
Bonito's native `CodeEditor` and `TerminalOutput` with a stateful Julia worker
implemented in `repl_worker.jl`.

Run or `Ctrl+Enter` evaluates code. Definitions, imports, and `ans` persist in
that browser session. Clearing affects only the transcript. Restarting replaces
the worker and namespace. Interrupt forcefully terminates the current worker
and starts a new one. Closing the Bonito session also terminates the worker.

The worker is process-isolated from the web server but is not an operating-
system sandbox. Entered code has the same account and filesystem permissions as
the server. Never expose this route through a public or organization-visible
port without adding a real container or VM security boundary.

## Shell behavior and author-visible facilities

- The sidebar is generated from deck groups, deck titles, and internal page
  titles. Search filters both deck and page entries.
- The top bar has previous/next navigation, global page count, fullscreen, and
  browser-print controls.
- Desktop presentation mode hides the sidebar and persists for the browser tab,
  including across cross-deck route loads. Fullscreen also survives intra-deck
  page changes because those use fragments rather than route reloads.
- `Left`/`Right`, `PageUp`/`PageDown`, and horizontal touch swipes navigate.
  Events from form controls, links, editors, and canvases are ignored.
- `Ctrl+/` or `Cmd+/` focuses search; `Escape` clears it.
- Print/PDF uses the browser print dialog. The print stylesheet emits each
  enabled page of the current deck on a separate sheet.
- The top-right LED has a reserved layout slot: green is connected, yellow is
  connecting, and red is disconnected.
- The landing page remains a real interactive Bonito page during numerical
  preparation. Its progress bar is determinate and phase-driven, not an
  indeterminate animation.

Julia Markdown is the prose source. Use double backticks for inline LaTeX and a
fenced `math` block for display mathematics. Bonito's bundled KaTeX is the only
math renderer; no Reveal math plugin or second DOM pass is present. Normal Julia
string escaping still applies.

## Running locally

From the repository root:

```sh
julia --project=showcase/bonito -e 'using Pkg; Pkg.instantiate()'
julia --project=showcase/bonito showcase/bonito/serve.jl
```

The one-shot Linux launcher starts the server and opens the default browser:

```sh
./showcase/bonito/launch.sh
```

It accepts a direct deck/page target:

```sh
./showcase/bonito/launch.sh '/parametric-ohl-ugc-transition#b5-impedance'
```

The runtime interface is `HOST` (default `127.0.0.1`), `PORT` (default `8080`),
and `PROXY_URL` (default `.`). A proxy must forward HTTP and WebSocket traffic.

Focused validation and the Cairo smoke export are:

```sh
julia --project=showcase/bonito showcase/bonito/test/runtests.jl
julia --project=showcase/bonito showcase/bonito/export.jl
```

## GitHub Codespaces

`.devcontainer/Dockerfile` installs Julia 1.12.7 on the standard Debian
devcontainer base. `.devcontainer/devcontainer.json` requests four CPUs, 8 GB
RAM, 32 GB storage, installs the Julia VS Code extension, instantiates and
precompiles the showcase after creation, and forwards port 8080 as **Bonito
showcase**.

In a codespace, run:

```sh
HOST=0.0.0.0 PORT=8080 PROXY_URL=. julia --project=showcase/bonito showcase/bonito/serve.jl
```

Open the forwarded port from the Ports panel and keep it **Private**, especially
because the live Julia page executes arbitrary code. PowerImpedance comes from
public GitHub over HTTPS and LineCableModels comes from the checkout, so no
GitLab SSH key or repository secret is required.

`launch.sh` is intended for a graphical Linux desktop. A codespace should use
the forwarded-port workflow instead of expecting `xdg-open` inside the
container to open the local browser.

## Validation state and known constraints

The focused showcase tests cover deck discovery and validation, visibility
flags, shared setup state, generated routes and DOM markers, hash navigation,
responsive/fullscreen shell hooks, math ownership, the cable data-model
boundary, persistent Makie plot identity, the local PowerImpedance projection,
response caching, explicit recaching, and the isolated Julia worker. The Cairo
command exercises saving an SVG artifact.

At the last handoff validation, the focused showcase tests, Cairo export, and a
server startup smoke test passed. Re-run the two commands above after future
changes. The repository-wide `Pkg.test()` was separately blocked before this
handoff because the root `Project.toml` declared Gmsh while the root
`Manifest.toml` did not contain it; that environment issue belongs to the
parallel engine work and was not repaired from this branch.

No browser automation or visual-regression harness exists. After material UI
changes, manually verify:

- direct links to every route and page fragment;
- slider dragging and repeated updates;
- two-tab session isolation;
- intra-deck and cross-deck previous/next navigation;
- sidebar search, keyboard navigation, and touch navigation;
- fullscreen/presentation persistence;
- moving the browser between differently sized monitors before maximizing;
- WGLMakie canvas resizing without pointer-coordinate drift;
- browser printing/PDF;
- numerical preparation and explicit recaching; and
- live Julia run, interrupt, restart, and session cleanup.

The Lato and JuliaMono font files are the only remote presentation assets. They
use the same pinned CDN versions as the generated documentation. If the browser
cannot reach the CDN, it will use fallback fonts and the precise visual metrics
will differ.

## Recommended continuation order

1. Confirm `git branch --show-current` prints
   `feat/codespaces-showcase` before touching the tree.
2. Read this report, `showcase/bonito/README.md`, and the deck file closest to
   the new application case.
3. Preserve all unrelated engine/formulation changes already present in the
   working tree.
4. Add presentation content as pages inside an existing deck when it shares
   state, or add another numbered deck file when it needs an independent
   module/session. Reuse `group = "Application cases"` for additional cases.
5. Keep expensive setup in the deck's `setup(session)` boundary and keep page
   builders focused on composing views from shared state.
6. Extend `PageAuthoring` only when a genuinely reusable authoring primitive is
   missing. Do not bypass it merely for ordinary page composition.
7. Run focused tests and the Cairo smoke export, then perform the relevant
   manual browser checks before committing.
