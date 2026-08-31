# LineCableModels playground continuation handoff

This report records the showcase architecture on `feat/codespaces-showcase`
after the publisher/content split. The complete authoring guide is
`showcase/bonito/README.md`.

## Worktree and branch

Showcase development belongs in the dedicated linked worktree:

```text
/home/amartins/Documents/KUL/LineCableModels-codespaces-showcase
```

It must remain on `feat/codespaces-showcase`. The original checkout at
`/home/amartins/Documents/KUL/LineCableModels` is reserved for release work.
The feature worktree's root `AGENTS.md` enforces the branch check.

## Product boundary

The application is now **LineCableModels playground**, a Bonito-native dynamic
publisher. It owns the Documenter-inspired shell, Home, presentation discovery,
namespaced routes, responsive layout, navigation, fullscreen, browser printing,
connection state, and authoring primitives.

It is not itself a live manual. A live manual, tutorial, study, or audience-
specific presentation is content below `presentations/<slug>/content`.
Reveal.js was removed permanently because its viewport and navigation ownership
conflicted with Bonito/WGLMakie and the accepted shell.

The visual baseline in `assets/theme.css` is accepted. Keep the dark sidebar,
top bar, typography, spacing, light Makie canvases, scoped CSS, native browser
coordinates, and public `resize_to = :parent` WGLMakie behavior.

## Publisher architecture

- `app.jl` discovers presentations, isolates their Julia namespaces, validates
  decks, creates namespaced routes, and owns shell behavior.
- `home.jl` is a first-class application page. It contains introductory text,
  the logo, discovered-presentation dropdown, explicit **Activate** step,
  bounded activity console, real progress, and gated **Open** button. Home is
  always the first sidebar entry.
- `PageAuthoring.jl` owns prose, webparts, controls, deck descriptors, raw
  `webgrid`, and the predefined master layouts.
- `serve.jl` activates WGLMakie and mounts `playground_routes()`.
- `launch.sh` starts the server and opens Home through `xdg-open`.
- `assets/theme.css` owns all presentation styling.

Routes are:

```text
/
/presentations/<slug>/<deck-id>#<page-id>
```

Presentation folder names are the only publication metadata. They must be
lowercase hyphenated slugs and are converted mechanically to display names.
There is deliberately no manifest, generator, database, or configuration DSL.

Each presentation is included inside a dedicated Julia namespace. Module names
and deck IDs therefore need to be unique only within one presentation. Module
definitions load during discovery, while deck `setup` remains route-lazy.
Expensive process-wide preparation is opt-in through `preflight_resource` and
runs only after Home activation or a defensive direct-route fallback.

Preflight resources expose independent observable `cold`, `preparing`, `hot`,
or `error` state. The publisher aggregates required resources without knowing
their domain. The current OHL/UGC deck declares the reference power flow and
linearization as its sole resource. Cable state, WGLMakie figures, and the Live
Julia worker are not preflight resources.

## Deck model

Every `content/*.jl` file is one module and one deck. Its final expression
returns a `DECK` descriptor. A deck may contain many `deck_page` rendering
delimiters. `setup(session)` runs once for the browser's deck session, and its
state is passed to every page builder in the file.

Page changes within a deck use URL fragments, preserving widgets, observables,
figures, caches, callbacks, fullscreen, and session identity. Cross-deck
navigation creates the target deck's session.

Deck-level and page-level `render = false` flags take effect after a server
restart. Deck IDs are unique within a presentation; page IDs are unique within
a deck.

## Master layouts and starter

`PageAuthoring` exports:

```julia
layout_single
layout_two_columns
layout_three_columns
layout_top_split
layout_split_bottom
layout_sidebar_main
layout_controls_plot
layout_two_rows
layout_quad
```

They provide named, validated regions with standard spacing, responsive
stacking, and print behavior. Horizontal semantic ratios are `:equal`,
`:narrow_wide`, and `:wide_narrow`; vertical ratios are `:equal`,
`:short_tall`, and `:tall_short`. `webgrid` remains the custom-layout escape
hatch.

`presentations/_template` is excluded from discovery by its invalid underscore
slug but loaded explicitly by tests. Its single executable deck demonstrates
Markdown, math, code, tables, sliders, dropdowns, buttons, observables, button
callbacks, a persistent WGLMakie figure, local images, and all nine layouts.

Create a presentation with:

```sh
cp -a showcase/bonito/presentations/_template \
  showcase/bonito/presentations/my-tutorial
```

Then rename the module/descriptor fields, delete unwanted examples, write the
content, and restart the server.

## Current live-manual content

The original showcase now lives under:

```text
presentations/live-manual/
├── content/
│   ├── 010_core_and_insulation.jl
│   ├── 020_ohl_ugc_transition.jl
│   └── 030_live_julia.jl
└── support/
    ├── PowerImpedanceDiagramExt.jl
    ├── PowerImpedanceDiagramExt/
    └── repl_worker.jl
```

The cable renderer is presentation-local and intentionally does not integrate
with `LineCableModels.PlotBuilder`. Data-model construction/tree inspection
remain confined to the cable deck because the main data model is still moving.

The OHL/UGC deck shares one session-owned network/model pair, diagram, and
impedance figure across three pages. Slider motion changes only OHL/UGC lengths
and recomputes a cached response. Power flow and linearization change only
through the explicit recache action. The not-yet-public PowerImpedance diagram
code remains narrowly extracted in this presentation's support directory.

The Live Julia deck uses Bonito's native editor and terminal with a stateful
worker process. It is not an OS sandbox. Keep the app local or the Codespaces
port private.

## Commands

```sh
julia --project=showcase/bonito -e 'using Pkg; Pkg.instantiate()'
./showcase/bonito/launch.sh
julia --project=showcase/bonito showcase/bonito/test/runtests.jl
julia --project=showcase/bonito showcase/bonito/export.jl
```

Direct launch example:

```sh
./showcase/bonito/launch.sh \
  '/presentations/live-manual/parametric-ohl-ugc-transition#b5-impedance'
```

For Codespaces:

```sh
HOST=0.0.0.0 PORT=8080 PROXY_URL=. \
  julia --project=showcase/bonito showcase/bonito/serve.jl
```

## Constraints and continuation

- Do not restore Reveal.js or add another page-owning frontend framework.
- Do not move PlotBuilder behavior into the playground.
- Keep content-specific helpers under their presentation, not beside `app.jl`.
- Keep session state in deck `setup`; keep only deliberately process-wide
  caches behind content-owned preflight resources. Home may orchestrate those
  resources but must not construct their domain objects itself.
- Extend `PageAuthoring` only for a genuinely reusable authoring pattern.
- Keep the shared Julia environment for now; per-presentation environments are
  intentionally out of scope.
- Run the focused suite after publisher or template changes and manually verify
  presentation discovery, direct routes, fullscreen persistence, monitor
  resizing, print/PDF, slider behavior, two-tab isolation, and worker cleanup.
