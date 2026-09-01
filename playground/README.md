# LineCableModels playground

This directory is the clean publisher foundation. Quarto owns the static home
source, document structure, responsive layout, navigation, and theme. Julia
renders that source and Bonito serves the generated site. NATS is reserved for
the later broker boundary; the home page starts no numerical workers.

## Bootstrap

From the repository root, run:

```sh
./playground/bootstrap.sh
```

The idempotent bootstrap:

- requires Julia 1.12 or newer;
- downloads the pinned Quarto CLI into the user's XDG data directory;
- verifies the official archive SHA-256 before extracting it;
- instantiates the Julia environment;
- registers the developed `linecablemodels` application;
- exposes it through `~/.local/bin`; and
- renders the initial static site.

It uses no `sudo` and does not modify shell startup files. Linux x86-64 and
ARM64 are supported. Rerunning it reuses a valid installation.

If `~/.local/bin` is not already on `PATH`, the bootstrap prints the exact
export line to add. The working-tree source remains live, so editing the Julia
package does not require reinstalling the command.

## Author and build

The home is authored in `index.qmd`; `_quarto.yml` owns site configuration and
`assets/theme.scss` owns the visual treatment. Render it with:

```sh
linecablemodels playground build
```

The generated `_site/` directory is disposable and ignored by Git.

## Reuse a page template

Open `/templates/` in the running site to inspect the rendered gallery. Every
entry is backed by the corresponding file in `templates/`; its **Source**
control shows the complete authoring source. Copy the closest `.qmd` file,
change its title and contents, and retain the `lc-layout-*` classes you need.

The gallery currently covers a full canvas, balanced columns, a feature with a
sidebar, a top band over two columns, a 2×2 dashboard, media with narrative,
and the intended responsive Bonito viewport boundary.

## Reuse a live widget

Open `/widgets/` in the running site for a catalog of session-local Bonito
controls. Quarto owns the surrounding explanation and layout; each Bonito app
is mounted through a responsive, same-origin iframe:

```markdown
{{< bonito route="/widgets/control-panel" title="Cable controls" height="34rem" >}}
```

The local shortcode is implemented in `_extensions/bonito/`. Reusable widget
factories and their registered routes live in `src/widgets.jl`. The examples
cover sliders, numeric spinners, checkboxes, dropdowns, text fields, actions,
progress feedback, a read-only console sink, persistent tabs, namespaced
toolbars, and a composed engineering control panel. Layout primitives are
consumed directly from `BonitoWidgets`; the console accepts observable Julia,
worker, or broker messages and starts no execution machinery of its own.

## Start

```sh
linecablemodels playground start
```

`start` renders the site, registers the generated static files beside Bonito's
live routes, starts the server, and opens the
default browser. Configuration can come from `HOST`, `PORT`, and `PROXY_URL`,
or explicit options:

```sh
linecablemodels playground start --host 127.0.0.1 --port 8080
linecablemodels playground start --no-open
linecablemodels playground start --no-render
```

Use `--no-render` only when `_site/index.html` already exists. Press `Ctrl+C`
to stop the server cleanly.
