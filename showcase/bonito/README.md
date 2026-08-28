# Bonito cable-design showcase

This standalone prototype presents a small live `CableDesign` as a technical
manual. A Documenter-inspired shell owns the typography, colors, sidebar, and
responsive layout; embedded reveal.js supplies slide state, keyboard/touch
navigation, overview, speaker notes, and browser printing. Bonito and WGLMakie
provide the live controls and figure.

The renderer is intentionally a showcase-local toy. It snapshots the current
material-preview colors and light Makie canvas, but it neither calls nor extends
LineCableModels' PlotBuilder machinery. Cairo parity is not a showcase goal.

Instantiate the pinned application environment from the repository root:

```sh
julia --project=showcase/bonito -e 'using Pkg; Pkg.instantiate()'
```

Start the live manual:

```sh
julia --project=showcase/bonito showcase/bonito/serve.jl
```

Open <http://127.0.0.1:8080>. Use the sidebar, focused arrow keys, or touch to
navigate. Search filters slide titles. Press `Ctrl+/` (or `Cmd+/`) while the app
has focus to jump to search, and press `Escape` to clear it. Reveal's overview
and speaker-note shortcuts remain available when the content pane is focused.

The launcher reads `HOST` (default `127.0.0.1`), `PORT` (default `8080`), and
`PROXY_URL` (default `.`). A relative proxy URL requires the proxy to forward
both HTTP and WebSocket traffic.

For browser/PDF printing, open
<http://127.0.0.1:8080/?view=print> and use the browser print dialog. The print
stylesheet removes the documentation sidebar and navbar from the output.

An existing CairoMakie smoke helper remains available, but its output is not a
visual acceptance target for this application:

```sh
julia --project=showcase/bonito showcase/bonito/export.jl
```

Run the focused tests:

```sh
julia --project=showcase/bonito showcase/bonito/test/runtests.jl
```

The reveal.js 6.0.1 JavaScript, structural CSS, and notes module are pinned
remote assets. The visual theme belongs exclusively to `assets/theme.css`,
which contains a small scoped snapshot of the repository's Documenter dark
theme. Lato and JuliaMono are also loaded from the same pinned CDN versions as
the generated documentation. Initial asset loading therefore requires network
access.
