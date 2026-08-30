# Bonito cable-design showcase

This standalone prototype presents a small live `CableDesign` as a technical
manual. A Documenter-inspired shell owns the typography, colors, sidebar, and
responsive layout. Bonito serves one application per deck route, and each deck
keeps one browser session while its hash-addressed pages change in place. The
shell adds keyboard/touch navigation, fullscreen, and browser printing. Bonito
and WGLMakie provide the live controls and figures; there is no page-owning
presentation framework.

The renderer is intentionally a showcase-local toy. It snapshots the current
material-preview colors and light Makie canvas, but it neither calls nor extends
LineCableModels' PlotBuilder machinery. Cairo parity is not a showcase goal.

The application-case page adapts PowerImpedance's `P2P_HVDC_ALT.jl` example.
PowerImpedance is pinned over public HTTPS to GitHub commit
`457ba27831a3841c97e14b9c832840390df946f4`. Until its GraphMakie work is
merged publicly, `PowerImpedanceDiagramExt.jl` carries only the topology
projection and direct Makie renderer from MR 36 commit
`71183fe1251cd178bc6c8704594d197d4d988414`; its PlotBuilder integration is
intentionally excluded. The impedance response likewise uses a showcase-local
Makie axis. No PowerImpedance package source is modified, and the local module
can be deleted when the upstream extension becomes public.

## GitHub Codespaces

The repository's `.devcontainer` configuration provides Julia 1.12 and prepares
the pinned showcase environment when a codespace is created. In GitHub, choose
**Code → Codespaces → New with options**, select the branch containing the
configuration, and create the codespace. The requested minimum machine is four
cores, 8 GB of memory, and 32 GB of storage.

The first creation downloads and precompiles the Julia dependencies. When the
post-create command finishes, start the application from the codespace terminal:

```sh
HOST=0.0.0.0 PORT=8080 PROXY_URL=. julia --project=showcase/bonito showcase/bonito/serve.jl
```

Open **Bonito showcase** from the VS Code **Ports** panel when port 8080 is
reported. Keep its visibility set to **Private**. The Live Julia workspace
evaluates arbitrary Julia code with the codespace user's permissions and must
not be exposed through an organization-visible or public port.

No GitLab credential or SSH secret is required: the showcase pins
PowerImpedance from its public GitHub repository over HTTPS, and
LineCableModels is loaded from the codespace checkout.

Instantiate the pinned application environment from the repository root:

```sh
julia --project=showcase/bonito -e 'using Pkg; Pkg.instantiate()'
```

## Deck discovery and authoring

Every `pages/*.jl` file is one self-contained deck module. At startup, `app.jl`
includes the files, validates their final `DECK` descriptor, sorts them by
`order` and filename, and generates the sidebar, deck routes, and internal page
navigation. A deck owns its Bonito widgets, observables, Makie figures,
callbacks, and helper methods. Its `setup(session)` function runs exactly once
for that browser session and returns the state shared by every page in the
file.

`PageAuthoring.jl` keeps presentation DOM construction out of deck files. Its
small public vocabulary is `slide`, `article`, `prose`, `webgrid`, `webpart`,
`control`, `value_list`, `diagnostic`, `status_line`, `color_key`,
`local_image`, `deck_page`, and `deck_descriptor`. Prose uses Julia's native
`md"""..."""` literals, which can interpolate live Bonito widgets,
Observables, and Makie content directly.

```julia
module MyDeck

using Bonito
using Markdown
using ..PageAuthoring

function setup(session::Session)
    slider = Bonito.Slider(1:10; value = 5)
    doubled = map(value -> 2value, session, slider.value)
    return (; slider, doubled)
end

function controls_page(::Session, state)
    canvas = webgrid(
        [:description :figure; :controls :figure];
        columns = "2fr 5fr",
        rows = "auto minmax(0, 1fr)",
        description = webpart(prose(md"""
        Explain the application here. Julia Markdown supports **formatting**,
        tables, images, code, and LaTeX such as ``Z(j\\omega)``.
        """)),
        controls = webpart(
            control("Example value", state.slider; value = state.doubled);
            kind = :controls
        ),
        figure = webpart("A WGLMakie figure goes here"; kind = :plot)
    )
    return (; body = slide("Controls", canvas; lede = md"An introduction."))
end

function interpretation_page(::Session, state)
    return (;
        body = slide(
            "Interpretation",
            article(md"The same deck state is still live: **$(state.doubled)**.")
        )
    )
end

const DECK = deck_descriptor(
    id = "my-deck",
    group = "Examples",
    title = "My deck",
    order = 40,
    render = true,
    setup = setup,
    pages = (
        deck_page(
            "Controls";
            id = "controls",
            class = "lc-my-controls-page",
            build = controls_page
        ),
        deck_page(
            "Interpretation";
            id = "interpretation",
            render = true,
            build = interpretation_page
        )
    )
)

end


MyDeck.DECK
```

Each `deck_page(...)` call is an explicit rendering delimiter inside the Julia
file. It changes only how the deck is divided on screen; it does not create a
new module, route, or Bonito session. A page builder receives `(session, state)`
and returns a named tuple containing `body`. Keep resource-heavy initialization
in `setup`, and let page builders compose views from the returned state.

The `webgrid` matrix is the canvas layout. Repeating an area name merges those
cells; use `nothing` for an empty cell. The named keyword arguments must match
the matrix areas, so malformed layouts fail during startup rather than becoming
broken CSS in the browser. On narrow screens, grids stack their named areas in
reading order unless `stack = false` is requested.

Set `render = false` on either a deck or an individual `deck_page` to leave it
out of that run. A disabled deck receives no route and is never constructed. A
disabled internal page does not enter navigation and its builder is never
called. Adding a file or changing either flag takes effect after restarting the
server; there is no runtime watcher. All `.jl` files in `pages/` are decks, so
keep private helpers inside their module. Deck IDs must be globally unique
lowercase hyphenated strings; page IDs must be unique within their deck.

Start the live manual:

```sh
julia --project=showcase/bonito showcase/bonito/serve.jl
```

On a Linux desktop, start the server and open the default `xdg` browser at the
landing page with one command:

```sh
./showcase/bonito/launch.sh
```

An optional first argument opens a particular deck page directly:

```sh
./showcase/bonito/launch.sh '/parametric-ohl-ugc-transition#b5-impedance'
```

Without the launcher, open <http://127.0.0.1:8080/#overview>. Enabled decks have
ordinary server routes: `/`, `/core-and-insulation`,
`/parametric-ohl-ugc-transition`, and `/live-julia-workspace`. Pages inside a
deck use URL fragments, such as
`/parametric-ohl-ugc-transition#linearization`. Fragment navigation does not
reload the document or create a new Bonito session, so deck state, WGL figures,
and native browser fullscreen survive the transition. Moving to another deck
loads that deck's route and creates its own session.

The top bar contains previous/next links across all enabled pages, the current
page number, fullscreen, and print/PDF controls. The menu button collapses the
sidebar on a desktop and opens it as a drawer on narrow screens. Presentation
mode hides the sidebar and lets the page occupy the complete browser surface.
That mode is stored for the current browser tab and is restored after a
cross-deck route load.

Typography scales gently with the viewport: the existing 1600-pixel-wide look
remains the reference size, smaller screens stop at 15 px, and large displays
grow to an 18 px base. This is ordinary scoped CSS rather than a presentation
transform, so page geometry and pointer coordinates remain native browser
coordinates. Browser zoom (`Ctrl`/`Cmd` with `+` or `-`, and `0` to reset) is
the explicit override when a room or projector needs still larger text.

The shell reconciles its exact CSS-pixel viewport immediately, after two
animation frames, and once more after the browser/window manager settles. It
listens to layout-viewport, visual-viewport, and display-resolution changes,
which covers moving an already-open window between monitors before maximizing.
WGLMakie remains configured through its public `resize_to = :parent` interface;
the shell does not inspect or modify WGLMakie canvas internals.

The branded landing shell is a complete Bonito page and remains interactive
while the numerical application case is prepared in a separate Julia process.
Its determinate progress bar reports the real network, power-flow,
linearization, and impedance-response phases. Opening the application-case
route before preparation finishes shows the same status and reloads that route
once the prepared model is available. The deck then constructs its diagram and
impedance canvas once for that session and reuses them across its internal
pages. The small status LED occupies a reserved slot at the right of the top
bar: green is connected, yellow is connecting, and red is disconnected.

Use the sidebar, `Left`/`Right`, `PageUp`/`PageDown`, or a horizontal touch
swipe to navigate. Keyboard navigation deliberately ignores controls, links,
and canvases. Search filters deck and page titles. Press `Ctrl+/` (or `Cmd+/`)
to focus search and `Escape` to clear it. Navigation uses ordinary links, so
every page has a stable, directly loadable route-and-fragment address.

The **Live Julia workspace** uses Bonito's native Julia `CodeEditor` and
`TerminalOutput`. `Run` or `Ctrl+Enter` evaluates the editor in a Julia worker
process owned by that browser session. Definitions, imports, and `ans` persist
between runs. `Clear console` changes only the transcript; `Restart worker`
discards the execution namespace. `Interrupt` forcefully terminates a running
evaluation and starts a fresh worker, so it also clears that worker's state.
Closing the Bonito session terminates its worker.

The worker is process-isolated from the Bonito server, not container-sandboxed.
Code entered in the page runs with the same operating-system account and file
permissions as the server. Keep this page on the default localhost binding; do
not expose it to a LAN or public interface unless the worker is moved behind an
appropriate container or virtual-machine boundary.

Julia-authored prose uses Julia Markdown's LaTeX syntax: double backticks for
inline expressions, such as ``` ``Z(j\omega)`` ```, and a fenced `math` block
for display expressions. Bonito renders those nodes with its bundled KaTeX.
Every `prose` region is marked `lc-bonito-markdown` and records Bonito as its
renderer. There is no second mathematics pass, so Julia string escaping is the
only syntax boundary authors need to consider.

In the **Application case - parametric OHL/UGC transition** deck, move the slider to
change the underground-cable share. The B4–B5 corridor changes continuously
from red (OHL) to blue (UGC), and the existing B5 impedance line is updated in
place. Server startup prepares one process-wide reference network, power flow,
linearization, and initial response in the isolated numerical worker. Each
application-case deck session receives its own deep copy of the prepared
network/model pair and creates exactly one network diagram and one impedance
figure shared by the deck's three pages. During
slider motion, the showcase mutates only that session's OHL and UGC lengths and
recomputes the frequency response; it does not rebuild the network, rerun the
power flow, or re-linearize. Completed slider positions are cached and
intermediate requests are dropped while a newer value is waiting.

The operating point is intentionally frozen at the default 50% split. Use
**Re-cache power flow + linearization** to force a fresh reference solve, rebuild
the current session's linearization at the selected split, clear its response
cache, and refresh the curve. Other already-open sessions retain their own
linearized models. `NetworkState` itself is not observable.

The launcher reads `HOST` (default `127.0.0.1`), `PORT` (default `8080`), and
`PROXY_URL` (default `.`). A relative proxy URL requires the proxy to forward
both HTTP and WebSocket traffic.

Use the print icon in the top bar to open the browser print dialog; choose
**Save as PDF** for a PDF artifact. Printing acts on the current deck route;
the print stylesheet removes the sidebar and navbar and prints each enabled
page in that deck on a separate sheet.

An existing CairoMakie smoke helper remains available, but its output is not a
visual acceptance target for this application:

```sh
julia --project=showcase/bonito showcase/bonito/export.jl
```

Run the focused tests:

```sh
julia --project=showcase/bonito showcase/bonito/test/runtests.jl
```

The visual theme belongs exclusively to `assets/theme.css`, which contains a
small scoped snapshot of the repository's Documenter dark theme. Its contents
are embedded in each page so the shell does not depend on Bonito's local asset
registration or the browser's asset cache. Lato and JuliaMono are loaded from
the same pinned CDN versions as the generated documentation; these fonts are
the only remote presentation assets.
