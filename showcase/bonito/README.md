# Bonito cable-design showcase

This standalone prototype presents a small live `CableDesign` as a technical
manual. A Documenter-inspired shell owns the typography, colors, sidebar, and
responsive layout. Bonito serves one application per real page route; the shell
adds keyboard/touch navigation, fullscreen, and browser printing. Bonito and
WGLMakie provide the live controls and figures; there is no page-owning
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

Instantiate the pinned application environment from the repository root:

```sh
julia --project=showcase/bonito -e 'using Pkg; Pkg.instantiate()'
```

## Page discovery

Every `pages/*.jl` file is one self-contained page module. At startup, `app.jl`
includes the files, validates their final `PAGE` descriptor, sorts them by
`order` and filename, and generates the sidebar and manual pages. A page owns
its Bonito widgets, observables, Makie figures, callbacks, and helper methods.

`PageAuthoring.jl` keeps presentation DOM construction out of page files. Its
small public vocabulary is `slide`, `article`, `prose`, `webgrid`, `webpart`,
`control`, `value_list`, `diagnostic`, `status_line`, `color_key`,
`local_image`, and `page_descriptor`. Page prose uses Julia's native
`md"""..."""` literals, which can interpolate live Bonito widgets,
Observables, and Makie content directly.

The visibility switch lives in the page descriptor:

```julia
module MyPage

using Bonito
using Markdown
using ..PageAuthoring

function build(::Session)
    content = webgrid(
        [:description :figure; :controls :figure];
        columns = "2fr 5fr",
        rows = "auto minmax(0, 1fr)",
        description = webpart(prose(md"""
        Explain the application here. Julia Markdown supports **formatting**,
        tables, images, code, and LaTeX such as ``Z(j\\omega)``.
        """)),
        controls = webpart("Bonito controls go here"; kind = :controls),
        figure = webpart("A WGLMakie figure goes here"; kind = :plot)
    )
    return (;
        body = slide("My page", content; lede = md"A concise introduction.")
    )
end

const PAGE = page_descriptor(
    id = "my-page",
    group = "Examples",
    title = "My page",
    order = 40,
    render = true,
    class = "lc-my-page",
    build = build
)

end


MyPage.PAGE
```

The `webgrid` matrix is the canvas layout. Repeating an area name merges those
cells; use `nothing` for an empty cell. The named keyword arguments must match
the matrix areas, so malformed layouts fail during startup rather than becoming
broken CSS in the browser. On narrow screens, grids stack their named areas in
reading order unless `stack = false` is requested.

Set `render = false` to leave a page out of that run. Disabled pages receive no
route, do not enter the navigation, and are never constructed in a browser
session. Adding a file or changing its flag takes effect after restarting the
server; there is no runtime watcher. All `.jl` files in `pages/` are treated as
pages, so keep any private helpers inside their page module rather than adding
helper files to the directory. Page IDs must be unique lowercase hyphenated
strings.

Start the live manual:

```sh
julia --project=showcase/bonito showcase/bonito/serve.jl
```

Open <http://127.0.0.1:8080>. The enabled pages have ordinary server routes:
`/`, `/core-and-insulation`, `/parametric-ohl-ugc-transition`, and
`/live-julia-workspace`. The top bar contains previous/next links, the current
page number, fullscreen, and print/PDF controls. The menu button collapses the
sidebar on a desktop and opens it as a drawer on narrow screens. Presentation
mode hides the sidebar and lets the page occupy the complete browser surface.
That mode is stored for the current browser tab, so it remains active when an
ordinary page link loads the next route. The button also requests native browser
fullscreen when available. Browsers leave native fullscreen when loading a new
document and do not permit the next document to re-enter it without another
user gesture; the persistent presentation layout therefore remains active even
when the browser chrome returns. Use the browser's own fullscreen command when
browser-level fullscreen must span route changes.

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
once the prepared model is available. Plot canvases are constructed only when
their own route is requested. The small status LED occupies a reserved slot at
the right of the top bar: green is connected, yellow is connecting, and red is
disconnected.

Use the sidebar, `Left`/`Right`, `PageUp`/`PageDown`, or a horizontal touch
swipe to navigate. Keyboard navigation deliberately ignores controls, links,
and canvases. Search filters page titles. Press `Ctrl+/` (or `Cmd+/`) to focus
search and `Escape` to clear it. Navigation uses ordinary links, so every page
has a stable, directly loadable path.

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

On **Application case - parametric OHL/UGC transition**, move the slider to
change the underground-cable share. The B4–B5 corridor changes continuously
from red (OHL) to blue (UGC), and the existing B5 impedance line is updated in
place. Server startup prepares one process-wide reference network, power flow,
linearization, and initial response in the isolated numerical worker. Each
application-case session receives its own deep copy of the
prepared network/model pair and creates only its own two Makie figures. During
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
**Save as PDF** for a PDF artifact. Printing acts on the current route; the
print stylesheet removes the sidebar and navbar and fits that page to the
sheet.

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
