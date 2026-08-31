# LineCableModels playground

The playground is a small Bonito-native dynamic publisher. It provides the
Documenter-inspired shell, responsive layout, navigation, fullscreen, printing,
connection state, and authoring primitives. The publisher does not decide
whether its content is a live manual, a targeted tutorial, an application
study, or a presentation assembled for one audience.

Presentation content lives below `presentations/<slug>/content`. Each Julia
file is one stateful deck, and each deck may contain several hash-addressed
pages that share one Bonito session. No page-owning JavaScript presentation
framework is used.

## Run it

Instantiate the pinned Julia 1.12 environment from the repository root:

```sh
julia --project=showcase/bonito -e 'using Pkg; Pkg.instantiate()'
```

Start the server:

```sh
julia --project=showcase/bonito showcase/bonito/serve.jl
```

On a graphical Linux desktop, start the server and open the default browser in
one command:

```sh
./showcase/bonito/launch.sh
```

An optional argument opens a presentation page directly:

```sh
./showcase/bonito/launch.sh \
  '/presentations/live-manual/parametric-ohl-ugc-transition#b5-impedance'
```

The launcher reads `HOST` (default `127.0.0.1`), `PORT` (default `8080`), and
`PROXY_URL` (default `.`). A proxy must forward HTTP and WebSocket traffic.

## Publisher and content boundary

The application-owned files are:

```text
showcase/bonito/
├── app.jl
├── home.jl
├── PageAuthoring.jl
├── serve.jl
├── launch.sh
├── assets/theme.css
└── presentations/
```

`home.jl` is a first-class Bonito landing page, not a deck. It introduces the
playground, displays the project logo, discovers publication folders, and
provides the presentation dropdown, explicit **Activate** step, bounded
preflight console, real progress, and gated **Open** button. Home is always the
first sidebar entry.

`app.jl` discovers URL-safe, case-sensitive directories containing `content/`.
Lowercase hyphenated names remain the default convention; uppercase letters are
supported when a presentation slug is an acronym such as `ICHQP2026`.
Each presentation is loaded into its own Julia namespace, so separate
presentations may reuse module names and deck IDs without colliding. The
publisher validates the descriptors, builds namespaced routes, and displays
only the active presentation in the sidebar.

Routes have this stable form:

```text
/
/presentations/<slug>/<deck-id>#<page-id>
```

The application includes deck modules during startup, but it does not run their
`setup` functions or heavy preparation. A deck session is constructed only
when its route is opened. Only resources declared with `preflight_resource`
participate in activation. Opening a cold heavy route directly invokes the
same idempotent activation as a defensive fallback.

## Create a presentation

Copy the executable starter:

```sh
cp -a showcase/bonito/presentations/_template \
  showcase/bonito/presentations/my-tutorial
```

The target directory name is the presentation slug. It must contain URL-safe
words separated by hyphens. Its display name is derived mechanically:
`my-tutorial` becomes **My tutorial**, while `ICHQP2026` remains **ICHQP2026**.
Deck and page identifiers remain lowercase and hyphenated. There is no
presentation manifest, generator, database, or configuration framework.

`_template` is deliberately ignored by normal discovery because its name is
not a valid slug. Restart the server after copying or renaming a presentation.

The template contains one executable multi-page deck demonstrating:

- Julia Markdown, code, tables, inline math, and display math;
- sliders, dropdowns, buttons, mapped values, and action callbacks;
- persistent WGLMakie figures updated through observables;
- an opt-in process-wide preflight resource;
- presentation-local images and captions; and
- all predefined page layouts.

The focused tests load the ignored template explicitly, so it cannot silently
drift away from the supported authoring API.

## Write a deck

Every `content/*.jl` file defines one module and returns one `DECK` descriptor:

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
    controls = webpart(
        control("Example value", state.slider; value = state.doubled);
        kind = :controls,
    )
    explanation = webpart(
        prose(md"The same state remains available to every page.");
        kind = :panel,
    )
    return (;
        body = slide(
            "Controls",
            layout_two_columns(
                controls,
                explanation;
                ratio = :narrow_wide,
            ),
        ),
    )
end

function interpretation_page(::Session, state)
    return (;
        body = slide(
            "Interpretation",
            layout_single(article(md"Result: **$(state.doubled)**.")),
        ),
    )
end

const DECK = deck_descriptor(
    id = "my-deck",
    group = "Examples",
    title = "My deck",
    order = 10,
    render = true,
    setup = setup,
    pages = (
        deck_page("Controls"; id = "controls", build = controls_page),
        deck_page(
            "Interpretation";
            id = "interpretation",
            build = interpretation_page,
        ),
    ),
)

end

MyDeck.DECK
```

`setup(session)` runs once for the browser's deck session. Keep shared widgets,
observables, Makie figures, callbacks, caches, and heavy models there. Each page
builder receives `(session, state)` and returns a named tuple containing `body`.
A `deck_page` is only a rendering delimiter; changing the page fragment does
not construct another module, route, or session.

Set `render = false` on a deck or page to exclude it after the next server
restart. Deck IDs must be unique inside their presentation. Page IDs must be
unique inside their deck. Number deck files consistently with their `order`,
for example `020_application_case.jl` with `order = 20`.

## Presentation preflight

Preflight is for expensive, process-wide resources that should be ready before
an audience enters a presentation. It does not construct deck sessions,
widgets, figures, or REPL workers.

Declare state and an idempotent activation function in the deck that owns the
resource:

```julia
const MODEL_STATE = preflight_state(label = "Model not activated")

function activate_model!()
    MODEL_STATE[].phase ∉ (:cold, :error) && return nothing
    set_preflight!(MODEL_STATE, :preparing, 0.1, "Loading inputs…")
    return @async begin
        try
            model = build_expensive_model()
            MODEL_CACHE[] = model
            set_preflight!(MODEL_STATE, :hot, 1.0, "Model hot")
        catch error
            set_preflight!(
                MODEL_STATE,
                :error,
                MODEL_STATE[].progress,
                sprint(showerror, error),
            )
        end
    end
end

const MODEL_RESOURCE = preflight_resource(
    id = "network-model",
    title = "Network model",
    state = MODEL_STATE,
    activate = activate_model!,
)
```

Attach it to the descriptor:

```julia
const DECK = deck_descriptor(
    # ...
    resources = (MODEL_RESOURCE,),
)
```

The publisher aggregates required resources from the selected presentation.
`Activate` is disabled while preparation runs, the console receives explicit
state transitions, and `Open` is enabled only when the aggregate phase is
`hot`. State and caches are server-process-wide, so concurrent browsers reuse
one activation. Restarting Julia returns resources to `cold`.

## Master layouts

`PageAuthoring.jl` owns the stable layout functions. They are ordinary Julia
composition helpers, not a second presentation language:

```julia
layout_single(content)
layout_two_columns(left, right; ratio = :equal)
layout_three_columns(left, center, right)
layout_top_split(top, left, right; ratio = :equal)
layout_split_bottom(left, right, bottom; ratio = :equal)
layout_sidebar_main(sidebar, main)
layout_controls_plot(controls, plot; footer = nothing)
layout_two_rows(top, bottom; ratio = :equal)
layout_quad(top_left, top_right, bottom_left, bottom_right)
```

Horizontal ratios are `:equal`, `:narrow_wide`, and `:wide_narrow`. Vertical
ratios are `:equal`, `:short_tall`, and `:tall_short`. Every master supplies
validated regions, common gaps, narrow-screen stacking, and print behavior.

For a page that must consume all space below its heading, give the page
descriptor `class = "lc-fill-page"` and use `height = "100%"` on its master
layout. Add `lc-fluid-type` when prose should scale with the available page
width and height. Both utilities retain bounded type sizes and fall back to a
scrolling document layout on compact screens.

Use `webgrid` directly only for a composition that genuinely falls outside the
masters. Its matrix repeats an area name to merge cells:

```julia
webgrid(
    [:summary :summary; :controls :figure];
    rows = "auto minmax(0, 1fr)",
    summary,
    controls,
    figure,
)
```

The other authoring primitives are `slide`, `title_canvas`, `article`, `prose`, `webpart`,
`control`, `value_list`, `diagnostic`, `status_line`, `color_key`, and
`local_image`.

Julia Markdown uses double backticks for inline LaTeX, such as
``` ``Z(j\omega)`` ```, and a fenced `math` block for display expressions.
Bonito's bundled KaTeX is the only math renderer.

## Current live-manual presentation

The existing showcase material was migrated without changing its state model:

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

The cable deck remains a small presentation-local renderer and does not call or
extend `LineCableModels.PlotBuilder`. The ongoing data-model coupling remains
confined to its construction and geometry functions.

The OHL/UGC deck shares one network/model pair, network diagram, and impedance
figure across three pages. Slider motion changes only the OHL/UGC corridor
lengths and recomputes the cached frequency response. It does not rebuild the
network, rerun power flow, or re-linearize. The explicit recache button rebuilds
the current session's operating point.

PowerImpedance is pinned to public GitHub commit
`457ba27831a3841c97e14b9c832840390df946f4`. Until equivalent diagram support
is public upstream, `support/PowerImpedanceDiagramExt.jl` contains only the
topology projection and direct Makie renderer extracted from MR 36 commit
`71183fe1251cd178bc6c8704594d197d4d988414`.

The Live Julia deck uses Bonito's native `CodeEditor` and `TerminalOutput` with
a stateful worker process. It is process-isolated from the Bonito server but not
container-sandboxed. Entered code has the server account's filesystem
permissions. Never expose that route through a public port without a real
container or VM boundary.

## Shell behavior

The sidebar, search, previous/next navigation, global page count, keyboard and
touch navigation, presentation mode, browser fullscreen, printing, and
connection LED are publisher-owned. Presentation mode is stored for the
browser tab. Page changes within a deck use fragments and preserve the current
session and fullscreen state.

Typography scales within bounded viewport limits. Browser zoom remains the
explicit projector override. The shell reconciles normal, visual, and display-
resolution viewport changes; WGLMakie continues to use only its public
`resize_to = :parent` behavior.

Use the print button and select **Save as PDF** in the browser. The print
stylesheet emits each enabled page in the current deck on a separate sheet.

## GitHub Codespaces

The repository's `.devcontainer` installs Julia 1.12, instantiates this
environment, and forwards port 8080. In a codespace run:

```sh
HOST=0.0.0.0 PORT=8080 PROXY_URL=. \
  julia --project=showcase/bonito showcase/bonito/serve.jl
```

Open **Bonito showcase** from the Ports panel and keep the port **Private**.
No GitLab key is required: PowerImpedance is fetched from public GitHub and
LineCableModels is loaded from the checkout.

## Validation

Run the focused suite:

```sh
julia --project=showcase/bonito showcase/bonito/test/runtests.jl
```

The legacy Cairo smoke export for the cable deck remains available:

```sh
julia --project=showcase/bonito showcase/bonito/export.jl
```

The visual system belongs to `assets/theme.css`. Do not add Reveal.js, another
page-owning framework, or selectors into Bonito/WGLMakie internals.
