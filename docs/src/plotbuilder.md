```@meta
EditURL = "../literate/plotbuilder.jl"
```

# PlotBuilder guide

`PlotBuilder` separates scientific preparation from backend drawing. A domain
owner completes detached pages; an explicitly loaded Makie extension draws
those pages inside one standard shell.

```text
source + plot definition
         │
         ▼
entitle → parse → resolve → fetch → finish
                                   │
                                   ▼
                      PlotRecipe(definition, pages)
                                   │
                                   ▼
                     loaded Makie extension
                                   │
                                   ▼
build_context → build_shell → draw!
    → place_legend! → place_colorbars! → build_widgets!
    → assemble → display!
```

Loading `LineCableModels` does not load Makie or select a backend. Recipe
construction cannot create Makie objects. Load CairoMakie, GLMakie, or
WGLMakie before calling `plot` or `preview`.

PlotBuilder is a developer API in v0.2 and may evolve before 1.0. The
user-facing scientific accessors and plotting entry points have their normal
SemVer guarantees.

## Ownership

PlotBuilder owns the fixed core operation sequence and four small declarations:
`LegendDefinition`, `ColorbarDefinition`, `ExportDefinition`, and
`AbstractWidgetDefinition`. Domain modules own plot-definition identities,
accepted inputs, keyword grammar, scientific publication, page cardinality,
and detached payload types. The Makie extension owns figures, axes, drawing,
interaction, responsive layout, widgets, and SVG replay.

`PlotRecipe` is a final product, not an intermediate compiler state. Its only
fields are `definition` and `pages`. A `PlotPage` owns `title`, `size`, `key`,
the definition-owned `payload`, and the legend, colorbar, widget, and export
declarations consumed by the renderer shell. Neither value retains the plotted
source, parsed request, resolved request, or backend object.

The payload is allowed to contain completed geometry, published observations,
styles, and captured display state. It contains drawing data, not shell
declarations, and must not require the Makie extension to reopen an owned
scientific result.

A rendered page returns `UIPlot(render, page, context)`. Its `UIContext` is the
only mutable owner of the figure, shell, canvas, panels, widgets, legend,
colorbars, status, observers, window, and current export settings. The familiar
`plot.figure`, `plot.panels`, and `plot.controls` properties forward to that
context; `UIPlot` does not store duplicate copies.

Public backend selectors accept symbols and immediately route to the same
`Val` method family used for availability, activation, restoration, and screen
construction. The backend specification for each tag holds its extension and
package names, so there is no parallel registry or symbol switch.

Legend and colorbar placement are renderer stages dispatched on the plot-
definition type. Their standard methods use the shell side dock. A definition
renderer may overload `place_legend!` or `place_colorbars!` and either delegate
to the two-argument standard method after adjusting the shell layout or replace
the stage when it needs another Makie arrangement. Placement is not inferred
from the definition-owned payload.

## Maintained families

| Owner | Definition | Accepted source | Entry point |
|:--|:--|:--|:--|
| `Engine` | `LineParameterPlotDefinition` | `SeriesImpedance`, `ShuntAdmittance`, or `LineParameters` | `plot(...)` |
| `Engine` | `LineParametersBenchmarkPlotDefinition` | two `LineParameters` results | benchmark plots |
| `DataModel` | `CablePreviewPlotDefinition` | `CableDesign` | `preview(design)` |
| `DataModel` | `CableCollectionPreviewPlotDefinition` | vector of `CableDesign`s | `preview(designs; layout)` |
| `DataModel` | `SystemPreviewPlotDefinition` | `LineCableSystem` | `preview(system)` |

Material scales are a qualified DataModel definition reused by previews. They
are not a separate root plotting entry point.

`MonteCarloResult` does not add a distribution definition or an operation-tag
grammar. `Makie.hist`, `Makie.stairs`, `Makie.ecdfplot`, `Makie.lines`, and
`Makie.qqplot` publish one explicit marginal and draw the corresponding native
primitive through the standard shell. Each non-mutating call returns `UIPlot`.

## Constructing a detached recipe

The example uses a deterministic two-conductor fixture. CairoMakie renders a
selected completed page after recipe construction.

````@example plotbuilder
using LineCableModels
using LineCableModels.PlotBuilder
using LineCableModels.Engine: LineParameterPlotDefinition
using CairoMakie #hide
CairoMakie.activate!() #hide
nothing; #hide

frequency = 10.0 .^ range(log10(50.0), log10(100_000.0); length = 80) #hide
frequency_scale = frequency ./ first(frequency) #hide
resistance_curve = 1.2e-4 .* (1 .+ 0.12 .* sqrt.(frequency_scale)) #hide
inductance_curve = 4.5e-7 .* (1 .- 0.025 .* log10.(frequency_scale)) #hide
conductance_curve = 2.0e-9 .* (1 .+ 0.08 .* log10.(frequency_scale)) #hide
capacitance_curve = 2.8e-10 .* (1 .- 0.01 .* log10.(frequency_scale)) #hide
Rvalues = zeros(2, 2, length(frequency)) #hide
Lvalues = zeros(2, 2, length(frequency)) #hide
Gvalues = zeros(2, 2, length(frequency)) #hide
Cvalues = zeros(2, 2, length(frequency)) #hide
for index in eachindex(frequency) #hide
    Rvalues[:, :, index] .= resistance_curve[index] .* [1.0 0.18; 0.18 1.15] #hide
    Lvalues[:, :, index] .= inductance_curve[index] .* [1.0 0.22; 0.22 1.08] #hide
    Gvalues[:, :, index] .= conductance_curve[index] .* [1.0 0.12; 0.12 1.1] #hide
    Cvalues[:, :, index] .= capacitance_curve[index] .* [1.0 0.16; 0.16 1.05] #hide
end #hide
omega = reshape(2π .* frequency, 1, 1, :) #hide
parameters = LineParameters( #hide
    complex.(Rvalues, Lvalues .* omega), #hide
    complex.(Gvalues, Cvalues .* omega), #hide
    frequency #hide
) #hide

function documentation_figure( #hide
        recipe::PlotRecipe, #hide
        page_index::Integer = 1; #hide
        export_mode::Bool = false, #hide
        export_theme = nothing #hide
) #hide
    extension = Base.get_extension(LineCableModels, :LineCableModelsMakieExt) #hide
    selected = PlotRecipe( #hide
        recipe.definition, #hide
        recipe.pages[[Int(page_index)]] #hide
    ) #hide
    handles = extension.UIComponents.build( #hide
        selected; #hide
        backend = :cairo, #hide
        display = false, #hide
        controls = false, #hide
        export_mode, #hide
        export_theme #hide
    ) #hide
    return only(handles).figure #hide
end; #hide
nothing #hide
````

`parse` merges definition-owned defaults, rejects unsupported keywords, and
separates semantic input from renderer settings. It does not create a recipe.

````@example plotbuilder
parsed = PlotBuilder.parse(
    LineParameterPlotDefinition,
    parameters;
    requests = (
        @observe((Z, abs)[:, :, :]),
        @observe((Z, angle)[:, :, :])
    ),
    export_theme = :publication
)

(;
    requests = parsed.input.requests,
    figure_size = parsed.renderer.fig_size,
    export_theme = parsed.renderer.export_theme
)
````

The keys returned by `input_defaults` declare semantic keywords. The keys from
`renderer_defaults` declare definition-specific renderer keywords. The only
common renderer keywords are `export_theme` and `open_export`.

`make_render` owns the complete core sequence. `entitle` rejects an unsupported
source before parsing can cause later work. `resolve` validates and completes
the request. `fetch` constructs detached pages. `finish` validates the page
collection and returns the final recipe.

````@example plotbuilder
recipe = make_render(
    LineParameterPlotDefinition,
    parameters;
    requests = (
        @observe((Z, abs)[:, :, :]),
        @observe((Z, angle)[:, :, :])
    )
)

(;
    recipe_fields = fieldnames(typeof(recipe)),
    page_fields = fieldnames(typeof(first(recipe.pages))),
    page_count = length(recipe.pages),
    page_titles = getproperty.(recipe.pages, :title),
    payload_types = typeof.(getproperty.(recipe.pages, :payload))
)
````

Recipe construction above created no figure. CairoMakie consumes one detached
page below.

````@example plotbuilder
documentation_figure(recipe) #hide
````

```@raw html
<br>
```

## Scientific observations and units

Line and Monte Carlo definitions publish scientific values with `observables`.
Selector functions such as `Z`, `Y`, `R`, `L`, `G`, and `C` identify the
request. `Units` supplies the physical quantity, display unit, label, symbol,
and scaling. The detached observation contract contains exactly `values`,
`quantity`, and `unit`.

````@example plotbuilder
observed = observables(
    parameters,
    (
        (frequencies, Colon()),
        (@observe Z[:, :, :]),
        (@observe Y[:, :, :])
    )
)
frequency_observation, impedance_observation, admittance_observation = observed
(;
    observation_keys = keys(impedance_observation),
    frequency_count = length(frequency_observation.values),
    impedance_size = size(impedance_observation.values),
    admittance_size = size(admittance_observation.values)
)
````

Domain code may derive display-ready curves or geometry from those published
values before constructing a page. The Makie extension receives those curves
and observations; it does not reconstruct scientific selectors or inspect the
fields of `LineParameters` or `MonteCarloResult`.

## Direct drawing in the standard shell

The extension implements one narrow drawing method per maintained definition:

```julia
function UIComponents.draw!(
        context::UIContext,
        ::Type{MyPlotDefinition},
        page::PlotPage
)
    payload = page.payload
    axis = PlotBuilder.axis!(
        context,
        context.canvas[1, 1],
        payload.x_observation,
        payload.y_observation;
        title = page.title
    )
    lines!(axis, payload.x, payload.y)
    return context
end
```

There is no package-owned axis, series, view, page-layout, or graphics-
primitive grammar between the payload and Makie. Missing drawing support fails
through ordinary dispatch.

`axis!` creates a scientific axis from observation payloads and registers it
with the common formatter, limits, scale controls, legend, interaction, and
export machinery. `register!` attaches those services to an ordinary Makie
axis created directly by the definition.

The shell then runs the fixed sequence shown at the top of this guide. It owns
theme application, the toolbar, responsive side dock, status row, widgets,
display, and SVG replay. `axis!` and `register!` install the shared formatter,
stable limits, and log-scale metadata when an axis enters the shell. Definition
drawing cannot replace that sequence.

## Monte Carlo primitives

A Monte Carlo request identifies the scientific quantity and, for matrix
results, its row, column, and frequency index. The Makie function identifies
the drawing operation:

```julia
marginal = @observe R[1, 1, 2]
Makie.hist(result, marginal; bins=20, normalization=:pdf)
Makie.stairs(result, marginal)
Makie.ecdfplot(result, marginal)
Makie.lines(result, marginal)
Makie.qqplot(result, marginal; qqline=:identity)
```

Each method checks that the result has one outer Gridspace point, calls
`observables` once for the required `samples` or `histograms` product, and
passes the detached observation to `axis!`. The UQ owner computes histogram
CDFs and Q–Q coordinates; the renderer only calls Makie primitives.

## Native canvas

Qualified `PlotBuilder.plotwindow` exposes the same shell when an advanced
caller needs arbitrary Makie layout. The callback receives the native canvas;
nested grids, spanning axes, heatmaps, legends, and colorbars can therefore be
built without adding another package graphics language.

```julia
using GLMakie

windows = PlotBuilder.plotwindow(title="Custom layout") do ui
    left = GridLayout(ui.canvas[1, 1])
    axis = Axis(left[1, 1])
    x = 1:10
    y = rand(10)
    curve = lines!(axis, x, y)
    PlotBuilder.register!(
        ui,
        axis;
        groups=(curve=(curve,)),
        labels=(curve="sample",),
        data=((xdata=x, ydata=y, group=:curve, label="sample"),)
    )
end
```

This is an advanced, qualified interface. `plotwindow`, `axis!`, and
`register!` are not exported from the package root.

## Definition-owned shell extensions

Legend and colorbar placement dispatch on the definition type. A renderer can
retain the standard shell and replace only the placement it owns:

```julia
function UIComponents.place_colorbars!(
        context::UIContext,
        ::Type{MyPlotDefinition},
        page::PlotPage
)
    # Populate context.shell.colorbars with the desired Makie layout.
    return context
end
```

A custom toolbar control is declared by a subtype of
`AbstractWidgetDefinition`. Its renderer method returns the next free toolbar
column and uses the shell helpers, so it shares icon sizing, callback lifetime,
status reporting, and exception handling:

```julia
struct RecomputeWidget <: AbstractWidgetDefinition end

function UIComponents.build_widget!(context, toolbar, column, ::RecomputeWidget)
    button = UIComponents.toolbar_button!(
        context, toolbar, column; key=:recompute, icon=""
    )
    UIComponents.bind_widget_callback!(
        context, button.clicks; success="Recomputed"
    ) do _
        recompute!()
    end
    return column + 1
end
```

A definition that needs a different shell geometry overloads `build_shell`
directly. Material-scale pages use this path for their colorbar-only layout;
ordinary definitions inherit the standard shell.

## Preview extension dispatch

DataModel owns preview traversal. A new layer adds local `preview_shapes` and
`preview_materials` methods beside its type:

```julia
DataModel.preview_materials(layer::MyLayer) = (layer.material_props,)

function DataModel.preview_shapes(layer::MyLayer, context)
    return DataModel.PreviewPolygon[
        make_polygon(layer, context),
    ]
end
```

The preview orchestrator assigns stable identities, traverses the cable, and
invokes these methods. It contains no switch over the known layer types.

## Responsive state and SVG export

Responsive legends retain the complete semantic series registry. When a
bounded interactive legend cannot fit, the shell shows the largest safe
prefix followed by an inert `(...)` entry. Colorbars are never shortened.
Headless and export rendering collapse interactive-only rows.

SVG export asks the definition-specific replay method for a new detached page
containing the current scale, limits, and visibility state. The shell redraws
that page through explicitly loaded CairoMakie. No backend is imported
dynamically.

- `:default` preserves interactive styling on a white background.
- `:publication` applies the publication font and font sizes.

```julia
plots = plot(parameters; export_theme=:publication)
export_svg(first(plots))
export_svg(first(plots); path="series_impedance.svg", open_file=false)
```

This call renders the publication theme without writing a file:

````@example plotbuilder
publication = make_render(
    LineParameterPlotDefinition,
    parameters;
    requests = (
        @observe(R[:, :, :]),
        @observe(L[:, :, :]),
        @observe(G[:, :, :]),
        @observe(C[:, :, :])
    ),
    export_theme = :publication
)
documentation_figure( #hide
    publication, #hide
    1; #hide
    export_mode = true, #hide
    export_theme = :publication #hide
) #hide
````

```@raw html
<br>
```

Without `path`, export uses `pwd()`. When that directory is inside the package
source, it falls back to `joinpath(tempdir(), "linecablemodels-exports")`.
Filenames are sanitised and timestamped, and existing files are not
overwritten.

## Testing a definition

Test a maintained definition at three boundaries:

1. Construct the recipe without Makie and assert page cardinality, order,
   payload data, observation keys, units, and labels, then assert shell
   declarations on the page.
2. Render with CairoMakie and exercise drawing, formatting, limits, visibility,
   scales, responsive docks, and current-state SVG replay.
3. Compare representative Cairo output with the immutable golden fixtures and
   exercise resize behavior with a real GL display.

Architecture tests also lock both fixed stage orders, inference for maintained
recipe owners, backend-free core loading, extension locality, root API absence,
and the absence of the deleted graphics-compiler vocabulary.
