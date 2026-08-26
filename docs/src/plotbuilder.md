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
build_context → build_shell → draw! → format_axes!
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
fields are `definition` and `pages`. A `PlotPage` contains only `title`,
`size`, `key`, and a definition-owned `payload`. Neither value retains the
plotted source, parsed request, resolved request, or backend object.

The payload is allowed to contain completed geometry, published observations,
styles, shell declarations, and captured display state. It must not require
the Makie extension to reopen an owned scientific result.

## Maintained families

| Owner | Definition | Accepted source | Entry point |
|:--|:--|:--|:--|
| `Engine` | `LineParameterPlotDefinition` | `SeriesImpedance`, `ShuntAdmittance`, or `LineParameters` | `plot(...)` |
| `Engine` | `LineParametersBenchmarkPlotDefinition` | two `LineParameters` results | benchmark plots |
| `UQ` | `MCDistributionPlotDefinition` | `MonteCarloResult` | `plot(result, quantity)` |
| `DataModel` | `CablePreviewPlotDefinition` | `CableDesign` | `preview(design)` |
| `DataModel` | `SystemPreviewPlotDefinition` | `LineCableSystem` | `preview(system)` |

Material scales are a qualified DataModel definition reused by previews. They
are not a separate root plotting entry point.

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
    quantities = (abs, angle),
    export_theme = :publication
)

(;
    quantities = parsed.input.quantities,
    figure_size = parsed.renderer.fig_size,
    export_theme = parsed.renderer.export_theme
)
````

Semantic keywords are declared through `input_kwargs` and `input_defaults`.
Definition-specific renderer keywords use `renderer_kwargs` and
`renderer_defaults`. The only common renderer keywords are `export_theme` and
`open_export`.

`make_render` owns the complete core sequence. `entitle` rejects an unsupported
source before parsing can cause later work. `resolve` validates and completes
the request. `fetch` constructs detached pages. `finish` validates the page
collection and returns the final recipe.

````@example plotbuilder
recipe = make_render(
    LineParameterPlotDefinition,
    parameters;
    quantities = (abs, angle)
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
        frequency = (frequencies, Colon()),
        series_impedance = Z,
        shunt_admittance = Y
    )
)
(;
    observation_keys = keys(observed.series_impedance),
    frequency_count = length(observed.frequency.values),
    impedance_size = size(observed.series_impedance.values),
    admittance_size = size(observed.shunt_admittance.values)
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
        title = payload.title
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
display, and SVG replay. Definition drawing cannot replace that sequence.

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
    quantities = (R, L, G, C),
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
   payload data, observation keys, units, labels, and shell declarations.
2. Render with CairoMakie and exercise drawing, formatting, limits, visibility,
   scales, responsive docks, and current-state SVG replay.
3. Compare representative Cairo output with the immutable golden fixtures and
   exercise resize behavior with a real GL display.

Architecture tests also lock both fixed stage orders, inference for maintained
recipe owners, backend-free core loading, extension locality, root API absence,
and the absence of the deleted graphics-compiler vocabulary.
