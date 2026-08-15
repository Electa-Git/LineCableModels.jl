# PlotBuilder guide

`PlotBuilder` is the backend-neutral recipe layer used by LineCableModels. It
turns scientific objects into declarative specifications; the Makie extension
renders those specifications and owns all interactive state.

```text
domain object
    │
    ▼
one generic make_render
    │  recipe accessors and Val dispatch
    ▼
RenderSpec → PageSpec → ViewSpec → AxisSpec / SeriesSpec
    │  named LayoutSpec and typed page components
    ▼
Makie extension → UIPlot
```

Loading `LineCableModels` does not load Makie or a backend. Recipes must not
construct Makie objects. Users explicitly load CairoMakie, GLMakie, or
WGLMakie before calling `plot`, `preview`, or `show_material_scale`.

PlotBuilder is a documented developer API in v0.2 and may evolve before 1.0.
The user-facing result accessors and plotting entry points have their normal
SemVer guarantees.

## Supported recipe families

| Module | Specification | Input | Entry point |
|:--|:--|:--|:--|
| `Engine` | `LineParameterPlotSpec` | `SeriesImpedance` and frequencies | `plot(series, frequencies)` |
| `Engine` | `LineParameterPlotSpec` | `ShuntAdmittance` and frequencies | `plot(shunt, frequencies)` |
| `Engine` | `LineParameterPlotSpec` | `LineParameters` | `plot(parameters)` |
| `UQ` | `MCDistributionPlotSpec` | `CableConstantsMC` | `plot(result, quantity)` |
| `UQ` | `MCDistributionPlotSpec` | `LineParametersMC` | `plot(result, quantity; ijk)` |
| `DataModel` | `CablePreviewPlotSpec` | `CableDesign` | `preview(design)` |
| `DataModel` | `SystemPreviewPlotSpec` | `LineCableSystem` | `preview(system)` |
| `DataModel` | `MaterialScalePlotSpec` | `nothing` | `show_material_scale()` |

Line-parameter recipes provide RLCG, Cartesian Z/Y, and polar Z/Y pages. UQ
recipes provide histogram, PDF, ECDF, and Q-Q views. DataModel recipes provide
cable and system previews plus the material-property scale.

## Recipe state and options

Every recipe is represented internally by:

```julia
PlotRecipe{O,I,R}(object, input, renderer)
```

`object` is the domain value. `input` is a typed `NamedTuple` of semantic
options. `renderer` is a typed `NamedTuple` of renderer options.

A recipe declares semantic options with `input_kwargs` and `input_defaults`,
and renderer options with `renderer_kwargs` and `renderer_defaults`. Defaults
must contain exactly the declared keys. A name cannot occur in both groups or
collide with the common renderer options `layout`, `export_theme`, and
`open_export`. Unsupported keywords are errors.

`resolve_input` validates and enriches a `PlotRecipe`. Expensive conversions
or repeated statistical transformations should be performed there once.

## Specification types

### Axes

`AxisSpec` stores its dimension, `UnitHandler.QuantityTag`, display units,
label, current scale, allowed scales, display exponent, and visual attributes.

- Current and allowed scales are `:linear` or `:log10`.
- The current scale must be present in `allowed_scales`.
- Declaring `:log10` requires positive finite plotted values and positive
  uncertainty lower bounds.
- The renderer derives scale toggles from `allowed_scales`; recipes do not
  declare separate log widgets.
- `exponent` applies compact base-ten scaling to linear tick labels. Logarithmic
  ticks are rendered at decades with typographic superscripts.

The attribute bag is for Makie axis styling. Semantic or renderer-owned keys,
including labels, scales, tick definitions, `title`, and `aspect`, are rejected
there. Axis and view visual attributes are merged when the Makie axis is built;
view attributes take precedence on duplicate visual keys.

### Series

`SeriesSpec` represents one primitive:

| `kind` | Required data | Renderer operation |
|:--|:--|:--|
| `:line` | equal-length `xdata`, `ydata` | `lines!`, plus measurement error bars |
| `:scatter` | equal-length `xdata`, `ydata` | `scatter!` |
| `:histogram` | `xdata` | `hist!` |
| `:stairs` | equal-length `xdata`, `ydata` | `stairs!` |
| `:heatmap` | `xdata`, `ydata`, matching `zdata` | `heatmap!` |
| `:polygon` | geometry in `zdata` | `poly!` |
| `:hline` | `ydata` | `hlines!` |

`label`, `group`, and `visible` are semantic fields. Several primitives with
the same `group` form one legend entry and visibility action. Declaration
order determines legend order. Visual Makie keywords such as `color`,
`linewidth`, `linestyle`, `markersize`, `bins`, and `normalization` belong in
`attributes`; `label`, `group`, and `visible` are rejected there.

The renderer dispatches through `draw!(axis, Val(series.kind), series)`. A new
primitive therefore requires a core validation entry, one renderer method,
and specification plus Cairo tests. Recipes using existing primitives need no
renderer changes.

### Views

`ViewSpec` owns up to three axes, a title, series, a semantic `key`, and:

- `placement::PlacementSpec`, selecting a named layout slot and optionally an
  explicit grid area;
- `aspect`, including `:data` for a physical 1:1 aspect ratio;
- `limits`, as `(xlimits, ylimits)` when the recipe defines bounds;
- an attribute bag forwarded to `Makie.Axis`.

`placement`, `aspect`, and `limits` are real fields and are rejected inside the
attribute bag. `aspect` accepts `nothing`, `:data`, or a positive finite ratio.
Explicit limits are finite, strictly increasing x/y pairs and must be positive
where logarithmic scaling is available. Without explicit limits, the renderer
uses Makie's normal margins and stabilizes effectively constant data, including
measurement uncertainty.

### Pages

`PageSpec` contains:

- `title` and pixel `size`;
- a semantic `key::NamedTuple`;
- a backend-neutral `LayoutSpec`;
- `views`;
- `ControlSpec`, `LegendSpec`, `Vector{ColorbarSpec}`, `StatusSpec`, and
  `ExportSpec`.

There is no general-purpose page keyword bag. Recipe identity belongs in
`key`; display and UI behavior belongs in the corresponding typed component.
The Julia field holding `ExportSpec` is named `export_spec` because `export` is
a Julia keyword.

`ControlSpec` declares reset and SVG buttons. Scale toggles come from axes.
`LegendSpec` controls legend visibility and legend-driven series visibility.
`ColorbarSpec` contains a label, colormap, limits, ticks, and destination slot.
Tick positions must be finite, lie within the color limits, and have one label
per position.
`StatusSpec` controls the status line. `ExportSpec` contains the theme, base
filename, and automatic-open preference.

`RenderSpec` contains the recipe type and its pages. `validate(render)` runs
after generic assembly and again before rendering.

## Named layouts

`LayoutSpec` is a named tree of `GridSpec` and `SlotSpec` values:

- `GridSpec` declares a unique name, parent grid, parent `GridArea`, row and
  column tracks, gaps, and padding.
- `SlotSpec` declares a unique content name, parent grid, `GridArea`, and
  horizontal and vertical alignment.
- `GridArea` uses one-based row and column ranges and supports spans.
- `PlacementSpec(:canvas)` requests automatic view placement.
- `PlacementSpec(:canvas, GridArea(...))` requests explicit placement inside
  the slot.

Track sizes are `FixedTrack(pixels)`, `RelativeTrack(weight)`, and
`ContentTrack()`. Built-in presets are `:single`, `:grid`, `:preview`, and
`:material_scale`. Automatic view placement uses
`ceil(Int, sqrt(number_of_views))` columns.

Callers may pass `layout=:grid` or a complete `LayoutSpec`. Precedence is:

1. caller `layout`;
2. recipe `layout_spec`;
3. `:single`.

Validation rejects duplicate names, missing parents, cycles, areas outside
declared tracks, sibling overlap, invalid tracks, missing content slots,
overlapping view placements, and mixed automatic/explicit placement in one
slot. Enabled toolbar, legend, colorbar, status, and view content must all have
destination slots.

Toolbar and status tracks collapse for headless and SVG rendering. The SVG
therefore preserves the declared content layout without interactive chrome.

## Accessor grammar

`make_render` is defined once in `src/plotbuilder/grammar.jl`. Recipes
specialize these accessors:

| Decision | Accessors |
|:--|:--|
| accepted object | `dispatch_on` |
| semantic options | `input_kwargs`, `input_defaults`, `resolve_input` |
| renderer options | `renderer_kwargs`, `renderer_defaults` |
| modes and facets | `recipe_mode`, `grouping_mode`, `page_facets`, `group_facets` |
| axes | `geom_axes`, `axis_quantity`, `axis_unit`, `axis_label`, `axis_scale`, `axis_scales`, `axis_exponent`, `axis_attributes` |
| primitives | `plot_kind`, `series_data`, `legend_label`, `series_group`, `series_visible`, `series_attributes` |
| views | `default_title`, `view_key`, `view_placement`, `view_aspect`, `view_limits`, `view_attributes` |
| pages | `default_figsize`, `layout_spec`, `page_identity`, `control_spec`, `legend_spec`, `colorbar_specs`, `status_spec`, `export_spec` |

`recipe_mode` and `grouping_mode` return `Val` values. Built-in grouping modes
are:

- `Val(:overlay)`: all group facets become series in one view;
- `Val(:panels)`: each group facet becomes a view in one page;
- `Val(:pages)`: each group facet becomes one page;
- `Val(:faceted_pages)`: page facets become pages and group facets become
  overlaid series;
- `Val(:empty)`: one page without axes.

`make_axes`, `make_series`, `make_views`, and `make_pages` remain advanced
backend-neutral hooks for geometry-heavy or unusual recipes. They may return
specifications but may not construct Makie objects or replace `make_render`.

Recipes should read domain objects through their public accessors. For current
result types these include `frequencies`, `basis`, `domain`, `Z`, `Y`, `R`,
`X`, `L`, `G`, `B`, `C`, `statistics`, `samples`, `distribution`, and
`surrogate`. Unit selection and scaling belong in the recipe through
`UnitHandler`; the renderer receives display-ready data.

## Recipe example: overlay, panels, and pages

The following complete example changes grouping by dispatch and never builds
the page hierarchy itself.

```@example plotbuilder
using LineCableModels

const PB = LineCableModels.PlotBuilder
const UH = LineCableModels.UnitHandler

struct ProfileResult
    frequency::Vector{Float64}
    response::Matrix{Float64}
end

struct ProfilePlotSpec <: PB.AbstractPlotSpec end

PB.dispatch_on(::Type{ProfilePlotSpec}) = ProfileResult
PB.input_kwargs(::Type{ProfilePlotSpec}) = (:grouping, :color)
PB.input_defaults(::Type{ProfilePlotSpec}, ::ProfileResult) =
    (; grouping=:overlay, color=:steelblue)
PB.renderer_kwargs(::Type{ProfilePlotSpec}) = (:size,)
PB.renderer_defaults(::Type{ProfilePlotSpec}, ::ProfileResult) =
    (; size=(800, 400))

PB.recipe_mode(::Type{ProfilePlotSpec}, recipe::PB.PlotRecipe) = Val(:profile)
PB.grouping_mode(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe,
) = Val(recipe.input.grouping)
PB.group_facets(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe, page_key,
) = axes(recipe.object.response, 2)

PB.axis_quantity(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x},
    recipe::PB.PlotRecipe, page_key, view_key,
) = UH.QuantityTag{:freq}()
PB.axis_quantity(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y},
    recipe::PB.PlotRecipe, page_key, view_key,
) = UH.QuantityTag{:dimensionless}()
PB.axis_unit(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x},
    quantity::UH.QuantityTag, recipe::PB.PlotRecipe, page_key, view_key,
) = UH.units(:base, :hertz)
PB.axis_label(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y},
    quantity::UH.QuantityTag, unit::UH.Units, recipe::PB.PlotRecipe,
    page_key, view_key,
) = "Response"
PB.axis_scales(
    ::Type{ProfilePlotSpec}, dim::Val, recipe::PB.PlotRecipe,
    series::Vector{PB.SeriesSpec},
) = (:linear, :log10)

PB.series_data(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x},
    recipe::PB.PlotRecipe, page_key, view_key, series_key::Int,
) = recipe.object.frequency
PB.series_data(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y},
    recipe::PB.PlotRecipe, page_key, view_key, series_key::Int,
) = recipe.object.response[:, series_key]
PB.legend_label(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe,
    page_key, view_key, series_key::Int,
) = "response $series_key"
PB.series_group(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe,
    page_key, view_key, series_key::Int,
) = Symbol("response_$series_key")

profile_linestyle(::Val) = :solid
profile_linestyle(::Val{2}) = :dash
PB.series_attributes(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe,
    page_key, view_key, series_key::Int,
) = (;
    color=recipe.input.color,
    linewidth=2,
    linestyle=profile_linestyle(Val(series_key)),
)

profile_title(::Val{:overlay}, page_key, view_key) = "Frequency responses"
profile_title(::Val{:panels}, page_key, view_key::Int) = "Response $view_key"
profile_title(::Val{:panels}, page_key, ::Nothing) = "Frequency responses"
profile_title(::Val{:pages}, page_key::Int, view_key) = "Response $page_key"
PB.default_title(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe,
    page_key, view_key,
) = profile_title(Val(recipe.input.grouping), page_key, view_key)

profile_layout(::Val{:overlay}) = :single
profile_layout(::Val{:panels}) = :grid
profile_layout(::Val{:pages}) = :single
PB.layout_spec(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe, page_key,
) = profile_layout(Val(recipe.input.grouping))
PB.default_figsize(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe, page_key,
) = recipe.renderer.size
PB.colorbar_specs(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PB.PlotRecipe, page_key,
) = [PB.ColorbarSpec(
    "Response scale",
    :viridis,
    (0.0, 2.5),
    ([0.0, 1.25, 2.5], ["0", "1.25", "2.5"]),
)]

result = ProfileResult(
    [50.0, 100.0, 500.0],
    [1.0 1.2; 1.5 1.6; 2.0 2.1],
)
overlay = PB.make_render(ProfilePlotSpec, result; color=:navy)
panels = PB.make_render(ProfilePlotSpec, result; grouping=:panels)
pages = PB.make_render(ProfilePlotSpec, result; grouping=:pages)

@assert length(only(only(overlay.figures).views).series) == 2
@assert length(only(panels.figures).views) == 2
@assert length(pages.figures) == 2
@assert basename(String(which(
    PB.make_render,
    (Type{ProfilePlotSpec}, ProfileResult),
).file)) == "grammar.jl"
nothing
```

The same recipe accepts a caller layout without another recipe method:

```@example plotbuilder
dashboard = PB.LayoutSpec(
    :dashboard,
    [
        PB.GridSpec(
            :root;
            rows=PB.AbstractTrackSize[
                PB.FixedTrack(36), PB.RelativeTrack(), PB.FixedTrack(20),
            ],
            columns=PB.AbstractTrackSize[PB.ContentTrack(), PB.ContentTrack()],
            rowgap=6,
            columngap=6,
            padding=(20, 20, 28, 28),
        ),
        PB.GridSpec(
            :plots;
            parent=:root,
            area=PB.GridArea(2, 1),
            rows=PB.AbstractTrackSize[PB.RelativeTrack()],
            columns=PB.AbstractTrackSize[PB.RelativeTrack()],
        ),
        PB.GridSpec(
            :side;
            parent=:root,
            area=PB.GridArea(2, 2),
            rows=PB.AbstractTrackSize[PB.ContentTrack(), PB.ContentTrack()],
            columns=PB.AbstractTrackSize[PB.ContentTrack()],
            rowgap=4,
        ),
    ],
    [
        PB.SlotSpec(:toolbar, :root, PB.GridArea(1, 1:2)),
        PB.SlotSpec(:canvas, :plots, PB.GridArea(1, 1)),
        PB.SlotSpec(:legend, :side, PB.GridArea(1, 1)),
        PB.SlotSpec(:colorbars, :side, PB.GridArea(2, 1)),
        PB.SlotSpec(:status, :root, PB.GridArea(3, 1:2)),
    ],
)

custom = PB.make_render(
    ProfilePlotSpec,
    result;
    grouping=:panels,
    layout=dashboard,
)
@assert only(custom.figures).layout === dashboard
nothing
```

For explicit faceting, specialize `view_placement` and return
`PlacementSpec(slot, GridArea(rows, columns))`. All views assigned to one slot
must use either automatic placement or explicit placement; mixing the two is
an error.

## SVG export

The save button reconstructs the current typed page state, including scales,
limits, visibility, layout, and placement, and renders it through explicitly
loaded CairoMakie. It never imports CairoMakie dynamically.

- `:default` preserves the interactive styling on a white background.
- `:publication` applies `Makie.theme_latexfonts()` and the established
  publication sizing.

```julia
plots = plot(parameters; export_theme=:publication)
export_svg(first(plots))
export_svg(first(plots); path="series_resistance.svg", open_file=false)
```

Without `path`, export uses `pwd()`. When `pwd()` is inside the package source,
it falls back to `joinpath(tempdir(), "linecablemodels-exports")`. Filenames
are sanitized and timestamped, and existing files are never overwritten.
Interactive export asks the operating system to open the result. Both the
status line and terminal report the absolute path.

## Testing a recipe

Test each recipe at three levels:

1. Build `RenderSpec` without Makie and assert data, units, labels, grouping,
   semantic keys, scales, layout, placements, and typed components.
2. Build with CairoMakie and test callbacks, visibility, scale changes,
   current-state SVG export, and backend restoration.
3. Compare representative Cairo output with tolerant golden images and add
   interactive recipes to the manual GL gallery.

The architecture tests additionally ensure that maintained recipes resolve to
the one generic `make_render`, core construction loads no Makie backend,
layouts affect renderer placement, and primitives render through `Val`
dispatch.
