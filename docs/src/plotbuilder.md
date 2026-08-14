# PlotBuilder guide

`PlotBuilder` is the declarative recipe layer between LineCableModels domain
objects and Makie. Its purpose is to keep physical selection, unit conversion,
and plot composition independent of a plotting backend.

The central boundary is:

```text
domain object and accessors
        │
        ▼
make_render(::Type{<:AbstractPlotSpec}, object; options...)  # defined once
        │
        ├─ parse_kwargs → resolve_input
        ├─ recipe_mode → grouping_mode
        └─ make_pages → make_views → make_series
        │
        ▼
RenderSpec → PageSpec → ViewSpec → SeriesSpec
        │
        ▼
Makie extension and UI renderer
        │
        ▼
Vector{UIPlot}
```

The first half is loaded by `using LineCableModels` and does not depend on
Makie. The second half becomes available only after the user explicitly loads
Makie and one of CairoMakie, GLMakie, or WGLMakie. Recipes must never import a
backend or construct Makie objects.

## Scope and stability

The general-purpose intent is retained at the recipe boundary: the renderer
knows about axes, series, pages, controls, and colorbars, but not about
`CableDesign`, `LineParameters`, Monte Carlo results, or their storage. A new
domain defines an `AbstractPlotSpec` and specializes recipe accessors without
replacing `make_render` or changing the renderer.

The generic pipeline owns option parsing and the complete
page/view/series hierarchy. Recipe types provide semantic decisions through
multiple dispatch. Plot modes and grouping strategies are `Val`-dispatched, so
adding a mode means adding methods for that value rather than branching inside
a recipe builder.

`PlotBuilder` is a developer API in v0.2. The stable user entry points remain
`plot`, `preview`, `show_material_scale`, and `export_svg`.

## Supported recipes

| Owning module | Plot specification | Accepted input | Public entry point |
|:--|:--|:--|:--|
| `Engine` | `LineParameterPlotSpec` | `SeriesImpedance` plus frequencies | `plot(series, frequencies)` |
| `Engine` | `LineParameterPlotSpec` | `ShuntAdmittance` plus frequencies | `plot(shunt, frequencies)` |
| `Engine` | `LineParameterPlotSpec` | `LineParameters` | `plot(parameters)` |
| `UQ` | `MCDistributionPlotSpec` | `CableConstantsMC` | `plot(result, quantity)` |
| `UQ` | `MCDistributionPlotSpec` | `LineParametersMC` | `plot(result, quantity; ijk)` |
| `DataModel` | `CablePreviewPlotSpec` | `CableDesign` | `preview(design)` |
| `DataModel` | `SystemPreviewPlotSpec` | `LineCableSystem` | `preview(system)` |
| `DataModel` | `MaterialScalePlotSpec` | `nothing` | `show_material_scale()` |

Line-parameter recipes support RLCG, Cartesian Z/Y, and polar Z/Y pages.
Monte Carlo recipes support histogram, probability-density, ECDF, and Q-Q
views. The preview recipes support cable geometry, cable-system cross-sections,
physical material colorbars, and layer visibility.

## Declarative object model

### `AbstractPlotSpec`

A zero-field subtype identifies a recipe family and provides the dispatch key
for `make_render`:

```julia
struct TemperatureProfilePlotSpec <: PlotBuilder.AbstractPlotSpec end
```

The type identifies the grammar. It should not contain plot data or mutable UI
state, and it should not overload `make_render`.

### `AxisSpec`

An `AxisSpec` contains:

| Field | Meaning |
|:--|:--|
| `dim` | `:x` or `:y` in the current two-dimensional renderer. `:z` is reserved. |
| `quantity` | A semantic `UnitHandler.QuantityTag`. |
| `units` | The display `UnitHandler.Units`; recipe data must already be scaled to it. |
| `label` | Complete human-readable label, normally including the unit. |
| `scale` | `:linear` or `:log10`. |

`quantity` and `units` remain available for inspection and future renderers;
the Makie renderer does not infer or convert units after receiving an
`AxisSpec`.

### `SeriesSpec`

A `SeriesSpec` describes one backend-neutral plotting primitive:

| `kind` | Data fields | Makie operation | Notes |
|:--|:--|:--|:--|
| `:line` | `xdata`, `ydata` | `lines!` | `Measurement` uncertainty produces x or y error bars. |
| `:scatter` | `xdata`, `ydata` | `scatter!` | Uses nominal values for measurements. |
| `:histogram` | `xdata` | `hist!` | `ydata` and `zdata` are unused. |
| `:stairs` | `xdata`, `ydata` | `stairs!` | Used for retained histogram PDFs. |
| `:heatmap` | `xdata`, `ydata`, `zdata` | `heatmap!` | `zdata` contains cell values. |
| `:polygon` | `zdata` | `poly!` | `zdata` contains GeometryBasics geometry, including polygons with holes. |
| `:hline` | `ydata` | `hlines!` | Used for physical reference lines. |

`label` supplies the legend label. All `attributes` except `group` are passed
to the corresponding Makie plotting function. Consequently, attributes such
as `color`, `linewidth`, `linestyle`, `markersize`, `normalization`, `bins`,
`strokecolor`, and `strokewidth` follow Makie's accepted values.

`group` is a PlotBuilder attribute. Series with the same group are represented
by one legend entry and are toggled together. Group order follows first
declaration order, so recipes should emit series in their intended physical
layer order.

### `ViewSpec`

A `ViewSpec` is one two-dimensional axis system. Its series share the same
axes. `key` is an arbitrary `NamedTuple` identifying the semantic facet, for
example `(; component=:R)` or `(; quantity=:R, selection=(1, 1, 2))`.

The renderer reserves two view attributes:

| Attribute | Meaning |
|:--|:--|
| `aspect=:data` | Use Makie's `DataAspect()` and preserve a physical 1:1 scale. |
| `limits=(xlimits, ylimits)` | Apply explicit axis limits. Each limit is a two-element tuple. |

Other view attributes are forwarded to `Makie.Axis`. Without explicit limits,
Makie supplies its normal margins for varying data. PlotBuilder adds stable
limits for effectively constant data and includes `Measurement` uncertainty
when determining those limits.

### `PageSpec`

A `PageSpec` represents one logical figure or window. It owns a title, pixel
size, layout label, views, and page-level keyword data. Current recipes use
`:single`, `:preview`, and `:material_scale` as semantic layout labels. The
renderer currently arranges multiple views in an automatic grid; `layout` is
not yet a general arbitrary-layout language.

Recognized page keywords are:

| Keyword | Meaning |
|:--|:--|
| `controls` | Result of `PlotBuilder.control_definitions`. |
| `display_legend` | Whether to build the side legend. Defaults to `true`. |
| `colorbars` | Tuple of declarative colorbar descriptors. |
| `export_name` | Base name used for timestamped SVG exports. |
| `export_theme` | `:default` for a faithful export or `:publication` for LaTeX-oriented typography. |
| `open_export` | Open a successful export with the registered system application. Defaults to `true`. |
| `x_exponent`, `y_exponent` | Base-ten display exponents for compact linear tick labels. |
| `configuration` | Recipe options retained for inspection and parity tests. |

A colorbar descriptor contains `label`, `colormap`, `limits`, and `ticks`.
`ticks` is `(positions, labels)`. Labels written as scientific `e` notation
are rendered with typographic superscripts.

The control declaration contains `reset`, `export_svg`, `xlog`, `ylog`,
`legend`, `visibility`, and `zoom` Boolean fields. Reset, SVG, log toggles, and
legend construction are explicit widgets. Visibility is currently provided by
clicking legend entries, and zoom uses the backend's native axis interaction;
there are no separate visibility or zoom widgets.

### `RenderSpec` and `UIPlot`

`RenderSpec` is the complete backend-neutral product of a recipe. Its public
data are:

```julia
render.spec       # the AbstractPlotSpec subtype
render.figures    # Vector{PageSpec}
```

The Makie extension builds one `UIPlot` per page. A `UIPlot` retains both
sides of the boundary:

```julia
handle.render     # complete RenderSpec
handle.page       # PageSpec represented by this window
handle.figure     # built Makie Figure
handle.panels     # built axes and plotted objects
handle.controls   # controls keyed by :reset, :export_svg, :xlog, ...
handle.context    # backend, status, and native-window context
```

Line-parameter plotting returns `Vector{UIPlot}` because each physical
quantity has its own page. Statistical plots and previews return one
`UIPlot`.

## SVG export

The save control reconstructs the current declarative page, including axis
scales, axis limits, and series visibility, and renders it with explicitly
loaded CairoMakie. It never captures the interactive window as a bitmap.

Recipes select one of two export themes:

- `export_theme=:default` keeps the current sizes, labels, styles, and limits
  on a clean white background;
- `export_theme=:publication` applies `Makie.theme_latexfonts()` together with
  the established title, axis-label, tick-label, legend, and colorbar sizing.

The selection is stored in each `PageSpec` and can be overridden for one call:

```julia
plots = plot(parameters; export_theme=:publication)
export_svg(first(plots))

# Batch or headless use
export_svg(first(plots); path="series_resistance.svg", open_file=false)
```

Without an explicit path, PlotBuilder saves a sanitized, timestamped SVG in
`pwd()`. If `pwd()` is the LineCableModels package directory or any directory
inside it, PlotBuilder instead uses
`joinpath(tempdir(), "linecablemodels-exports")`. An explicit path is always
honored and an existing file is never overwritten.

After a successful interactive export, PlotBuilder asks the operating system
to open the SVG with its registered application: `open` on macOS, the Windows
shell association on Windows, and `xdg-open` or `gio open` on Linux. The status
bar and terminal report the absolute path. Set `open_export=false` in the
recipe or `open_file=false` in `export_svg` for batch jobs.

## Recipe interface and accessors

`make_render` is defined once by PlotBuilder. It checks `dispatch_on`, parses
the declared keywords, resolves the recipe, dispatches its plotting and
grouping modes, and materializes the hierarchy. A recipe specializes only the
decisions it owns.

The accessor families are:

| Decision | Accessors |
|:--|:--|
| accepted object | `dispatch_on` |
| semantic options | `input_kwargs`, `input_defaults`, `resolve_input` |
| renderer options | `renderer_kwargs`, `renderer_defaults` |
| plot specialization | `recipe_mode` |
| hierarchy | `grouping_mode`, `page_facets`, `group_facets` |
| axes | `geom_axes`, `axis_quantity`, `axis_unit`, `axis_label`, `axis_scale` |
| primitives | `plot_kind`, `series_data`, `legend_label`, `series_attributes` |
| views | `default_title`, `view_key`, `view_attributes` |
| pages | `default_figsize`, `figure_layout`, `enable_logscale`, `page_kwargs` |

`recipe_mode` and `grouping_mode` return `Val` values. The built-in grouping
modes are:

- `Val(:overlay)`: all `group_facets` become series in one view and one page;
- `Val(:panels)`: every group facet becomes a view in one page;
- `Val(:pages)`: every group facet becomes a page containing one view;
- `Val(:faceted_pages)`: `page_facets` become pages and `group_facets` become
  the overlaid series inside each page;
- `Val(:empty)`: a page without axes, used by the standalone material scale.

The grouping implementation lives in PlotBuilder. A recipe supplies semantic
facet keys and data accessors; it does not construct pages or branch over
layout choices.

Recipes should use domain accessors instead of depending on container storage.
For maintained result types, the relevant accessors include:

- `frequencies`, `basis`, `domain`, `nconductors`, and `nfrequencies`;
- `Z`, `Y`, `R`, `X`, `L`, `G`, `B`, and `C`;
- `statistics`, `samples`, `distribution`, and `surrogate` for Monte Carlo
  results;
- `has_samples` and `has_distributions` before requesting optional data.

Use `UnitHandler.QuantityTag`, `units`, `default_unit`, `display_unit`,
`scale_factor`, `get_label`, and `get_symbol` for semantic units and labels.
Unit conversion belongs in the recipe, not in the renderer and not in Makie
tick-formatting callbacks.

## Defining a new recipe

The following recipe supports one panel with multiple lines, one panel per
line, or one page per line. It does not construct an `AxisSpec`, `SeriesSpec`,
`ViewSpec`, `PageSpec`, or `RenderSpec`:

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

PB.recipe_mode(::Type{ProfilePlotSpec}, recipe::NamedTuple) = Val(:profile)
PB.grouping_mode(
    ::Type{ProfilePlotSpec},
    ::Val{:profile},
    recipe::NamedTuple,
) = Val(recipe.input.grouping)

PB.group_facets(
    ::Type{ProfilePlotSpec},
    ::Val{:profile},
    recipe::NamedTuple,
    page_key,
) = axes(recipe.object.response, 2)

PB.axis_quantity(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x},
    recipe::NamedTuple, page_key, view_key,
) = UH.QuantityTag{:freq}()
PB.axis_quantity(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y},
    recipe::NamedTuple, page_key, view_key,
) = UH.QuantityTag{:dimensionless}()
PB.axis_unit(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x},
    quantity::UH.QuantityTag, recipe::NamedTuple, page_key, view_key,
) = UH.units(:base, :hertz)
PB.axis_label(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y},
    quantity::UH.QuantityTag, unit::UH.Units, recipe::NamedTuple,
    page_key, view_key,
) = "Response"

PB.series_data(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x},
    recipe::NamedTuple, page_key, view_key, series_key::Int,
) = recipe.object.frequency
PB.series_data(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y},
    recipe::NamedTuple, page_key, view_key, series_key::Int,
) = recipe.object.response[:, series_key]
PB.legend_label(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple,
    page_key, view_key, series_key::Int,
) = "response $series_key"

profile_linestyle(::Val) = :solid
profile_linestyle(::Val{2}) = :dash
PB.series_attributes(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple,
    page_key, view_key, series_key::Int,
) = (;
    color=recipe.input.color,
    linewidth=2,
    linestyle=profile_linestyle(Val(series_key)),
    group=Symbol("response_$series_key"),
)

profile_title(::Val{:overlay}, page_key, view_key) = "Frequency responses"
profile_title(::Val{:panels}, page_key, view_key::Int) = "Response $view_key"
profile_title(::Val{:panels}, page_key, ::Nothing) = "Frequency responses"
profile_title(::Val{:pages}, page_key::Int, view_key) = "Response $page_key"
PB.default_title(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple,
    page_key, view_key,
) = profile_title(Val(recipe.input.grouping), page_key, view_key)

profile_layout(::Val{:overlay}) = :single
profile_layout(::Val{:panels}) = :grid
profile_layout(::Val{:pages}) = :single
PB.figure_layout(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple, page_key,
) = profile_layout(Val(recipe.input.grouping))
PB.default_figsize(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple, page_key,
) = recipe.renderer.size
PB.enable_logscale(
    ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::NamedTuple, page_key,
) = (:x, :y)

result = ProfileResult(
    [50.0, 100.0, 500.0],
    [1.0 1.2; 1.5 1.6; 2.0 2.1],
)
overlay = PB.make_render(ProfilePlotSpec, result; color=:navy)
panels = PB.make_render(ProfilePlotSpec, result; grouping=:panels)
pages = PB.make_render(
    ProfilePlotSpec,
    result;
    grouping=:pages,
    export_theme=:publication,
)
@assert length(only(only(overlay.figures).views).series) == 2
@assert length(only(panels.figures).views) == 2
@assert length(pages.figures) == 2
@assert basename(String(which(PB.make_render, (Type{ProfilePlotSpec}, ProfileResult)).file)) ==
        "grammar.jl"
nothing
```

The generic hierarchy applies the grouping mode:

- `:overlay` puts every response in one view;
- `:panels` puts every response in its own view in one page;
- `:pages` puts every response in its own page;
- equal `attributes.group` values combine several primitives into one legend
  entry and one visibility action.

The recipe owns no hierarchy construction. To add a new plot mode, return a
new `Val` from `recipe_mode` and add methods specialized on that value for the
facets, axes, data, labels, or styles that differ. Existing modes do not need
editing.

Inside LineCableModels, connect the recipe to the Makie extension by adding a
narrow `plot` or `preview` method in `ext/LineCableModelsMakieExt.jl`. That
method should only normalize Makie-facing arguments, call `make_render`, and
pass the resulting specification to `UIComponents.build`.

`UIComponents.build` is currently an internal extension entry point. New
LineCableModels recipes connect to it through a narrow method in
`ext/LineCableModelsMakieExt.jl`.

## Customizing an existing recipe

Prefer semantic recipe options over post-build Makie mutation. Existing
customization points include:

- line parameters: `mode`, `coord`, `freq_unit`, `length_unit`,
  `quantity_units`, `con`, `fig_size`, `xscale`, `yscale`, `export_theme`, and
  `open_export`;
- Monte Carlo plots: `quantity`, `ijk`, `mode`, `data`, `length_unit`,
  `quantity_units`, `nbins`, `normalization`, `fig_size`, `export_theme`, and
  `open_export`;
- cable previews: offsets, size, ID, legend, and colorbar visibility;
- system previews: earth model, zoom factor, size, ID, legend, and colorbar
  visibility.

When several recipes need the same behavior, generalize at the narrowest
declarative level:

1. Move physical selection or scaling into a shared domain/UnitHandler helper.
2. Move repeated series construction into a helper returning `SeriesSpec`
   values.
3. Add a renderer feature only when it is independent of a particular domain.
4. Keep UI state out of the recipe; SVG export reconstructs a page from the
   current UI state.

Adding a new primitive currently requires a renderer change in
`UIComponents._draw!`, because primitive dispatch is a closed `kind` switch.
Such a change must define the required data fields, forward only valid Makie
attributes, and add backend-neutral specification tests plus Cairo and GL
coverage. Existing primitives do not require renderer changes.

## Tests for a recipe

Every recipe should be tested at three levels:

1. **Specification semantics:** call `make_render` without Makie and compare
   data, quantities, units, labels, grouping, controls, and page counts.
2. **Renderer semantics:** build without display and exercise controls,
   visibility, axis scaling, and SVG export programmatically.
3. **Visual parity:** compare representative Cairo output to tolerant golden
   images and include the recipe in the manual GL gallery when it is
   interactive.

This separation is what keeps PlotBuilder useful beyond the current cable
plots: new scientific domains can reuse the declarative and rendering layers
without teaching either layer about their object model.
