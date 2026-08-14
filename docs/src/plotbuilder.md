# PlotBuilder

`PlotBuilder` is the declarative recipe layer between LineCableModels domain
objects and Makie. Its purpose is to keep physical selection, unit conversion,
and plot composition independent of a plotting backend.

The central boundary is:

```text
domain object and accessors
        │
        ▼
make_render(::Type{<:AbstractPlotSpec}, object; options...)
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
domain can define an `AbstractPlotSpec` and a `make_render` method without
changing the renderer, provided its output uses the supported declarative
vocabulary.

Recipe decisions use Julia dispatch and explicit, inspectable specification
values. Selection rules remain in each recipe, while rendering behavior remains
in the Makie extension.

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
state.

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

There is no hidden list of required trait methods. A recipe adds one method:

```julia
PlotBuilder.make_render(::Type{MyPlotSpec}, object::MyType; options...)
```

The method states every plotting decision through the specification hierarchy:

| Recipe decision | Where it is defined |
|:--|:--|
| accepted domain type | The `object::MyType` dispatch signature. |
| semantic options and defaults | Typed `make_render` keywords and validation. |
| primitive | `SeriesSpec.kind`. |
| quantity and unit | `AxisSpec.quantity` and `AxisSpec.units`. |
| log-scale availability | `AxisSpec.scale` plus `controls.xlog` or `controls.ylog`. |
| title and legend text | `ViewSpec.title`, `PageSpec.title`, and `SeriesSpec.label`. |
| overlay, panels, or pages | How series are partitioned into views and views into pages. |
| data selection | Calls to domain accessors inside `make_render`. |
| Cartesian or polar values | UnitHandler selection and native `real`, `imag`, `abs`, and `angle`. |
| Makie styling | `SeriesSpec.attributes`, `ViewSpec.attributes`, and `PageSpec.kwargs`. |

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

The following is the minimum shape of a new recipe using existing primitives:

```@example plotbuilder
using LineCableModels

const PB = LineCableModels.PlotBuilder
const UH = LineCableModels.UnitHandler

struct ProfileResult
    frequency::Vector{Float64}
    response::Matrix{Float64}
end

struct ProfilePlotSpec <: PB.AbstractPlotSpec end

function PB.make_render(
    ::Type{ProfilePlotSpec},
    result::ProfileResult;
    color=:steelblue,
    grouping::Symbol=:overlay,
    size=(800, 400),
    export_theme::Symbol=:default,
    open_export::Bool=true,
)
    grouping in (:overlay, :panels, :pages) ||
        throw(ArgumentError("grouping must be :overlay, :panels, or :pages"))
    export_theme in (:default, :publication) ||
        throw(ArgumentError("export_theme must be :default or :publication"))
    xunit = UH.units(:base, :hertz)
    xaxis = PB.AxisSpec(
        :x,
        UH.QuantityTag{:freq}(),
        xunit,
        "Frequency [$(UH.get_label(xunit))]",
        :linear,
    )
    yaxis = PB.AxisSpec(
        :y,
        UH.QuantityTag{:dimensionless}(),
        UH.Units(),
        "Response",
        :linear,
    )
    series = [
        PB.SeriesSpec(
            :line,
            result.frequency,
            result.response[:, index],
            nothing,
            "response $index";
            attributes=(;
                color,
                linewidth=2,
                linestyle=index == 2 ? :dash : :solid,
                group=Symbol("response_$index"),
            ),
        ) for index in axes(result.response, 2)
    ]
    make_view(items, title, key) = PB.ViewSpec(
        xaxis,
        yaxis,
        nothing,
        title,
        items,
        key,
    )
    views = grouping === :overlay ?
            [make_view(series, "Frequency responses", (; group=:all))] :
            [make_view([item], "Response $index", (; response=index))
             for (index, item) in enumerate(series)]
    make_page(page_views, title, export_name) = PB.PageSpec(
        title,
        size,
        grouping === :panels ? :grid : :single,
        page_views,
        (;
            export_theme,
            open_export,
            controls=PB.control_definitions(xlog=true, ylog=true),
            export_name,
            configuration=(; color, grouping),
        ),
    )
    pages = grouping === :pages ?
            [make_page([view], view.title, "response_$index")
             for (index, view) in enumerate(views)] :
            [make_page(views, "Frequency responses", "frequency_responses")]
    return PB.RenderSpec(ProfilePlotSpec, pages)
end

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
nothing
```

The hierarchy defines grouping without a renderer-specific grouping trait:

- multiple `SeriesSpec` values in one `ViewSpec` overlay in one panel;
- multiple `ViewSpec` values in one `PageSpec` produce separate panels in one
  figure;
- multiple `PageSpec` values in one `RenderSpec` produce separate figures or
  native windows;
- equal `attributes.group` values combine several primitives into one legend
  entry and one visibility action.

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
