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

The old experimental layer of generic tuple parsing and many independent trait
functions is not part of the current implementation. It was replaced by Julia
dispatch plus explicit, inspectable specification values. This removed a large
amount of unproven machinery; it did not couple the renderer to the currently
supported domains.

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

| Field | Contract |
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

| Keyword | Contract |
|:--|:--|
| `controls` | Result of `PlotBuilder.control_definitions`. |
| `display_legend` | Whether to build the side legend. Defaults to `true`. |
| `colorbars` | Tuple of declarative colorbar descriptors. |
| `export_name` | Base name used for timestamped SVG exports. |
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

## Traits and accessors

There is no hidden list of required trait methods. The current extension
contract has one required method:

```julia
PlotBuilder.make_render(::Type{MyPlotSpec}, object::MyType; options...)
```

Earlier PlotBuilder prototypes expressed the same decisions through many trait
functions. Their current equivalents are explicit:

| Earlier trait concept | Current contract |
|:--|:--|
| accepted domain type | The `object::MyType` dispatch signature. |
| plot kind | `SeriesSpec.kind`. |
| semantic input options and defaults | Typed `make_render` keywords and validation. |
| axis quantity and unit | `AxisSpec.quantity` and `AxisSpec.units`. |
| log-scale capability | `AxisSpec.scale` plus `controls.xlog` or `controls.ylog`. |
| figure size | `PageSpec.size`. |
| title and legend labels | `ViewSpec.title`, `PageSpec.title`, and `SeriesSpec.label`. |
| grouping mode | How the method partitions series into views and views into pages. |
| data container and field selection | Calls to domain accessors inside `make_render`. |
| complex representation | Explicit component selection through UnitHandler and native `real`, `imag`, `abs`, and `angle`. |
| renderer keywords | `SeriesSpec.attributes`, `ViewSpec.attributes`, and `PageSpec.kwargs`. |

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
    response::Vector{Float64}
end

struct ProfilePlotSpec <: PB.AbstractPlotSpec end

function PB.make_render(
    ::Type{ProfilePlotSpec},
    result::ProfileResult;
    color=:steelblue,
    size=(800, 400),
)
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
    series = PB.SeriesSpec(
        :line,
        result.frequency,
        result.response,
        nothing,
        "response";
        attributes=(; color, linewidth=2),
    )
    view = PB.ViewSpec(
        xaxis,
        yaxis,
        nothing,
        "Frequency response",
        [series],
        (; quantity=:response),
    )
    page = PB.PageSpec(
        "Frequency response",
        size,
        :single,
        [view],
        (;
            controls=PB.control_definitions(xlog=true, ylog=true),
            export_name="frequency_response",
            configuration=(; color),
        ),
    )
    return PB.RenderSpec(ProfilePlotSpec, [page])
end

result = ProfileResult([50.0, 100.0, 500.0], [1.0, 1.5, 2.0])
render = PB.make_render(ProfilePlotSpec, result; color=:navy)
@assert render isa PB.RenderSpec{ProfilePlotSpec}
@assert only(only(render.figures).views).series[1].attributes.color === :navy
nothing
```

Inside LineCableModels, connect the recipe to the Makie extension by adding a
narrow `plot` or `preview` method in `ext/LineCableModelsMakieExt.jl`. That
method should only normalize Makie-facing arguments, call `make_render`, and
pass the resulting specification to `UIComponents.build`.

At v0.2, `UIComponents.build` is an internal extension entry point rather than
a public PlotBuilder function. External packages can safely define and inspect
custom `RenderSpec` values, but should not reach through `Base.get_extension`
to render them. A public forwarding build hook should be introduced before
advertising PlotBuilder as a third-party recipe ecosystem.

## Customizing an existing recipe

Prefer semantic recipe options over post-build Makie mutation. Existing
customization points include:

- line parameters: `mode`, `coord`, `freq_unit`, `length_unit`,
  `quantity_units`, `con`, `fig_size`, `xscale`, and `yscale`;
- Monte Carlo plots: `quantity`, `ijk`, `mode`, `data`, `length_unit`,
  `quantity_units`, `nbins`, `normalization`, and `fig_size`;
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
