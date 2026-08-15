# # PlotBuilder guide
#
# `PlotBuilder` is the backend-neutral recipe layer used by LineCableModels. It
# turns scientific objects into declarative specifications; the Makie extension
# renders those specifications and owns all interactive state.
#
# ```text
# domain object
#     │
#     ▼
# one generic make_render
#     │  recipe accessors and Val dispatch
#     ▼
# RenderSpec → PageSpec → ViewSpec → AxisSpec / SeriesSpec
#     │  named LayoutSpec and typed page components
#     ▼
# Makie extension → UIPlot
# ```
#
# Loading `LineCableModels` does not load Makie or a backend. Recipes must not
# construct Makie objects. Users explicitly load CairoMakie, GLMakie, or
# WGLMakie before calling `plot` or `preview`.
#
# PlotBuilder is a documented developer API in v0.2 and may evolve before 1.0.
# The user-facing result accessors and plotting entry points have their normal
# SemVer guarantees.
#
# ## Architecture
#
# PlotBuilder follows five non-negotiable rules:
#
# 1. `make_render` is one generic pipeline. Domain recipes specialize accessors;
#    they do not replace the pipeline.
# 2. Recipe variation uses Julia dispatch, including `Val` dispatch for modes,
#    grouping, and primitive rendering.
# 3. `PlotRecipe`, `RenderSpec`, layouts, axes, series, views, and page components
#    are backend-neutral. Makie objects exist only in the Makie extension.
# 4. Scientific selection, units, labels, grouping, and page identity remain
#    declarative. Interactive state is reconstructed into the same typed model
#    for export.
# 5. Layouts are named grid trees. Recipes provide defaults and callers may
#    supply a preset or `LayoutSpec` without changing a recipe implementation.
#
# These rules apply to every maintained plot family. The user-facing `plot`,
# `preview`, and `export_svg` methods remain narrow adapters over the same recipe
# and renderer path. Material-scale rendering is an internal, reusable preview
# component rather than a separate public entry point.
#
# The developer API may evolve before LineCableModels 1.0, but changes must keep
# the separation between domain recipes, backend-neutral specifications, and
# backend rendering. Architectural tests enforce the single pipeline and the
# absence of Makie from core specification construction.
#
# ## Supported recipe families
#
# | Module | Specification | Input | Entry point |
# |:--|:--|:--|:--|
# | `Engine` | `LineParameterPlotSpec` | `SeriesImpedance` and frequencies | `plot(series, frequencies)` |
# | `Engine` | `LineParameterPlotSpec` | `ShuntAdmittance` and frequencies | `plot(shunt, frequencies)` |
# | `Engine` | `LineParameterPlotSpec` | `LineParameters` | `plot(parameters)` |
# | `UQ` | `MCDistributionPlotSpec` | `CableConstantsMC` | `plot(result, quantity)` |
# | `UQ` | `MCDistributionPlotSpec` | `LineParametersMC` | `plot(result, quantity; ijk)` |
# | `DataModel` | `CablePreviewPlotSpec` | `CableDesign` | `preview(design)` |
# | `DataModel` | `SystemPreviewPlotSpec` | `LineCableSystem` | `preview(system)` |
#
# Line-parameter recipes provide RLCG, Cartesian Z/Y, and polar Z/Y pages. UQ
# recipes provide histogram, PDF, ECDF, and Q-Q views. DataModel recipes provide
# cable and system previews. `MaterialScalePlotSpec` is an internal reusable
# component for composing and testing the material-property scales used by those
# previews; it is not a public plotting entry point.
#
# ## Recipe state and options
#
# Every recipe is represented internally by `PlotRecipe{O,I,R}`. `object` is the
# domain value, `input` is a typed `NamedTuple` of semantic options, and
# `renderer` is a typed `NamedTuple` of rendering options. The following example
# uses the maintained line-parameter recipe. Specification construction itself
# remains backend-neutral; this documentation loads CairoMakie only so later
# sections can show what the specifications render.
#
using LineCableModels
using LineCableModels.PlotBuilder
using LineCableModels.UnitHandler
using LineCableModels.Engine: LineParameterPlotSpec
import LineCableModels.PlotBuilder: axis_label #hide
import LineCableModels.PlotBuilder: axis_quantity #hide
import LineCableModels.PlotBuilder: axis_scales #hide
import LineCableModels.PlotBuilder: axis_unit #hide
import LineCableModels.PlotBuilder: colorbar_specs #hide
import LineCableModels.PlotBuilder: default_figsize #hide
import LineCableModels.PlotBuilder: default_title #hide
import LineCableModels.PlotBuilder: dispatch_on #hide
import LineCableModels.PlotBuilder: group_facets #hide
import LineCableModels.PlotBuilder: grouping_mode #hide
import LineCableModels.PlotBuilder: input_defaults #hide
import LineCableModels.PlotBuilder: input_kwargs #hide
import LineCableModels.PlotBuilder: layout_spec #hide
import LineCableModels.PlotBuilder: legend_label #hide
import LineCableModels.PlotBuilder: recipe_mode #hide
import LineCableModels.PlotBuilder: renderer_defaults #hide
import LineCableModels.PlotBuilder: renderer_kwargs #hide
import LineCableModels.PlotBuilder: resolve_input #hide
import LineCableModels.PlotBuilder: series_attributes #hide
import LineCableModels.PlotBuilder: series_data #hide
import LineCableModels.PlotBuilder: series_group #hide
using CairoMakie #hide
CairoMakie.activate!() #hide
nothing; #hide

# The guide uses one deterministic two-conductor fixture throughout. Its
# construction is hidden in the generated page so the examples stay focused
# on PlotBuilder decisions rather than synthetic matrix bookkeeping.
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
        render_spec, #hide
        page_index::Integer = 1; #hide
        export_mode::Bool = false, #hide
        export_theme = nothing #hide
) #hide
    extension = Base.get_extension(LineCableModels, :LineCableModelsMakieExt) #hide
    selected = RenderSpec( #hide
        render_spec.spec, #hide
        [render_spec.figures[Int(page_index)]] #hide
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

# parse_kwargs applies the recipe defaults, validates the caller's keywords,
# and separates scientific choices from renderer choices.
recipe = parse_kwargs(
    LineParameterPlotSpec,
    parameters;
    mode = :ZY,
    coord = :polar,
    export_theme = :publication
)

(;
    object_type = typeof(recipe.object),
    mode = recipe.input.mode,
    coordinates = recipe.input.coord,
    export_theme = recipe.renderer.export_theme
)
#
# A recipe declares semantic options with `input_kwargs` and `input_defaults`,
# and renderer options with `renderer_kwargs` and `renderer_defaults`. Defaults
# must contain exactly the declared keys. A name cannot occur in both groups or
# collide with the common renderer options `layout`, `export_theme`, and
# `open_export`. Unsupported keywords are errors.
#
# `resolve_input` validates and enriches a `PlotRecipe`. Expensive conversions
# or repeated statistical transformations should be performed there once.
#
# Calling the generic `make_render` resolves that recipe and assembles the typed
# page tree. Inspecting the result is useful when developing a recipe because it
# tests the complete declarative path without opening a window:
#
render = make_render(
    LineParameterPlotSpec,
    parameters;
    mode = :ZY,
    coord = :polar
)

# Polar Z/Y produces magnitude and angle pages for both matrix families.
[(page.title, only(page.views).yaxis.label) for page in render.figures]
#
# The value above proves that the declarative page tree contains the expected
# semantics. The first page below is that same `RenderSpec` materialized by the
# Cairo renderer. The recipe did not construct this figure.

documentation_figure(render, 1) #hide
#md # ```@raw html
#md # <br>
#md # ```
#
# ## Specification types
#
# ### Axes
#
# `AxisSpec` stores its dimension, `UnitHandler.QuantityTag`, display units,
# label, current scale, allowed scales, display exponent, and visual attributes.
#
# - Current and allowed scales are `:linear` or `:log10`.
# - The current scale must be present in `allowed_scales`.
# - Declaring `:log10` requires positive finite plotted values and positive
#   uncertainty lower bounds.
# - The renderer derives scale toggles from `allowed_scales`; recipes do not
#   declare separate log widgets.
# - `exponent` applies compact base-ten scaling to linear tick labels. Logarithmic
#   ticks are rendered at decades with typographic superscripts.
#
# The attribute bag is for Makie axis styling. Semantic or renderer-owned keys,
# including labels, scales, tick definitions, `title`, and `aspect`, are rejected
# there. Axis and view visual attributes are merged when the Makie axis is built;
# view attributes take precedence on duplicate visual keys.
#
# ### Series
#
# `SeriesSpec` represents one primitive:
#
# | `kind` | Required data | Renderer operation |
# |:--|:--|:--|
# | `:line` | equal-length `xdata`, `ydata` | `lines!`, plus measurement error bars |
# | `:scatter` | equal-length `xdata`, `ydata` | `scatter!` |
# | `:histogram` | `xdata` | `hist!` |
# | `:stairs` | equal-length `xdata`, `ydata` | `stairs!` |
# | `:heatmap` | `xdata`, `ydata`, matching `zdata` | `heatmap!` |
# | `:polygon` | geometry in `zdata` | `poly!` |
# | `:hline` | `ydata` | `hlines!` |
#
# `label`, `group`, and `visible` are semantic fields. Several primitives with
# the same `group` form one legend entry and visibility action. Declaration
# order determines legend order. Visual Makie keywords such as `color`,
# `linewidth`, `linestyle`, `markersize`, `bins`, and `normalization` belong in
# `attributes`; `label`, `group`, and `visible` are rejected there.
#
# The renderer dispatches through `draw!(axis, Val(series.kind), series)`. A new
# primitive therefore requires a core validation entry, one renderer method,
# and specification plus Cairo tests. Recipes using existing primitives need no
# renderer changes.
#
# ### Views
#
# `ViewSpec` owns up to three axes, a title, series, a semantic `key`, and:
#
# - `placement::PlacementSpec`, selecting a named layout slot and optionally an
#   explicit grid area;
# - `aspect`, including `:data` for a physical 1:1 aspect ratio;
# - `limits`, as `(xlimits, ylimits)` when the recipe defines bounds;
# - an attribute bag forwarded to `Makie.Axis`.
#
# `placement`, `aspect`, and `limits` are real fields and are rejected inside the
# attribute bag. `aspect` accepts `nothing`, `:data`, or a positive finite ratio.
# Explicit limits are finite, strictly increasing x/y pairs and must be positive
# where logarithmic scaling is available. Without explicit limits, the renderer
# uses Makie's normal margins and stabilizes effectively constant data, including
# measurement uncertainty.
#
# ### Pages
#
# `PageSpec` contains:
#
# - `title` and pixel `size`;
# - a semantic `key::NamedTuple`;
# - a backend-neutral `LayoutSpec`;
# - `views`;
# - `ControlSpec`, `LegendSpec`, `Vector{ColorbarSpec}`, `StatusSpec`, and
#   `ExportSpec`.
#
# There is no general-purpose page keyword bag. Recipe identity belongs in
# `key`; display and UI behavior belongs in the corresponding typed component.
# The Julia field holding `ExportSpec` is named `export_spec` because `export` is
# a Julia keyword.
#
# `ControlSpec` declares reset and SVG buttons. Scale toggles come from axes.
# `LegendSpec` controls legend visibility and legend-driven series visibility.
# `ColorbarSpec` contains a label, colormap, limits, ticks, and destination slot.
# Tick positions must be finite, lie within the color limits, and have one label
# per position.
# `StatusSpec` controls the status line. `ExportSpec` contains the theme, base
# filename, and automatic-open preference.
#
# `RenderSpec` contains the recipe type and its pages. `validate(render)` runs
# after generic assembly and again before rendering.
#
# ## Named layouts
#
# `LayoutSpec` is a named tree of `GridSpec` and `SlotSpec` values:
#
# - `GridSpec` declares a unique name, parent grid, parent `GridArea`, row and
#   column tracks, gaps, and padding.
# - `SlotSpec` declares a unique content name, parent grid, `GridArea`, and
#   horizontal and vertical alignment.
# - `GridArea` uses one-based row and column ranges and supports spans.
# - `PlacementSpec(:canvas)` requests automatic view placement.
# - `PlacementSpec(:canvas, GridArea(...))` requests explicit placement inside
#   the slot.
#
# Track sizes are `FixedTrack(pixels)`, `RelativeTrack(weight)`, and
# `ContentTrack()`. Built-in presets are `:single`, `:grid`, `:preview`, and
# `:material_scale`. Automatic view placement uses
# `ceil(Int, sqrt(number_of_views))` columns.
#
# Callers may pass `layout=:grid` or a complete `LayoutSpec`. Precedence is:
#
# 1. caller `layout`;
# 2. recipe `layout_spec`;
# 3. `:single`.
#
# Validation rejects duplicate names, missing parents, cycles, areas outside
# declared tracks, sibling overlap, invalid tracks, missing content slots,
# overlapping view placements, and mixed automatic/explicit placement in one
# slot. Enabled toolbar, legend, colorbar, status, and view content must all have
# destination slots.
#
# Toolbar and status tracks collapse for headless and SVG rendering. The SVG
# therefore preserves the declared content layout without interactive chrome.
#
# ## Accessor grammar
#
# `make_render` is defined once in
# [`src/plotbuilder/grammar.jl`](https://github.com/Electa-Git/LineCableModels.jl/blob/main/src/plotbuilder/grammar.jl).
# Recipes
# specialize these accessors:
#
# | Decision | Accessors |
# |:--|:--|
# | accepted object | `dispatch_on` |
# | semantic options | `input_kwargs`, `input_defaults`, `resolve_input` |
# | renderer options | `renderer_kwargs`, `renderer_defaults` |
# | modes and facets | `recipe_mode`, `grouping_mode`, `page_facets`, `group_facets` |
# | axes | `geom_axes`, `axis_quantity`, `axis_unit`, `axis_label`, `axis_scale`, `axis_scales`, `axis_exponent`, `axis_attributes` |
# | primitives | `plot_kind`, `series_data`, `legend_label`, `series_group`, `series_visible`, `series_attributes` |
# | views | `default_title`, `view_key`, `view_placement`, `view_aspect`, `view_limits`, `view_attributes` |
# | pages | `default_figsize`, `layout_spec`, `page_identity`, `control_spec`, `legend_spec`, `colorbar_specs`, `status_spec`, `export_spec` |
#
# `recipe_mode` and `grouping_mode` return `Val` values. Built-in grouping modes
# are:
#
# - `Val(:overlay)`: all group facets become series in one view;
# - `Val(:panels)`: each group facet becomes a view in one page;
# - `Val(:pages)`: each group facet becomes one page;
# - `Val(:faceted_pages)`: page facets become pages and group facets become
#   overlaid series;
# - `Val(:empty)`: one page without axes.
#
# `make_axes`, `make_series`, `make_views`, and `make_pages` remain advanced
# backend-neutral hooks for geometry-heavy or unusual recipes. They may return
# specifications but may not construct Makie objects or replace `make_render`.
#
# Recipes should read domain objects through their public accessors. For current
# result types these include `frequencies`, `basis`, `domain`, `Z`, `Y`, `R`,
# `X`, `L`, `G`, `B`, `C`, `statistics`, `samples`, `distribution`, and
# `surrogate`. `UnitHandler.quantity` binds the physical accessors to semantic
# quantity tags; `UnitHandler` then supplies their units, labels, symbols, and
# scaling. It does not extract values from result containers. The renderer
# receives display-ready data.
#
# ## Build a recipe with accessors
#
# This complete recipe has one frequency vector and any number of response
# columns. Every block is executed in one Literate/Documenter session, so the
# example develops in the same order as a real recipe and fails the docs build
# if an accessor signature becomes stale.
#
# ### 1. Declare the domain type and its options
#
# The domain type contains scientific data only. The empty specification type
# identifies the recipe. Its accessors declare the exact accepted object and
# keywords; PlotBuilder rejects anything else before it assembles a page.
#
# The rows of response correspond to frequency samples. Each column will
# become a separate series, panel, or page depending on the grouping mode.
struct ProfileResult
    frequency::Vector{Float64}
    response::Matrix{Float64}
end
nothing; #hide

# A specification type carries no state. Accessor methods supply its behavior.
struct ProfilePlotSpec <: AbstractPlotSpec end
nothing; #hide

# This recipe accepts ProfileResult and exactly two scientific options.
dispatch_on(::Type{ProfilePlotSpec}) = ProfileResult
input_kwargs(::Type{ProfilePlotSpec}) = (:grouping, :color)
function input_defaults(::Type{ProfilePlotSpec}, ::ProfileResult)
    (; grouping = :overlay, color = :steelblue)
end
nothing; #hide

# Figure size affects rendering, so it belongs to renderer options instead.
renderer_kwargs(::Type{ProfilePlotSpec}) = (:size,)
renderer_defaults(::Type{ProfilePlotSpec}, ::ProfileResult) = (; size = (800, 400))

nothing; #hide
#
# Common renderer keywords—`layout`, `export_theme`, and `open_export`—are
# provided by PlotBuilder and must not be redeclared. A real recipe should use
# `resolve_input` to validate its values and to cache any expensive derived data.
# This example validates its array dimensions and grouping through dispatch:
#
# Each supported grouping has a method. The fallback produces an actionable
# error without putting a grouping conditional in the generic pipeline.
profile_grouping(::Val{:overlay}) = Val(:overlay)
profile_grouping(::Val{:panels}) = Val(:panels)
profile_grouping(::Val{:pages}) = Val(:pages)
function profile_grouping(::Val{G}) where {G}
    throw(ArgumentError("unsupported profile grouping :$G"))
end

function resolve_input(
        ::Type{ProfilePlotSpec},
        recipe::PlotRecipe
)
    length(recipe.object.frequency) == size(recipe.object.response, 1) ||
        throw(DimensionMismatch("each response row needs one frequency"))
    all(isfinite, recipe.object.frequency) ||
        throw(ArgumentError("frequencies must be finite"))
    all(isfinite, recipe.object.response) ||
        throw(ArgumentError("responses must be finite"))
    profile_grouping(Val(recipe.input.grouping))
    return recipe
end

nothing; #hide
#
# ### 2. Select the mode, grouping, and axes
#
# `recipe_mode` gives this family a semantic mode. `grouping_mode` selects one of
# the generic overlay, panel, or page assemblers, while `group_facets` identifies
# the response columns to assemble.
#
recipe_mode(::Type{ProfilePlotSpec}, recipe::PlotRecipe) = Val(:profile)
function grouping_mode(
        ::Type{ProfilePlotSpec},
        ::Val{:profile},
        recipe::PlotRecipe
)
    profile_grouping(Val(recipe.input.grouping))
end
function group_facets(
        ::Type{ProfilePlotSpec},
        ::Val{:profile},
        recipe::PlotRecipe,
        page_key
)
    axes(recipe.object.response, 2)
end
nothing; #hide

# Axis accessors describe physical meaning and display units. No Makie object
# is created here; UnitHandler remains the authority for unit labels.
function axis_quantity(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x},
        recipe::PlotRecipe, page_key, view_key
)
    QuantityTag{:freq}()
end
function axis_quantity(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y},
        recipe::PlotRecipe, page_key, view_key
)
    QuantityTag{:dimensionless}()
end
function axis_unit(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x},
        quantity::QuantityTag, recipe::PlotRecipe, page_key, view_key
)
    units(:base, :hertz)
end
function axis_label(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y},
        quantity::QuantityTag, unit::Units, recipe::PlotRecipe,
        page_key, view_key
)
    "Response"
end
nothing; #hide

# Both dimensions are positive in this example, so the renderer may expose a
# linear/logarithmic toggle for each axis.
function axis_scales(
        ::Type{ProfilePlotSpec}, dim::Val, recipe::PlotRecipe,
        series::Vector{SeriesSpec}
)
    (:linear, :log10)
end

nothing; #hide
#
# ### 3. Describe series data and appearance
#
# The generic assembler calls these methods once for every facet. Data, labels,
# groups, and visibility are semantic accessors. Only visual Makie properties
# belong in `series_attributes`.
#
function series_data(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:x},
        recipe::PlotRecipe, page_key, view_key, series_key::Int
)
    recipe.object.frequency
end
function series_data(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, ::Val{:y},
        recipe::PlotRecipe, page_key, view_key, series_key::Int
)
    recipe.object.response[:, series_key]
end
function legend_label(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PlotRecipe,
        page_key, view_key, series_key::Int
)
    "response $series_key"
end
function series_group(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PlotRecipe,
        page_key, view_key, series_key::Int
)
    Symbol("response_$series_key")
end
nothing; #hide

# Dispatch selects a special line style for the second response. Adding a
# third specialization does not change the recipe assembler.
profile_linestyle(::Val) = :solid
profile_linestyle(::Val{2}) = :dash
function series_attributes(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PlotRecipe,
        page_key, view_key, series_key::Int
)
    (;
        color = recipe.input.color,
        linewidth = 2,
        linestyle = profile_linestyle(Val(series_key))
    )
end

nothing; #hide
#
# ### 4. Describe views and pages
#
# Titles and layouts vary by grouping through small dispatch functions. The
# recipe still does not construct `ViewSpec` or `PageSpec`; the one generic
# `make_render` pipeline does that from these decisions.
#
profile_title(::Val{:overlay}, page_key, view_key) = "Frequency responses"
profile_title(::Val{:panels}, page_key, view_key::Int) = "Response $view_key"
profile_title(::Val{:panels}, page_key, ::Nothing) = "Frequency responses"
profile_title(::Val{:pages}, page_key::Int, view_key) = "Response $page_key"
function default_title(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PlotRecipe,
        page_key, view_key
)
    profile_title(Val(recipe.input.grouping), page_key, view_key)
end

profile_layout(::Val{:overlay}) = :single
profile_layout(::Val{:panels}) = :grid
profile_layout(::Val{:pages}) = :single
function layout_spec(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PlotRecipe, page_key
)
    profile_layout(Val(recipe.input.grouping))
end
function default_figsize(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PlotRecipe, page_key
)
    recipe.renderer.size
end
nothing; #hide

# Typed page components are also supplied by accessors. This colorbar docks in
# the preset's :colorbars slot and can be replaced by another specialization.
function colorbar_specs(
        ::Type{ProfilePlotSpec}, ::Val{:profile}, recipe::PlotRecipe, page_key
)
    [ColorbarSpec(
        "Response scale",
        :viridis,
        (0.0, 2.5),
        ([0.0, 1.25, 2.5], ["0", "1.25", "2.5"])
    )]
end

nothing; #hide
#
# ### 5. Assemble overlay, panel, and page variants
#
# All three variants use the same data and recipe methods. Only the grouping
# value changes. These assertions execute during the documentation build.
#
result = ProfileResult(
    [50.0, 100.0, 500.0],
    [1.0 1.2; 1.5 1.6; 2.0 2.1]
)

overlay = make_render(ProfilePlotSpec, result; color = :navy)
panels = make_render(ProfilePlotSpec, result; grouping = :panels)
pages = make_render(ProfilePlotSpec, result; grouping = :pages)

(;
    overlay_series = length(only(only(overlay.figures).views).series),
    panel_views = length(only(panels.figures).views),
    separate_pages = length(pages.figures),
    pipeline_file = basename(String(which(
        make_render,
        (Type{ProfilePlotSpec}, ProfileResult)
    ).file))
)
#
# #### Overlay
#
# `Val(:overlay)` places every response in one `ViewSpec`. The legend remains
# outside the canvas because its destination is declared separately by the
# page's `LegendSpec`.

documentation_figure(overlay) #hide
#md # ```@raw html
#md # <br>
#md # ```
#
# #### Panels
#
# `Val(:panels)` reuses the same two `SeriesSpec` values but gives each one its
# own view. The `:grid` layout preset automatically arranges those views.

documentation_figure(panels) #hide
#md # ```@raw html
#md # <br>
#md # ```
#
# #### Pages
#
# `Val(:pages)` produces one page per response. Showing the second page makes
# the difference from panels explicit without duplicating both figures here.

documentation_figure(pages, 2) #hide
#md # ```@raw html
#md # <br>
#md # ```
#
# ## Override the layout without changing the recipe
#
# The same recipe accepts a caller-provided named layout. This example creates a
# root grid with toolbar and status rows, a plot grid, and a side grid where the
# legend and colorbar are docked.
#
# The root reserves content-sized outer tracks and a flexible center.
root_grid = GridSpec(
    :root;
    rows = AbstractTrackSize[
        FixedTrack(36), RelativeTrack(), FixedTrack(20)
    ],
    columns = AbstractTrackSize[ContentTrack(), ContentTrack()],
    rowgap = 6,
    columngap = 12,
    padding = (20, 20, 28, 28)
)
#
# The canvas occupies the flexible left side of the center row.
plot_grid = GridSpec(
    :plots;
    parent = :root,
    area = GridArea(2, 1),
    rows = AbstractTrackSize[RelativeTrack()],
    columns = AbstractTrackSize[RelativeTrack()]
)
#
# The right side stacks legend and colorbar content.
side_grid = GridSpec(
    :side;
    parent = :root,
    area = GridArea(2, 2),
    rows = AbstractTrackSize[ContentTrack(), ContentTrack()],
    columns = AbstractTrackSize[ContentTrack()],
    rowgap = 4
)
#
# Components target slots by name; the renderer materializes them.
dashboard_slots = [
    SlotSpec(:toolbar, :root, GridArea(1, 1:2)),
    SlotSpec(:canvas, :plots, GridArea(1, 1)),
    SlotSpec(:legend, :side, GridArea(1, 1); halign = :left, valign = :top),
    SlotSpec(:colorbars, :side, GridArea(2, 1); halign = :left, valign = :top),
    SlotSpec(:status, :root, GridArea(3, 1:2))
]
#
# `LayoutSpec` validates the complete named tree before a renderer sees it.
dashboard = LayoutSpec(
    :dashboard,
    [root_grid, plot_grid, side_grid],
    dashboard_slots
)

custom = make_render(
    ProfilePlotSpec,
    result;
    grouping = :panels,
    layout = dashboard
)

(;
    layout = only(custom.figures).layout.name,
    slots = getproperty.(dashboard.slots, :name),
    views = length(only(custom.figures).views)
)
#
# The recipe accessors are unchanged. Only the caller-provided named grid moves
# the canvas, legend, colorbar, toolbar, and status destinations.

documentation_figure(custom) #hide
#md # ```@raw html
#md # <br>
#md # ```
#
# For explicit faceting, specialize `view_placement` and return
# `PlacementSpec(slot, GridArea(rows, columns))`. All views assigned to one slot
# must use either automatic placement or explicit placement; mixing the two is
# an error.
#
# ## SVG export
#
# The save button reconstructs the current typed page state, including scales,
# limits, visibility, layout, and placement, and renders it through explicitly
# loaded CairoMakie. It never imports CairoMakie dynamically.
#
# - `:default` preserves the interactive styling on a white background.
# - `:publication` applies `Makie.theme_latexfonts()` and the established
#   publication sizing.
#
# ```julia
# plots = plot(parameters; export_theme=:publication)
# export_svg(first(plots))
# export_svg(first(plots); path="series_resistance.svg", open_file=false)
# ```
#
# This headless rendering uses the same export path without writing a file. The
# toolbar and status tracks collapse, and `:publication` applies the configured
# LaTeX-font theme to the declarative page state.

publication = make_render(
    LineParameterPlotSpec,
    parameters;
    mode = :RLCG,
    export_theme = :publication
)
documentation_figure( #hide
    publication, #hide
    1; #hide
    export_mode = true, #hide
    export_theme = :publication #hide
) #hide
#md # ```@raw html
#md # <br>
#md # ```
#
# Without `path`, export uses `pwd()`. When `pwd()` is inside the package source,
# it falls back to `joinpath(tempdir(), "linecablemodels-exports")`. Filenames
# are sanitized and timestamped, and existing files are never overwritten.
# Interactive export asks the operating system to open the result. Both the
# status line and terminal report the absolute path.
#
# ## Testing a recipe
#
# Test each recipe at three levels:
#
# 1. Build `RenderSpec` without Makie and assert data, units, labels, grouping,
#    semantic keys, scales, layout, placements, and typed components.
# 2. Build with CairoMakie and test callbacks, visibility, scale changes,
#    current-state SVG export, and backend restoration.
# 3. Compare representative Cairo output with tolerant golden images and add
#    interactive recipes to the manual GL gallery.
#
# The architecture tests additionally ensure that maintained recipes resolve to
# the one generic `make_render`, core construction loads no Makie backend,
# layouts affect renderer placement, and primitives render through `Val`
# dispatch.
