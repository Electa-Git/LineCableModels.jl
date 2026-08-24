# # PlotBuilder guide
#
# `PlotBuilder` is the backend-neutral recipe layer used by LineCableModels. A
# plot definition selects scientific observations and completes one
# `PlotRecipe`; an explicitly loaded Makie extension renders that recipe and
# owns interactive state.
#
# ```text
# domain object + plot definition
#     │
#     ▼
# entitle → parse → resolve → observe
#     │
#     ▼
# axes → series → views → pages → layout → decorate → finish
#     │
#     ▼
# PlotRecipe → Makie extension → UIPlot
# ```
#
# Loading `LineCableModels` does not load Makie or select a backend. Core recipe
# construction cannot create Makie objects. Users load CairoMakie, GLMakie, or
# WGLMakie before calling `plot` or `preview`.
#
# PlotBuilder is a developer API in v0.2 and may evolve before 1.0. The
# user-facing scientific accessors and plotting entry points have their normal
# SemVer guarantees.
#
# ## Architecture
#
# Every maintained plot family follows four rules:
#
# 1. `make_render` is the sole orchestration method. Definitions specialize
#    stage hooks; they do not replace the sequence.
# 2. `PlotRecipe` is the completed backend-neutral representation. Its pages,
#    axes, series, layouts, and controls are implementation details of that
#    representation rather than parallel public render models.
# 3. Definitions read completed calculations through `observables`. Explicit
#    mathematical accessors such as `Z`, `Y`, `R`, `L`, `G`, and `C` are used
#    only for derived selections requested by the caller.
# 4. Makie figures, observables, widgets, and callbacks exist only in the Makie
#    extension. SVG export reconstructs the current typed recipe state.
#
# Public mode, facet, and key-enumeration hook vocabularies are intentionally absent.
# Variation belongs to the definition that owns a recipe family and is resolved
# within the common stages. This keeps presentation branching from becoming a
# second calculation grammar.
#
# ## Maintained recipe families
#
# | Owner | Definition | Input | Entry point |
# |:--|:--|:--|:--|
# | `Engine` | `LineParameterPlotDefinition` | `SeriesImpedance`, `ShuntAdmittance`, or `LineParameters` | `plot(...)` |
# | `Engine` | `LineParametersBenchmarkPlotDefinition` | two `LineParameters` results | benchmark plots |
# | `UQ` | `MCDistributionPlotDefinition` | `MonteCarloResult` | `plot(result, quantity)` |
# | `DataModel` | `CablePreviewPlotDefinition` | `CableDesign` | `preview(design)` |
# | `DataModel` | `SystemPreviewPlotDefinition` | `LineCableSystem` | `preview(system)` |
#
# Material scales are an internal reusable definition used by previews. They
# are not a separate user-facing plotting entry point.
#
# `plot(parameters)` produces separate impedance and admittance pages with real
# and imaginary parts in adjacent views. An accessor tuple selects a different
# representation, for example `(R, L, G, C)` or `(abs, angle)`. Monte Carlo
# plots provide histogram, density, empirical-CDF, histogram-CDF, and Q-Q views
# when the required retained observations are present.
#
# ## Recipe options and completed state
#
# The following deterministic two-conductor fixture keeps the examples focused
# on recipe behavior. CairoMakie is loaded only to render selected completed
# pages later in the guide.
#
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
        recipe.spec, #hide
        recipe.object, #hide
        recipe.input, #hide
        recipe.renderer, #hide
        [recipe.figures[Int(page_index)]] #hide
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

# `parse_kwargs` applies definition-owned defaults, validates caller keywords,
# and separates scientific input from renderer options. `make_render` runs the
# entitlement check before parsing. Unsupported keywords are errors.
#
parsed = parse_kwargs(
    LineParameterPlotDefinition,
    Z(parameters);
    frequencies = frequencies(parameters),
    quantities = (abs, angle),
    export_theme = :publication
)

(;
    object_type = typeof(parsed.object),
    quantities = parsed.input.quantities,
    export_theme = parsed.renderer.export_theme,
    completed_pages = length(parsed.figures)
)
#
# Semantic options are declared by `input_kwargs` and `input_defaults`.
# Definition-specific renderer options use `renderer_kwargs` and
# `renderer_defaults`. The common renderer options are `layout`,
# `export_theme`, and `open_export`. A name cannot occur in both groups, and
# defaults must contain exactly their declared keys.
#
# `resolve_input` validates and enriches parsed input. `observe` then resolves
# the scientific observations consumed by the remaining stages. Expensive
# statistical transformations belong in one of those stages so they are not
# repeated by views or by the renderer.
#
# `make_render` runs the complete fixed sequence and returns another
# `PlotRecipe`, now with validated pages:
#
recipe = make_render(
    LineParameterPlotDefinition,
    Z(parameters);
    frequencies = frequencies(parameters),
    quantities = (abs, angle)
)

(;
    recipe_type = typeof(recipe),
    page_count = length(recipe.figures),
    page_titles = getproperty.(recipe.figures, :title),
    view_titles = [getproperty.(page.views, :title) for page in recipe.figures]
)
#
# The same completed recipe is materialized below by CairoMakie. Recipe
# construction itself created no figure.

documentation_figure(recipe) #hide
#md # ```@raw html
#md # <br>
#md # ```
#
# ## The fixed operation grammar
#
# Definitions specialize the narrowest stage that owns a decision:
#
# | Stage | Responsibility | Main hooks |
# |:--|:--|:--|
# | entitlement | accept one domain type | `dispatch_on`, `entitle` |
# | parsing | split and validate keywords | `input_kwargs`, `input_defaults`, `renderer_kwargs`, `renderer_defaults`, `parse_kwargs` |
# | resolution | normalize semantic input once | `resolve_input` |
# | observation | obtain immutable scientific values | `observe`, `observables` |
# | axes | quantities, units, labels, scales | `geom_axes`, `axis_quantity`, `axis_unit`, `axis_label`, `axis_scale`, `axis_scales`, `axis_exponent`, `axis_attributes`, `make_axes` |
# | series | primitive data and appearance | `plot_kind`, `series_data`, `legend_label`, `series_group`, `series_visible`, `series_attributes`, `make_series` |
# | views | titles, placement, aspect, limits | `default_title`, `view_key`, `view_placement`, `view_aspect`, `view_limits`, `view_attributes`, `make_views` |
# | pages and layout | size, identity, named layout, controls | `default_figsize`, `page_identity`, `layout_spec`, `control_spec`, `legend_spec`, `colorbar_specs`, `status_spec`, `export_spec`, `make_pages` |
# | completion | final decoration and validation | `decorate`, `finish` |
#
# `make_axes`, `make_series`, `make_views`, and `make_pages` are advanced
# backend-neutral hooks for definitions whose geometry cannot be expressed by
# the narrower accessors. They may return recipe components but cannot create
# Makie objects or replace `make_render`.
#
# A new drawing primitive requires core validation and one renderer method
# dispatched on its primitive symbol. Definitions that reuse an existing
# primitive require no renderer change.
#
# ## Scientific observations and units
#
# Plot definitions consume the immutable named tuples returned by
# `observables`. Presentation code does not read result fields as an alternate
# result protocol. `UnitHandler` maps the resulting scientific keys to quantity
# tags, display units, labels, symbols, and scaling. The renderer receives
# display-ready series and does not interpret calculation containers.
#
# For line parameters, the published quantities are independent of their
# presentation:
#
observed = observables(parameters)
(;
    keys = keys(observed),
    frequency_count = length(observed.frequency),
    impedance_size = size(observed.series_impedance),
    admittance_size = size(observed.shunt_admittance)
)
#
# Mathematical accessors remain available when the caller explicitly asks for
# a derived view. They are not replaced by plot-only quantity keys.
#
# ## Layout and responsive state
#
# Layouts are named grid trees owned by the completed recipe. Callers may
# select a maintained preset with `layout=:single`, `:grid`, `:preview`, or
# `:material_scale`; definitions may supply a complete structured layout.
# Caller selection takes precedence over the definition default.
#
# The common renderer validates layouts before rendering. It rejects missing
# destinations, overlapping areas, invalid tracks, and mixed automatic and
# explicit placement. Toolbars and status rows collapse for headless and SVG
# rendering.
#
# Responsive legends preserve semantic series and visibility state. When a
# bounded legend cannot fit, the interactive renderer shows the largest safe
# prefix followed by an inert `(...)` entry. Material scales are not shortened.
# SVG export reconstructs the complete legend and expands the output height if
# required without resizing the interactive figure.
#
# Passing another maintained preset changes placement without changing the
# definition or scientific observations:
#
grid_recipe = make_render(
    LineParameterPlotDefinition,
    parameters;
    quantities = (R, L, G, C),
    layout = :grid
)

(;
    page_layouts = getproperty.(getproperty.(grid_recipe.figures, :layout), :name),
    pages = length(grid_recipe.figures)
)

documentation_figure(grid_recipe) #hide
#md # ```@raw html
#md # <br>
#md # ```
#
# ## SVG export
#
# The export control reconstructs scales, limits, visibility, layout, and
# placement from the current typed state and renders it through explicitly
# loaded CairoMakie. It never imports CairoMakie dynamically.
#
# - `:default` preserves interactive styling on a white background.
# - `:publication` applies the established publication font and sizing theme.
#
# ```julia
# plots = plot(parameters; export_theme=:publication)
# export_svg(first(plots))
# export_svg(first(plots); path="series_resistance.svg", open_file=false)
# ```
#
# The next figure exercises that publication path without writing a file:
#
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
#md # ```@raw html
#md # <br>
#md # ```
#
# Without `path`, export uses `pwd()`. When that directory is inside the
# package source, it falls back to
# `joinpath(tempdir(), "linecablemodels-exports")`. Filenames are sanitized and
# timestamped, and existing files are never overwritten.
#
# ## Testing a definition
#
# Test a maintained definition at three boundaries:
#
# 1. Complete a `PlotRecipe` without Makie and assert scientific data, units,
#    labels, scales, layout, and semantic identities.
# 2. Render with CairoMakie and test callbacks, visibility, scale changes,
#    current-state SVG export, and backend restoration.
# 3. Compare representative Cairo output with tolerant golden images and add
#    interactive cases to the manual GL gallery.
#
# Architecture tests also keep `make_render` singular, ensure core construction
# loads no Makie backend, verify that removed public mode/facet/key-enumeration
# hooks remain absent, and exercise renderer primitive dispatch.
