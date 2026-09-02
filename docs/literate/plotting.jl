# # Makie plotting: implementation guide and complete gallery

# Loading a Makie backend adds a small, opinionated convenience layer to the
# ordinary Makie API. The convenience call should get a scientific plot almost
# finished; the returned `Figure`, `Axis`, `Legend`, `Colorbar`, and controls are
# nevertheless native Makie objects owned by the caller.
#
# This page is both an executable gallery and the implementation guide for the
# plotting surface. Every gallery image below is produced by the call printed
# immediately above it.

using LineCableModels
using CairoMakie
using LinearAlgebra: diag
using Measurements: measurement

# ## The implementation contract

# Plotting follows one path rather than maintaining a second declarative plot
# model:
#
# 1. Scientific owners expose public observations, geometry, material values,
#    and units. They contain no plot preparation, colors, legends, or Makie
#    blocks.
# 2. `LineCableModelsMakieExt` normalizes plotting requests into physical
#    quantity/coordinate facets and constructs ordinary Makie blocks and plot
#    primitives. Matrix coordinates identify subplots; result containers
#    identify overlaid series.
# 3. The common addon shell contributes the toolbar, observable-driven scale
#    controls, scientific labels, limits, docks, and export action.
# 4. The call returns a [`UIPlot`](@ref) containing the live native objects. It
#    does not retain a shadow plot specification and does not replay the plot
#    after the caller changes it.
#
# The relevant implementation files are:
#
# | Responsibility | Implementation |
# |:--|:--|
# | Optional entry points and `UIPlot` | `src/plotbuilder/` |
# | Line request normalization | `ext/LineCableModelsMakieExt/recipes/line_data.jl` and `comparison_data.jl` |
# | Domain-only preview geometry and ranges | `src/datamodel/preview/geometry.jl` and `materials.jl` |
# | Preview presentation data | `ext/LineCableModelsMakieExt/recipes/preview_data.jl` |
# | Independent material palettes | `ext/LineCableModelsMakieExt/material_colors.jl` |
# | Shell, axes, limits, docks, widgets | `ext/LineCableModelsMakieExt/shell.jl` |
# | Line, preview, and publication drawing | `ext/LineCableModelsMakieExt/recipes/*_render.jl` |
# | Monte Carlo verbs | `ext/LineCableModelsMakieExt/montecarlo.jl` |
# | Publication SVG snapshot | `ext/LineCableModelsMakieExt/native_export.jl` |
#
# `UIPlot` is deliberately small. Its `figure`, `axes`, `controls`, `legend`,
# `panel_legends`, and `colorbars` fields point at the exact objects on screen.
# A one-page recipe returns that handle directly. A request that genuinely
# produces several figures returns `Vector{UIPlot}`. No `only(...)` wrapper is
# required for a one-page call.

# ### Naming titles and legends

# Layout text has one name per scope. These keywords are presentation metadata;
# none of them changes a result tensor or geometry tag.
#
# | Keyword | Scope |
# |:--|:--|
# | `figure_title` | One visible title above the whole native figure |
# | `title_attributes` | Native Makie `Label` attributes for `figure_title` |
# | `panel_titles` | Axis-title overrides, positionally or by semantic key |
# | `legend_title` | Heading of the controlled figure legend |
# | `series_labels` | Names of overlaid result containers; these are the legend entries |
# | `legend_position` | `:inside`, a named outer dock, or a positive dock grid position |
# | `legend_anchor` | Inside alignment such as `:lt`, `:cc`, or `:rb` |
#
# Some established recipes also accept `title` as their window/export or
# single-recipe heading. Use the explicit scoped names above when composing a
# dashboard. The old line/publication keyword `labels` remains a compatibility
# alias for `panel_titles`, but new code should not use it.
#
# The same scopes are mutable after construction. `figuretitle!` and
# `paneltitle!` replace titles. `figurelegend!` and `panellegend!` rebuild a
# native legend from the retained semantic plot handles, so relabeling does not
# destroy grouped visibility behavior.

# ### Observation, facet, page, and legend semantics

# `LineParameters` owns dense ``Z`` and ``Y`` tensors. The observation grammar
# is what turns those tensors into displayable physical quantities: ``Z``
# expands to ``R`` then ``X``; ``Y`` expands to ``G`` then ``B``. An exact
# request such as `@observe Z[1,1,:]` keeps its conductor coordinates while
# selecting the complete frequency range. Plotting consumes that resolved
# request; it does not ask the caller to repeat `(R, X)`.
#
# Every `(quantity, row, column)` is one facet and therefore one axis. Its
# default title comes from the unit registry plus the coordinate relation, for
# example `Self series resistance — conductor 1` or
# `Mutual series reactance — conductor 1 → conductor 2`. Coordinates never
# become legend entries. A legend exists to distinguish overlaid result sets;
# a single unlabeled result therefore has no legend by default.
#
# `layout` chooses how those already-selected facets are grouped into pages:
#
# | `layout` | Page grouping |
# |:--|:--|
# | `nothing` | Exact complex coordinate requests are paired; matrix selections become one dashboard per quantity |
# | `(1, 1)` | One facet per figure |
# | `(1, 2)` | Real/imaginary or derived family components paired side by side per coordinate |
# | `(2, 1)` | The same pair stacked per coordinate |
# | `(N, N)` | One matrix dashboard per quantity for an `N×N` result |
#
# Layout is applied after observation selection. It cannot add an excluded
# conductor, merge unlike physical units onto one axis, or turn a coordinate
# into a data-series label.

# ## Gallery data

# The frequency responses are deliberately small but non-constant so that
# logarithmic axes, engineering units, legends, and automatic limits are all
# visible in the generated documentation.

frequency = collect(10.0 .^ range(1, 4; length = 24));
angular_frequency = reshape(2π .* frequency, 1, 1, :);
resistance = cat((
    [1.0 0.22; 0.22 1.8] .* 1.0e-4 .* (1 + 0.12 * log10(f / first(frequency)))
    for f in frequency
)...; dims = 3);
inductance = cat((
    [2.0 0.28; 0.28 2.5] .* 1.0e-7 .* (1 - 0.04 * log10(f / first(frequency)))
    for f in frequency
)...; dims = 3);
conductance = cat((
    [3.0 -0.45; -0.45 4.0] .* 1.0e-9 .* (1 + 0.08 * log10(f / first(frequency)))
    for f in frequency
)...; dims = 3);
capacitance = repeat([4.0 -0.7; -0.7 5.0] .* 1.0e-10, 1, 1, length(frequency));
parameters = LineParameters(
    complex.(resistance, inductance .* angular_frequency),
    complex.(conductance, capacitance .* angular_frequency),
    frequency
);

# The geometry gallery uses the example library shipped with the repository.
# A system is assembled from two placements without introducing a plotting-only
# geometry representation.

cable_library = CablesLibrary();
LineCableModels.load!(
    cable_library;
    file_name = joinpath(pkgdir(LineCableModels), "examples", "cables_library.json")
);
mv_design = get(cable_library, "18kV_1000mm2");
hv_design = get(cable_library, "525kV_1600mm2");
earth = EarthModel(100.0, 10.0, 1.0);
cable_system = build(
    LineCableSystem,
    [mv_design, mv_design],
    [(-0.06, -0.20), (0.06, -0.20)];
    environment = earth,
    system_id = "two-cable-gallery",
    line_length = 1_000.0
);

# Monte Carlo plots consume retained result products. The documentation fixture
# therefore contains samples, their summary, and a normalized histogram model;
# none of the plotting calls below reruns a simulation.

sample_values = collect(range(1.0, 5.0; length = 81));
sample_summary = SampleSummary(sample_values);
sample_model = HistogramDensity(
    collect(range(1.0, 5.0; length = 9)),
    fill(0.25, 8)
);
mc_formulation = MonteCarlo(
    Formulation();
    trials = length(sample_values),
    seed = 41,
    return_samples = true,
    return_histograms = true
);
mc_result = MonteCarloResult(
    mc_formulation,
    [CableConstants(3.0, 3.0, 3.0, 3.0)],
    [(R = [sample_summary], L = [sample_summary], C = [sample_summary],
      G = [sample_summary])],
    [(R = reshape(sample_values, 1, :), L = reshape(sample_values, 1, :),
      C = reshape(sample_values, 1, :), G = reshape(sample_values, 1, :))],
    [(R = [sample_model], L = [sample_model], C = [sample_model],
      G = [sample_model])],
    UInt64(41),
    UInt64[41],
    [length(sample_values)]
);

# ## Line-parameter recipes

# ### Complete default view

# `Makie.plot(parameters)` is the minimal call. It observes everything in the
# canonical order ``Z`` then ``Y``, expands those families to ``R``, ``X``,
# ``G``, and ``B``, and returns four matrix-dashboard pages. For this 2×2
# result, every page contains four axes. The calls below render the full default
# gallery in exactly that order.

default_line_pages = Makie.plot(
    parameters;
    backend = :cairo,          # choose the already-loaded native backend
    display_plot = false,      # Documenter owns display; interactive use may omit this
    controls = false,          # omit toolbar chrome from this static gallery
    xscale = :log10,           # initial scale; interactive controls may change it
    fig_size = (900, 680)      # size of each generated page
)
default_line_pages[1].figure #hide
default_line_pages[2].figure #hide
default_line_pages[3].figure #hide
default_line_pages[4].figure #hide

# To change only appearance, mutate the returned native axes. To change which
# physical values exist, pass an observation selector or exact `@observe`
# request. To change pagination, pass one of the layouts in the table above.

# ### Cartesian series impedance

# This exact observation asks for the self impedance of conductor 1 over the
# full frequency range. The observation layer expands ``Z`` into resistance and
# reactance; the inferred one-page layout places the two incompatible physical
# quantities on separate axes.

## `xscale=:log10` is an initial state; the live toolbar can still change it.
## `figure_title` labels the whole figure; `panel_titles` labels its two axes.
## `series_labels` names result sets, so supplying one explicitly opts into a legend.
## `:inside` overlays the union of axis viewports; `:rb` anchors at its corner.
series_cartesian = Makie.plot(
    parameters,
    @observe Z[1, 1, :];       # self impedance, conductor 1, every frequency
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    fig_size = (1000, 500),
    figure_title = "Conductor 1 self impedance",
    panel_titles = ("Resistance", "Reactance"),
    legend_title = "Result set",
    series_labels = ("reference",),
    legend_position = :inside,
    legend_anchor = :rb,
    legend_attributes = (; backgroundcolor = (:white, 0.92)),
    legend_overflow = :show_all
)
series_cartesian.figure #hide

# To modify this recipe, select other matrix coordinates with an explicit
# observation request, pass native legend attributes, or mutate either returned
# axis. For example, `series_cartesian.axes[1].title[] = "Measured resistance"`
# changes the live title without asking LineCableModels to rebuild anything.

# ### Cartesian shunt admittance

# Conductance and susceptance are the real and imaginary parts of ``Y``. They
# use the same dashboard implementation, unit publication, limits, legend
# groups, and controls as impedance; only the scientific requests differ.

## A bottom dock is horizontal by default and remains a native Makie Legend.
shunt_cartesian = Makie.plot(
    parameters,
    @observe Y[1, 1, :];
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    fig_size = (900, 480),
    panel_titles = Dict(:G => "Self conductance", :B => "Self susceptance")
)
shunt_cartesian.figure #hide

# Any accepted Makie `Legend` keyword belongs in `legend_attributes`. For a
# labeled source or comparison, moving the
# block uses `legend_position = :inside`, `:left`, `:right`, `:top`, `:bottom`,
# or a positive `(row, column)` dock coordinate. `(2, 2)` remains the plot
# canvas. `legend_anchor` is used only for an inside legend.

# ### Figure-wide and panel-scoped legends

# A controlled comparison registers one group per result set, then exposes two
# views of those groups. The figure legend spans every panel. A panel legend
# filters the same source handles to one logical plot position. Hiding a global
# result entry changes that source's visibility in both the R and X facets.

legend_scope_demo = Makie.plot(
    parameters,
    parameters,
    @observe Z[1, 1, :];
    series_labels = ("reference", "candidate"),
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    fig_size = (1100, 600),
    figure_title = "Draft impedance dashboard",
    panel_titles = ("Self resistance", "Self reactance"),
    legend_position = nothing
)
# A figure-scoped legend can be placed in any addon dock and given any native
# Makie Legend attributes.
figurelegend!(
    legend_scope_demo;
    position = :top,
    title = "Result set",
    legend_labels = Dict(
        :result_1 => "baseline",
        :result_2 => "alternative"
    ),
    orientation = :horizontal,
    nbanks = 2,
    overflow = :show_all
)
# Logical position `(1, 1)` is the resistance panel.
panellegend!(
    legend_scope_demo,
    (1, 1);
    position = :inside,
    anchor = :lb,
    title = "Resistance result",
    legend_labels = ("base R", "alternative R"),
    backgroundcolor = (:white, 0.92),
    overflow = :show_all
)
# Titles use the same logical panel address as panel legends.
figuretitle!(legend_scope_demo, "Series impedance dashboard"; fontsize = 20)
paneltitle!(legend_scope_demo, (1, 2), "Reactance matrix")
legend_scope_demo.figure #hide

# The same panel legends can be requested at construction time with
# `panel_legends=Dict((1, 1) => (position=:inside, anchor=:lb,
# legend_labels=("base R", "alternative R")))`. At construction, prefer
# `series_labels`; `legend_labels` remains a compatibility alias. Runtime
# dictionaries target stable source keys such as `:result_1`. These labels do
# not alter result tensors or observation coordinates.

# ### Derived R/L/G/C dashboard

# `L` and `C` are not aliases for stored arrays: the publication layer derives
# them from reactance or susceptance and angular frequency, rejects zero
# frequency where the quantities are undefined, and assigns their real physical
# units. For this 2×2 result, `layout=(2,2)` requests one matrix dashboard per
# physical quantity.

## `layout` fixes the requested panel grid; it does not restrict later Makie edits.
rlgc_dashboards = Makie.plot(
    parameters,
    (R, L, G, C);
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    layout = (2, 2),
    fig_size = (950, 700),
    legend_position = nothing
)
rlgc_dashboards[1].figure #hide
rlgc_dashboards[2].figure #hide
rlgc_dashboards[3].figure #hide
rlgc_dashboards[4].figure #hide

# Changing `layout` never merges or reinterprets physical quantities. `(1,1)`
# would instead return one figure for every quantity/coordinate facet.

# ### Polar impedance and admittance

# Function selectors are resolved through the same observation grammar.
# `abs` and `angle` publish magnitude and phase for both ``Z`` and ``Y``; the
# Makie extension does not infer those transforms from plot labels.

## The four transformed quantities become four matrix dashboards.
polar_dashboards = Makie.plot(
    parameters,
    (abs, angle);
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    layout = (2, 2),
    fig_size = (950, 700),
    legend_position = nothing
)
polar_dashboards[1].figure #hide
polar_dashboards[2].figure #hide
polar_dashboards[3].figure #hide
polar_dashboards[4].figure #hide

# `real`, `imag`, `Z`, or `Y` are valid selectors as well. New scientific
# transforms belong in the observation/publication grammar; plot-specific
# transforms should not be hidden in the Makie extension.

# ### Standalone series-impedance result

# `SeriesImpedance` can be plotted before it is bundled into `LineParameters`.
# Because that object does not own a frequency vector, the frequency samples are
# an explicit positional argument. The rest of the native dashboard behavior is
# identical.

standalone_impedance = Makie.plot(
    parameters.Z,
    parameters.f,
    (Z, 1, 1, Colon());
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    fig_size = (900, 440),
    legend_overflow = :show_all
)
standalone_impedance.figure #hide

# ### Standalone shunt-admittance result

# `ShuntAdmittance` follows the same rule: frequencies are explicit and its
# default scientific family is `(G, B)`.

standalone_admittance = Makie.plot(
    parameters.Y,
    parameters.f,
    (Y, 1, 1, Colon());
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    fig_size = (900, 440),
    legend_overflow = :show_all
)
standalone_admittance.figure #hide

# ### Exact observation requests and modal coordinates

# An `@observe` request is the precise extension point for matrix coordinates.
# It passes through `Grammar.observation_request`, so selection rules stay shared
# by plotting, tables, and reports. This example keeps only one diagonal entry
# from each quantity.

resistance_request = @observe R[1, 1, :];
inductance_request = @observe L[1, 1, :];
selected_response = Makie.plot(
    parameters,
    (resistance_request, inductance_request);
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    layout = (1, 2),
    fig_size = (900, 420),
    legend_overflow = :show_all
)
selected_response.figure #hide

# ### Intent inference and a stacked layout

# A complex impedance request cannot be drawn honestly on one real-valued axis.
# The plotting adapter therefore expands `Z[1,1,:]` into its `R` and `X`
# observables and infers two axes. `layout=(2,1)` only stacks those already
# inferred axes; `(1,2)` would place the identical data side by side.

self_impedance_request = @observe Z[1, 1, :];
stacked_self_impedance = Makie.plot(
    parameters,
    self_impedance_request;
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    layout = (2, 1),              # two inferred quantities, one column
    panel_titles = ("Self resistance", "Self reactance"),
    legend_position = :inside,    # native overlay; use any outer dock instead
    legend_anchor = :rt,
    legend_title = "Result set",
    series_labels = ("solution",), # explicitly opt into a one-source legend
    legend_overflow = :show_all,
    fig_size = (720, 720)
)
stacked_self_impedance.figure #hide

# Select more coordinates to add subplots, or more quantities to add pages or
# paired axes. Overlaying curves means passing more result containers, not
# smuggling matrix coordinates into a legend. No layout changes request meaning.

# Modal results retain their physical domain. A diagonal request is therefore
# explicit rather than being smuggled through an `i`, `j`, or plotting adapter
# keyword.

modal_parameters = compute(
    ModalTransformationProblem(parameters),
    ModalTransformationFormulation(:Fortescue; tolerance = 1.0)
);
modal_inductance = Makie.plot(
    modal_parameters,
    (@observe((L, diag)[:, :]),);
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    fig_size = (820, 430),
    legend_overflow = :show_all
)
modal_inductance.figure #hide

# Extend coordinate behavior in the grammar/domain owner. Extend only visual
# appearance by mutating the returned axes or forwarding Makie plot attributes.

# ### Measurement uncertainty

# If the published scalar type carries uncertainty, `_addon_line!` draws native
# Makie error bars around the nominal line. Limits include the uncertainty
# bounds, so error bars are not clipped by a nominal-only autolimit pass.

measured_parameters = LineParameters(
    complex.(
        measurement.(resistance, 0.05 .* resistance),
        measurement.(inductance .* angular_frequency,
            0.05 .* inductance .* angular_frequency)
    ),
    complex.(
        measurement.(conductance, abs.(0.05 .* conductance)),
        measurement.(capacitance .* angular_frequency,
            abs.(0.05 .* capacitance .* angular_frequency))
    ),
    frequency
);
uncertainty_plot = Makie.plot(
    measured_parameters,
    (@observe(R[1, 1, :]), @observe(L[1, 1, :]));
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    fig_size = (900, 440),
    legend_overflow = :show_all
)
uncertainty_plot.figure #hide

# A new uncertainty scalar owner extends nominal/error extraction. The plotting
# code remains unchanged and should not inspect a package-specific storage type.

# ### Comparing completed line results

# A named tuple supplies both completed results and default legend labels. The
# comparison preparation validates compatible matrix shape, basis, and domain,
# then overlays sources in one axis for every selected matrix position. Each
# result may retain its own frequency samples.

candidate = LineParameters(
    parameters.Z.values .* (1.06 + 0im),
    parameters.Y.values .* (0.94 + 0im),
    parameters.f
);
comparison_plot = Makie.plot(
    (; reference = parameters, candidate),
    (R,);
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    layout = (2, 2),
    fig_size = (900, 700),
    legend_position = :top,
    legend_attributes = (; orientation = :horizontal),
    legend_overflow = :show_all
)
comparison_plot.figure #hide

# The exact-request positional form is
# `Makie.plot(reference, candidate, @observe(Z[1,1,:]);
# series_labels=("reference", "candidate"))`.
# Change line styling after construction through the native plot objects in each
# axis; change source identity through `series_labels` or named-tuple keys. The
# older `legend_labels` and `legend` keywords remain compatibility aliases.

# ### Generic observation publication

# The generic publication plot is the fallback for detached scientific
# observations. Each distinct published quantity becomes one axis. A published
# frequency column is the abscissa, while row/column coordinates define series;
# only publications with no domain coordinate fall back to sample index. This is
# also the path used by report illustrations.

publication = observables(
    parameters,
    (resistance_request, inductance_request)
);
publication_plot = Makie.plot(
    publication;
    title = "Published diagonal observations",
    figure_title = "Published diagonal observations",
    panel_titles = ("R[1,1]", "L[1,1]"),
    backend = :cairo,
    display_plot = false,
    controls = false,
    layout = (1, 2),
    fig_size = (900, 480),
    display_legend = true,
    legend_title = "Published quantity",
    ## R and L share the same Z[1,1] coordinate group across both panels.
    legend_labels = ("self impedance",),
    legend_position = :bottom,
    legend_attributes = (; orientation = :horizontal),
    legend_overflow = :show_all
)
publication_plot.figure #hide

# Add an owner-specific recipe only when a publication needs visual semantics
# beyond its declared quantity and coordinate columns.

# ## Geometry preview recipes

# ### Cable-design cross-section

# `DataModel.preview_shapes` exposes detached physical polygons with only their
# material and construction tag. The Makie preview adapter derives optional
# presentation groups and `_addon_preview_axis!` draws the polygons with native
# `poly!`, locks the axis to `DataAspect`, and computes
# geometry limits. The aspect canvas and its dock are centered as one responsive
# group.

## `display_id=true` promotes the design identifier into the panel title.
## Legend and colorbars share the right dock without replacing one another.
cable_preview = preview(
    mv_design;
    backend = :cairo,
    display_plot = false,
    controls = false,
    display_id = true,
    size = (950, 700),
    ## Group the physical `core_strands_*` regions into one presentation entry.
    ## Geometry tags, terminals, and electrical construction remain untouched.
    legend_group = region -> startswith(
        String(region.source.tag), "core_strands_") ?
        :stranded_core : region.source.tag,
    legend_labels = Dict(:stranded_core => "Stranded core"),
    legend_position = :right,
    legend_attributes = (; nbanks = 2),
    colorbar_position = :right,
    colorbar_attributes = (; vertical = false)
)
cable_preview.figure #hide

# `legend_group` accepts either the shown callback over a `PlacedRegion` or a
# tag-to-group dictionary. `legend_labels` maps the resulting presentation group
# to text. This keeps display grouping independent of physical tags. For detailed
# annotation, mutate `cable_preview.axes[1]` and add ordinary Makie plots. Each
# object in `cable_preview.colorbars` is a native `Colorbar`.

# ### Cable-design collection

# A collection preview repeats the same detached geometry path for every design
# and assigns the caller-requested layout. The material ranges are aggregated
# once, so every panel uses comparable colors and one shared set of scales.

design_collection = preview(
    [mv_design, hv_design, hv_design, mv_design];
    layout = (2, 2),
    backend = :cairo,
    display_plot = false,
    controls = false,
    size = (1000, 850),
    figure_title = "Cable design family",
    panel_titles = ("MV option A", "HV option A", "HV option B", "MV option B"),
    colorbar_position = :bottom,
    colorbar_attributes = (; vertical = false)
)
design_collection.figure #hide

# Use any sufficient `(rows, columns)` layout. The default omits layer legends;
# request selected local legends with `panel_legends`, using the same logical
# grid-position contract as line dashboards.

# ### Cable-system cross-section

# A system preview resolves every placed region into the system frame, adds
# reference geometry such as the earth interface, and derives limits from the
# physical placement or `zoom_factor`. Earth properties use the same atomic
# color-scheme contract as cable materials.

system_preview = preview(
    cable_system;
    earth_model = earth,
    zoom_factor = 1.35,
    backend = :cairo,
    display_plot = false,
    controls = false,
    display_id = true,
    size = (1000, 700),
    legend_position = :right,
    colorbar_position = :right,
    colorbar_attributes = (; vertical = false)
)
system_preview.figure #hide

# Modify physical contents before previewing; modify visual annotations after
# previewing. `zoom_factor` changes only the initial view, and the reset control
# returns to those computed limits.

# ## Material colors and native colorbars

# ### One reusable material scheme

# A color scheme is the atom. `material_property_ranges` obtains values from a
# design, a design collection, or the material defaults. `materialcolors`
# turns exactly one property and range into a Makie-compatible named tuple:
# `label`, `colormap`, `limits`, and `ticks`. `materialscale!` places exactly
# that one scheme wherever the caller chooses.

material_ranges = LineCableModels.DataModel.material_property_ranges(mv_design);
rho_scheme = materialcolors(
    :rho,
    material_ranges.rho
);
rho_scale_figure = Figure(size = (850, 180), figure_padding = 24)
materialscale!(
    rho_scale_figure[1, 1],
    rho_scheme;
    vertical = false,
    width = Relative(0.85)
)
rho_scale_figure

# The scheme contains no placement. Put that `Colorbar` in any `GridPosition`,
# combine it with a heatmap, or reuse the same scheme in another figure. Define a
# new property by defining another palette producer next to `materialcolors`;
# physical range collection remains a DataModel concern.

# ### High-level material-scale reference

# `show_material_scale` is preview sugar that privately composes the three
# default property schemes. Its result still exposes three independent native
# colorbars; it is not the reusable unit of the color API.

material_scale = show_material_scale(
    backend = :cairo,
    display_plot = false,
    controls = false,
    size = (850, 360),
    figure_title = "Reusable material-property schemes",
    colorbar_attributes = (; vertical = false)
)
material_scale.figure #hide

# For one property, use the preceding `materialscale!(position, scheme)` pattern.
# For a different high-level combination, compose schemes in the caller's own
# Makie layout.

# ## Monte Carlo result recipes

# Monte Carlo methods first publish the requested marginal from retained
# samples and/or `HistogramDensity`. Cable-constant requests accept `R`, `L`,
# `C`, `G`, or `(selector, assembly)`. Matrix-valued line results require an
# exact request such as `@observe R[1, 1, 3]`. Native Makie function identity
# chooses the visual primitive.

# ### Sample histogram

# `Makie.hist` uses retained samples and calls native `hist!`. `bins` and
# `normalization` retain their Makie meanings; unit options are consumed while
# publishing the marginal.

sample_histogram = Makie.hist(
    mc_result,
    R;
    bins = 12,
    normalization = :pdf,
    backend = :cairo,
    display_plot = false,
    controls = false,
    fig_size = (820, 400),
    figure_title = "Retained Monte Carlo samples",
    legend_overflow = :show_all,
    color = :steelblue
)
sample_histogram.figure #hide

# Pass ordinary `hist!` attributes such as `color`, `strokewidth`, or
# `transparency` through the remaining keywords. Use `(R, assembly)` when a
# `CableConstants` result contains more than one assembly.

# ### Retained-model probability density

# `Makie.stairs` consumes the retained histogram model and calls native
# `stairs!` with post-step edges. It does not regenerate a stochastic result.

model_density = Makie.stairs(
    mc_result,
    R;
    backend = :cairo,
    display_plot = false,
    controls = false,
    fig_size = (820, 400),
    legend_overflow = :show_all,
    color = :darkorange,
    linewidth = 3
)
model_density.figure #hide

# Change the model upstream or select another retained marginal to change the
# data. Use native `stairs!` keywords for appearance.

# ### Empirical cumulative distribution

# `Makie.ecdfplot` consumes retained samples and delegates the empirical curve
# to Makie's `ecdfplot!` recipe.

empirical_cdf = Makie.ecdfplot(
    mc_result,
    R;
    backend = :cairo,
    display_plot = false,
    controls = false,
    fig_size = (820, 400),
    legend_overflow = :show_all,
    color = :seagreen,
    linewidth = 3
)
empirical_cdf.figure #hide

# The returned axis can be combined with confidence bands or additional native
# curves; the addon only owns the initial empirical series and its legend group.

# ### Retained-model cumulative distribution

# `Makie.lines` evaluates the retained histogram model through the UQ owner's
# cumulative-probability function and draws the resulting grid with `lines!`.

model_cdf = Makie.lines(
    mc_result,
    R;
    backend = :cairo,
    display_plot = false,
    controls = false,
    fig_size = (820, 400),
    legend_overflow = :show_all,
    color = :firebrick,
    linewidth = 3
)
model_cdf.figure #hide

# Increase visual resolution by extending the owner-side model grid policy; add
# purely visual reference curves directly to `model_cdf.axes[1]`.

# ### Sample/model Q-Q plot

# `Makie.qqplot` requests both retained products, lets the UQ owner calculate
# matching quantile pairs, draws native scatter points, and optionally adds the
# identity reference line.

quantile_plot = Makie.qqplot(
    mc_result,
    R;
    qqline = :identity,
    backend = :cairo,
    display_plot = false,
    controls = false,
    fig_size = (820, 440),
    panel_titles = ("Sample versus retained model",),
    legend_title = "Q–Q elements",
    legend_labels = ("sample quantiles", "identity reference"),
    legend_position = :inside,
    legend_anchor = :lt,
    legend_attributes = (; backgroundcolor = (:white, 0.92)),
    legend_overflow = :show_all,
    color = :purple,
    markersize = 10
)
quantile_plot.figure #hide

# Set `qqline=:none` to remove the reference. Other scatter attributes are
# forwarded to `scatter!`; data and units remain publication concerns.

# ## Native composition with the addon shell

# ### A caller-owned 2×2 plot

# `plotwindow` is the escape hatch when no high-level recipe is appropriate. It
# creates the common shell and passes only its content `GridLayout` to the
# callback. The callback uses normal Makie constructors; afterward the shell
# discovers the native axes and attaches reset and export controls.

custom_dashboard = LineCableModels.plotwindow(
    title = "Caller-owned diagnostics",
    figure_title = "Four caller-owned Makie axes",
    size = (900, 650),
    layout = (2, 2),
    backend = :cairo,
    display_plot = false,
    controls = false
) do grid
    for row in 1:2, column in 1:2
        axis = Axis(
            grid[row, column];
            title = "Response $row,$column",
            xlabel = "Frequency [Hz]",
            ylabel = "Amplitude"
        )
        lines!(
            axis,
            frequency,
            @. (row + column) * sin(log10(frequency));
            color = Makie.wong_colors()[(row - 1) * 2 + column],
            linewidth = 2
        )
        axis.xscale = log10
    end
end
custom_dashboard.figure #hide

# Nothing prevents nested layouts, `Axis3`, `Colorbar`, `Slider`, custom Makie
# recipes, or arbitrary plot primitives in the callback. If a native composition
# needs no LineCableModels controls or export policy, use `Figure` directly.

# ## What the addons do

# ### Scientific axes, scale controls, and limits

# Each managed scientific axis receives labels from its published quantity and
# unit. On a linear axis, large or small values are divided by a power of ten and
# that exponent is rendered once in the axis label; tick labels remain readable.
# Logarithmic axes instead receive explicit decade ticks. The x/y toggles mutate
# the axes' native `scale` observables and run the same reset/format callbacks.
#
# Limits are calculated from finite visible data, including measurement error
# bounds. Constant series receive nonzero padding. Legend visibility changes
# trigger another limit pass, so hiding a dominant curve exposes the remaining
# data instead of leaving a stale range.

# ### Legends and docks

# `legend_position` selects only placement. `:inside` overlays the union of the
# figure's axis viewports or the single axis viewport for a panel legend;
# `legend_anchor` selects the corner or center. Named and grid-position docks
# remain outside the plot area. `legend_attributes` is merged into the native
# `Legend` constructor, so orientation, bank count, padding, background,
# alignment, and other Makie options remain available. `legend_overflow` is the
# one addon policy: `:ellipsis` fits entries to the current bounding box and
# restores them when space returns; `:show_all` always retains all entries.
# `figurelegend!` and `panellegend!` rebuild native legends from the retained
# semantic handle registry, so placement, title, and semantic labels can also be
# changed after construction. Clicking/toggling a grouped Makie legend entry
# continues to affect every plot handle in that group.
#
# High-level recipes may put a legend and colorbars in the same dock. The shell
# creates a nested `GridLayout` and gives each block its own cell. Resizing moves
# the entire canvas/dock group; it does not anchor a block to the raw window
# while stretching only the outer figure.

# ### Colorbars

# `colorbar_position` and `colorbar_attributes` mirror the legend placement
# contract. A preview can privately request several schemes, but
# `_addon_colorbar!` always consumes one scheme and creates one native
# `Colorbar`. The reusable public atom remains
# `materialcolors(property, range)` plus `materialscale!(position, scheme)`.

# ### Observables and ownership

# Two distinct concepts are intentionally present. `@observe` belongs to the
# scientific observation grammar: it selects a physically meaningful view of a
# result vault, including matrix and frequency coordinates. Makie's
# `Observable` type drives live UI state such as controls, scales, limits,
# visibility, and layout bounds. The adapter consumes the former and wires the
# latter; it does not replace either one with an `AxisBehavior` aggregate.
# Keeping the `UIPlot` alive keeps the native figure and its subscriptions alive.
#
# Native mutation is therefore the normal extension mechanism:

owned_plot = Makie.plot(
    parameters,
    @observe R[1, 1, :];
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    fig_size = (820, 400),
    legend_overflow = :show_all
);
owned_axis = only(owned_plot.axes);
owned_axis.title[] = "Caller-owned resistance";
vlines!(owned_axis, [100.0, 1_000.0]; color = :black, linestyle = :dash);
owned_plot.figure #hide

# ### Responsive layout

# Inferred two-component pages can reflow between side-by-side and stacked axes
# as a live window changes aspect. Explicit layouts and matrix dashboards retain
# their semantic grid coordinates. GL windows also retain a scientific minimum
# size. Data-aspect previews size the square canvas first, place
# legend/colorbar docks next to it, and center the resulting group.

# ### Current-state SVG export

# [`export_svg`](@ref) saves the current live figure through CairoMakie. For a
# publication export it temporarily hides the toolbar and status row, switches
# the figure's font roles to Makie's LaTeX font theme, uses a white background,
# and then restores every changed observable and the previously active backend.
# Caller-added plots, visibility, scales, limits, and annotations are saved
# because export does not reconstruct an earlier specification.

export_directory = mktempdir();
export_svg(
    owned_plot;
    path = joinpath(export_directory, "caller_owned_resistance.svg"),
    theme = :publication,
    open_file = false
); #hide

# ## Adding or changing a managed recipe

# A new high-level recipe should preserve the same ownership boundary:
#
# 1. Expose numerical observations or physical geometry through the owning
#    module's existing public protocol. Do not add plot preparation, labels,
#    colors, layout, or Makie types to scientific owners.
# 2. Add request normalization and the narrow public dispatch method in
#    `LineCableModelsMakieExt`, then call
#    native Makie constructors or primitives there.
# 3. Reuse only the addon services the recipe needs: `_addon_shell`,
#    `_addon_axis!`, `_addon_finish!`, or `plotwindow`. Do not create a second
#    plot specification or an optional adapter hierarchy.
# 4. Return `UIPlot` with the actual native objects and leave further mutation to
#    the caller.
# 5. Add the real call to this literate gallery, a Cairo rendering assertion,
#    and an interactive GL inspection fixture when resizing or widgets matter.
#
# A purely visual variation usually needs no new managed recipe: pass a native
# attribute, mutate the returned block, add a Makie primitive, or start from
# `plotwindow`. A new recipe is justified when LineCableModels owns meaningful
# publication, units, coordinate selection, initial limits, or a reusable piece
# of scientific interaction.
