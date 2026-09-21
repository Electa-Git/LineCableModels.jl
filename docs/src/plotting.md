```@meta
EditURL = "../literate/plotting.jl"
```

# Makie plotting: implementation guide and complete gallery

Plot cable geometry, electrical parameters, and uncertainty results with Makie.
Each example shows the call that produces its figure. The returned `UIPlot`
contains native Makie objects that you can modify directly.

````@example plotting
using LineCableModels
using CairoMakie
using LinearAlgebra: diag
using Measurements: measurement
using Statistics: mean
````

## Figures, titles, and legends

A call that produces one figure returns a [`UIPlot`](@ref); a call that
produces several figures returns `Vector{UIPlot}`. The `figure`, `axes`,
`controls`, `legend`, `panel_legends`, and `colorbars` fields expose the
displayed Makie objects.

### Naming titles and legends

Layout text has one name per scope. These keywords are presentation metadata;
none of them changes a result tensor or geometry tag.

| Keyword | Scope |
|:--|:--|
| `figure_title` | One visible title above the whole native figure |
| `title_attributes` | Native Makie `Label` attributes for `figure_title` |
| `panel_titles` | Axis-title overrides, positionally or by semantic key |
| `legend_title` | Heading of the controlled figure legend |
| `series_labels` | Names of overlaid result containers; these are the legend entries |
| `series_attributes` | Native Makie attributes for all legend groups, or one named tuple per group in legend order |
| `legend_position` | `:inside`, a named outer dock, or a positive dock grid position |
| `legend_attributes` | Native `halign`/`valign`, orientation, banks, fonts, and padding |

Some established recipes also accept `title` as their window/export or
single-recipe heading. Use the explicit scoped names above when composing a
dashboard.
Benchmark windows default to `case ID — quantity`, followed by block indices
when split. Deterministic, mean ± std, and explicit-statistic comparisons use
this same rule. An explicit `title` overrides the default case prefix; subplot
titles and `figure_title` remain independent.

`series_attributes` uses the shared PlotBuilder controls for matrix and benchmark
plots, observed results, statistical plots, and geometry previews.
A named tuple applies to every group; a tuple or vector of named tuples styles
each group separately. For example,
`series_attributes=((marker=:circle, markersize=8), (;), (linestyle=:dash,))`
adds markers to the first series, keeps the second's defaults, and dashes the
third. Styles apply across facets and pages, including their legends and
visibility controls. Attributes must be supported by the group's native plots.
In `plotwindow`, groups follow native plot insertion order across its axes.

The same scopes are mutable after construction. `figuretitle!` and
`paneltitle!` replace titles. `figurelegend!` and `panellegend!` rebuild a
native legend from the retained semantic plot handles, so relabeling does not
destroy grouped visibility behavior.

### Observation, facet, page, and legend semantics

`LineParameters` owns dense ``Z`` and ``Y`` tensors. The observation grammar
is what turns those tensors into displayable physical quantities: ``Z``
expands to ``R`` then ``X``; ``Y`` expands to ``G`` then ``B``. An exact
request such as `@observe Z[1,1,:]` keeps its conductor coordinates while
selecting the complete frequency range. Plotting consumes that resolved
request; it does not ask the caller to repeat `(R, X)`.
The plot-facing name for this ordinate selection is `ydata`; it may be passed
positionally or as a keyword, for example `plot(result; ydata=(R, L))`.

Every `(quantity, row, column)` is one facet and therefore one axis. Its
default title comes from the unit registry plus the coordinate relation, for
example `Self series resistance — conductor 1` or
`Mutual series reactance — conductor 1 → 2`. Coordinates never
become legend entries. A legend exists to distinguish overlaid result sets;
a single unlabeled result therefore has no legend by default.

`layout` is the sole nominal panel capacity. Every selected quantity or
statistical meaning has its own figure family; multiple observations add traces.
Before: `blocks` paginated matrices and `layout` could pair quantities. Now
`layout` alone sets capacity within each quantity; `blocks` has been removed.

| `layout` | Page grouping |
|:--|:--|
| `nothing` | Full declared matrix extent, or a near-square flow capacity |
| `(1, 1)` | One selected coefficient per figure, per quantity |
| `(1, 2)` | Up to two matrix columns or two flow panels per page |
| `(2, 1)` | Up to two matrix rows or two flow panels per page |
| `(N, N)` | One full matrix page per quantity for an `N×N` result |

Layout follows observation selection. It cannot add an excluded coefficient,
merge different quantities, or turn coordinates into series labels.

### Matrix pagination and overlays

`plot(results; ydata=(R,L), layout=(2,2))` produces four R pages and four L
pages for a full 3×3 matrix. Their actual extents are 2×2, 2×1, 1×2, and 1×1.
Pages follow request order and then row-major matrix-page order. Original
coefficient identities remain the addresses for titles, legends, scales, and reset.
Positional `panel_titles` bind to the complete selected population before paging.

Residual pages contain only actual domain tracks. In-domain selection holes
remain part of the matrix frame region. Native decoration measurements establish
equal initial data frames; residual outer windows fit compactly around them.
`fig_size` describes nominal capacity and `figure.size` takes precedence. Later
manual resizing affects only that page. There is no window aspect lock.
Explicit diagonal observations and preview collections use the same capacity,
preserving their source order and original identities in compact flow pages.

Result overlays use solid lines with sparse, staggered markers on saved
sample points: black curves and hollow circles for references, colored curves
and filled shapes for candidates. References may themselves carry uncertainty.
Default routes and explicit implementations have the same styling semantics.
Deterministic references mark both endpoints.
Colors and marker identities remain stable
when formulations are filtered. No curves are merged because they agree.
Use `series_attributes=(marker=nothing,)` for lines only, or an explicit native
`marker` for all-sample placement. References remain comparison operands, not
declarations of physical truth.

With multiple displayed series, `errorbar_sampling=:staggered` selects sparse error
bars at retained samples separately from automatic markers. Both references
and candidates retain their full mean curves. Coordinates, uncertainties,
comparison calculations, and full-data axis limits are unchanged. X and Y
intervals use the same indices. When very few samples are available, intervals
take priority over conflicting automatic markers; the legend retains identity.
Use `errorbar_sampling=:all` to inspect every interval, with automatic markers
omitted on uncertain series. A single displayed series uses `:all` by default.
Explicit native markers still use
all samples. Native `whiskerwidth` and `linewidth` overrides take priority.
Sparse glyphs are an overview: use full intervals or separate standard-deviation
curves to inspect uncertainty variation. Mean ± std is not mean ± standard error.
Many exactly coincident methods cannot all remain distinguishable at finite
screen resolution; use legend visibility to inspect them separately.

Within one identified physical point, observation-owned groups may share a trace
for equivalent quantity-relevant selections: impedance choices on Z/R/L/X pages,
admittance choices on Y/G/C/B pages. Every relevant composite route, control,
coordinate and uncertainty meaning participates; equal curves or descriptions
alone never merge cases. Different physical points remain separate. Conflicting
observations under the same selection raise an error. Saved results remain intact.
Default labels show only relevant differences across candidates and reference,
omitting common physical inputs and individual controls. Names are captured from
owner-dispatched `description(...; compact=true)` methods used by report tables:
`FEM (reference)`, `PSCAD (reference)`, `Monte Carlo (reference)`, and candidates
such as `LEP` or `earth Z=Saad`, without candidate numbering.
`formulations=[3,1]` selects recorded original formulation identities, preserving
order and colors. Explicit `series_labels` and styles follow the retained source;
they do not create new formula identities. When nothing varies, labels use the
applicable compact owner description. Automatic candidates are chromatic;
the separate reference is black and does not change candidate colors.
Full scientific explanations and common settings remain in `formula_details`.

Top/bottom legends fit a measured row-major grid and wrap long labels without
dropping fields. Explicit `legend_attributes=(orientation=..., nbanks=...)`
retains native manual layout control.

## Gallery data

The frequency responses are deliberately small but non-constant so that
logarithmic axes, engineering units, legends, and automatic limits are all
visible in the generated documentation.

````@example plotting
frequency = collect(10.0 .^ range(1, 4; length = 24));
angular_frequency = reshape(2π .* frequency, 1, 1, :);
resistance = cat(
    (
        [1.0 0.22; 0.22 1.8] .* 1.0e-4 .* (1 + 0.12 * log10(f / first(frequency)))
    for f in frequency
    )...;
    dims = 3);
inductance = cat(
    (
        [2.0 0.28; 0.28 2.5] .* 1.0e-7 .* (1 - 0.04 * log10(f / first(frequency)))
    for f in frequency
    )...;
    dims = 3);
conductance = cat(
    (
        [3.0 -0.45; -0.45 4.0] .* 1.0e-9 .* (1 + 0.08 * log10(f / first(frequency)))
    for f in frequency
    )...;
    dims = 3);
capacitance = repeat([4.0 -0.7; -0.7 5.0] .* 1.0e-10, 1, 1, length(frequency));
parameters = LineParameters(
    complex.(resistance, inductance .* angular_frequency),
    complex.(conductance, capacitance .* angular_frequency),
    frequency
);
nothing #hide
````

The geometry gallery uses the example library shipped with the repository.
A system is assembled from two placements without introducing a plotting-only
geometry representation.

````@example plotting
cable_library = CablesLibrary();
LineCableModels.load!(
    cable_library;
    file_name = joinpath(pkgdir(LineCableModels), "examples", "cables_library.json")
);
mv_design = cable_library["18kV_1000mm2"];
hv_design = cable_library["525kV_1600mm2"];
earth = EarthModel(100.0, 10.0, 1.0);
cable_system = build(
    LineCableSystem,
    [mv_design, mv_design],
    [(-0.06, -0.20), (0.06, -0.20)];
    environment = earth,
    system_id = "two-cable-gallery",
    line_length = 1_000.0
);
nothing #hide
````

A small synthetic retained fixture demonstrates the completed UQ storage used
below. It does not run a Monte Carlo campaign or claim a physical cable study.
Moments, samples, and histograms belong to the UQ owner before plotting begins.

````@example plotting
retained_samples = (
    R=reshape([2.,3.,5.,8.].*1e-4,1,:),
    L=reshape([11.,13.,17.,19.].*1e-7,1,:),
    C=reshape([23.,29.,31.,37.].*1e-11,1,:),
    G=reshape([41.,43.,47.,53.].*1e-10,1,:),
);
retained_statistics = map(x -> [SampleSummary(vec(x))], retained_samples);
retained_histograms = map(x -> [HistogramDensity(vec(x);bins=2)], retained_samples);
retained_core = CableConstants(mean(retained_samples.R),mean(retained_samples.L),
    mean(retained_samples.C),mean(retained_samples.G));
retained_core = LineCableModels.materialize(retained_core,retained_statistics);
mc_formulation = MonteCarlo(Formulation();trials=4,seed=41,
    return_samples=true,return_histograms=true);
mc_result = MonteCarloResult(mc_formulation,[retained_core],[retained_statistics],
    [retained_samples],[retained_histograms],UInt64(41),UInt64[42],[4]);
nothing #hide
````

## Line-parameter recipes

### Complete default view

`Makie.plot(parameters)` is the minimal call. It observes everything in the
order ``Z`` then ``Y``, expands those families to ``R``, ``X``,
``G``, and ``B``, and returns four matrix-dashboard pages. For this 2×2
result, every page contains four axes. The calls below render the full default
gallery in exactly that order.

````@example plotting
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
````

To change only appearance, mutate the returned native axes. To change which
physical values exist, pass an observation selector or exact `@observe`
request. To change pagination, pass one of the layouts in the table above.

### Cartesian series impedance

This exact observation asks for the self impedance of conductor 1 over the
full frequency range. The observation layer expands ``Z`` into resistance and
reactance. They remain separate quantity figures, even with a multi-panel layout.

````@example plotting
# `xscale=:log10` is an initial state; the live toolbar can still change it.
# `figure_title` labels each figure; `panel_titles` follows the selected quantity order.
# `series_labels` names result sets, so supplying one explicitly opts into a legend.
# `:inside` overlays the logical frame region with native alignment attributes.
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
    legend_attributes = (; halign=:right,valign=:bottom,backgroundcolor=(:white,0.92)),
    legend_overflow = :show_all
)
series_cartesian[1].figure #hide
````

````@example plotting
series_cartesian[2].figure #hide
````

To modify this recipe, select other matrix coordinates with an explicit
observation request, pass native legend attributes, or mutate either returned
axis. For example, `series_cartesian[1].axes[1].title[] = "Measured resistance"`
changes the live title without asking LineCableModels to rebuild anything.

### Cartesian shunt admittance

Conductance and susceptance are the real and imaginary parts of ``Y``. They
use the same dashboard implementation, retained units, limits, legend
groups, and controls as impedance; only the scientific requests differ.

````@example plotting
# A bottom dock is horizontal by default and remains a native Makie Legend.
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
shunt_cartesian[1].figure #hide
````

````@example plotting
shunt_cartesian[2].figure #hide
````

Any accepted Makie `Legend` keyword belongs in `legend_attributes`. For a
labeled source or comparison, moving the
block uses `legend_position = :inside`, `:left`, `:right`, `:top`, `:bottom`,
or a positive `(row, column)` dock coordinate. `(2, 2)` remains the plot
canvas. Native `halign` and `valign` also position inside legends.

### Figure-wide and panel-scoped legends

A controlled comparison registers one group per result set, then exposes two
views of those groups. The figure legend spans every panel. A panel legend
filters the same source handles to one logical plot position. Hiding a global
result entry changes that source's visibility across both resistance coefficients.

````@example plotting
legend_scope_demo = Makie.plot(
    parameters,
    parameters,
    @observe R[1, 1:2, :];
    series_labels = ("reference", "candidate"),
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    fig_size = (1100, 600),
    figure_title = "Resistance dashboard",
    panel_titles = ("Self resistance", "Mutual resistance"),
    legend_position = nothing
)
````

A figure-scoped legend can be placed in any addon dock and given any native
Makie Legend attributes.

````@example plotting
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
````

Logical position `(1, 1)` is the resistance panel.

````@example plotting
panellegend!(
    legend_scope_demo,
    (1, 1);
    position = :inside,
    halign=:left,valign=:bottom,
    title = "Resistance result",
    legend_labels = ("base R", "alternative R"),
    backgroundcolor = (:white, 0.92),
    overflow = :show_all
)
````

Titles use the same logical panel address as panel legends.

````@example plotting
figuretitle!(legend_scope_demo, "Resistance dashboard"; fontsize = 20)
paneltitle!(legend_scope_demo, (1, 2), "Mutual resistance")
legend_scope_demo.figure #hide
````

The same panel legends can be requested at construction time with
`panel_legends=Dict((1, 1) => (position=:inside, halign=:left, valign=:bottom,
legend_labels=("base R", "alternative R")))`. Runtime dictionaries target
stable source keys such as `:result_1`. These labels do not alter result
tensors or observation coordinates.

### Derived R/L/G/C dashboard

`L` and `C` are observation-owned proxies derived from reactance or susceptance
and angular frequency. Their DC samples are retained as unavailable; independently
valid components remain available. The observation owner assigns their physical
units. For this 2×2 result, `layout=(2,2)` requests one matrix dashboard per
physical quantity.

````@example plotting
# `layout` fixes the requested panel grid; it does not restrict later Makie edits.
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
````

Changing `layout` never merges or reinterprets physical quantities. `(1,1)`
would instead return one figure for every quantity/coordinate facet.

### Polar impedance and admittance

Function selectors are resolved through the same observation grammar.
`abs` and `angle` publish magnitude and phase for both ``Z`` and ``Y``; the
Makie extension does not infer those transforms from plot labels.

````@example plotting
# The four transformed quantities become four matrix dashboards.
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
````

`real`, `imag`, `Z`, or `Y` are valid selectors as well. New scientific
transforms belong in the observation grammar; plot-specific
transforms should not be hidden in the Makie extension.

### Standalone series-impedance result

`SeriesImpedance` can be plotted before it is bundled into `LineParameters`.
Because that object does not own a frequency vector, the frequency samples are
an explicit positional argument. The rest of the native dashboard behavior is
identical.

````@example plotting
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
standalone_impedance[1].figure #hide
````

````@example plotting
standalone_impedance[2].figure #hide
````

### Standalone shunt-admittance result

`ShuntAdmittance` follows the same rule: frequencies are explicit and its
default scientific family is `(G, B)`.

````@example plotting
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
standalone_admittance[1].figure #hide
````

````@example plotting
standalone_admittance[2].figure #hide
````

### Exact observation requests and modal coordinates

An `@observe` request is the precise extension point for matrix coordinates.
It passes through `Grammar.observation_request`, so selection rules stay shared
by plotting, tables, and reports. This example keeps only one diagonal entry
from each quantity.

````@example plotting
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
selected_response[1].figure #hide
````

````@example plotting
selected_response[2].figure #hide
````

### Family requests and a vertical capacity

The observation owner expands `Z[1,1,:]` into retained R and X products.
`layout=(2,1)` gives each quantity its own nominal two-row capacity; it never
puts different quantities into the same figure. Only selected original panels
are drawn, and no out-of-domain row is allocated.

````@example plotting
self_impedance_request = @observe Z[1, 1, :];
stacked_self_impedance = Makie.plot(
    parameters,
    self_impedance_request;
    backend = :cairo,
    display_plot = false,
    controls = false,
    xscale = :log10,
    layout = (2, 1),              # two-row nominal capacity per quantity
    panel_titles = ("Self resistance", "Self reactance"),
    legend_position = :inside,    # native overlay; use any outer dock instead
    legend_attributes=(halign=:right,valign=:top),
    legend_title = "Result set",
    series_labels = ("solution",), # explicitly opt into a one-source legend
    legend_overflow = :show_all,
    fig_size = (720, 720)
)
stacked_self_impedance[1].figure #hide
````

````@example plotting
stacked_self_impedance[2].figure #hide
````

Select more coordinates to add panels, or more quantities to add figure families. Overlaying curves means passing more result containers, not
smuggling matrix coordinates into a legend. No layout changes request meaning.

Modal results retain their physical domain. A diagonal request is therefore
explicit rather than being smuggled through an `i`, `j`, or plotting adapter
keyword.

````@example plotting
modal_parameters = compute(
    ModalTransformationProblem(parameters),
    ModalTransformationFormulation(:default);
    options = (offdiagonal_tolerance = 1.0,)
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
````

Extend coordinate behavior in the grammar/domain owner. Extend only visual
appearance by mutating the returned axes or forwarding Makie plot attributes.

### Measurement uncertainty

Uncertain retained values draw native intervals around their nominal curves.
Full uncertainty support participates in limits even when intervals are sparse.
`ObservedResult` classifies engineering zero using absolute nominal magnitude
and owner-defined cutoffs. With `clip=true`, it recenters those nominal values
without discarding their uncertainty dependencies. Phase unavailability and
linked X/L or B/C thresholds are observation concerns. Plotting applies no
further clipping and never applies an operand floor to retained RMS errors.
`clip=false` and `atol` are raw acquisition options. Compatible display units
also work directly on retained inputs: `plot(observed; length_unit=:base)`.
Plotting delegates re-expression to ObservedResult; it preserves recorded
clipping decisions and does not rerun comparison or change saved RMS units.

Legend actions hide or restore the nominal line, its markers, and its x/y
error bars together. Figure legends act across the figure; panel legends act
only on their panel. This also applies to overlays and report illustrations.
Native error-bar handles remain independently editable; the next series
hide/show action restores the complete set of components.

````@example plotting
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
uncertainty_plot[1].figure #hide
````

````@example plotting
uncertainty_plot[2].figure #hide
````

A new uncertainty scalar owner extends nominal/error extraction. The plotting
code remains unchanged and should not inspect a package-specific storage type.

### Comparing completed line results

A named tuple supplies both completed results and default legend labels. The
observation owner validates usable units and coordinates. Plotting overlays sources in one axis for every selected matrix position. Each
result may retain its own frequency samples.

````@example plotting
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
````

The exact-request positional form is
`Makie.plot(reference, candidate, @observe(Z[1,1,:]);
series_labels=("reference", "candidate"))`.
Change line styling after construction through the native plot objects in each
axis; change source identity through `series_labels` or named-tuple keys.

### Detached observed results

`ObservedResult` retains complete primary representations and their original
matrix coordinates. The existing matrix renderer selects these records for
display. Report illustrations use the same observed-input plotting method.

````@example plotting
observed = ObservedResult(
    parameters,
    (resistance_request, inductance_request)
);
observed_plot = Makie.plot(
    observed;
    ydata=(resistance_request,inductance_request),
    title = "Retained coefficient observations",
    figure_title = "Retained coefficient observations",
    panel_titles = ("R[1,1]", "L[1,1]"),
    backend = :cairo,
    display_plot = false,
    controls = false,
    layout = (1, 2),
    fig_size = (900, 480),
    legend_title = "Observed result",
    # R and L retain the same original coefficient in separate figures.
    series_labels = ("self impedance",),
    legend_position = :bottom,
    legend_attributes = (; orientation = :horizontal),
    legend_overflow = :show_all
)
observed_plot[1].figure #hide
````

````@example plotting
observed_plot[2].figure #hide
````

A new primary result owner supplies observation methods; the same retained
quantity records reuse table and plot consumption.

## Geometry preview recipes

### Cable-design cross-section

`DataModel.preview_shapes` exposes detached physical polygons with only their
material and construction tag. The Makie preview adapter derives optional
presentation groups and `_addon_preview_axis!` draws the polygons with native
`poly!`, locks the axis to `DataAspect`, and computes
geometry limits. The aspect canvas and its dock are centered as one responsive
group.

````@example plotting
# `display_id=true` promotes the design identifier into the panel title.
# Legend and colorbars share the right dock without replacing one another.
cable_preview = preview(
    mv_design;
    backend = :cairo,
    display_plot = false,
    controls = false,
    display_id = true,
    size = (950, 700),
    # Group the physical core-wire regions into one presentation entry.
    # Geometry tags, terminals, and electrical construction remain untouched.
    legend_group = region -> region.source.tag === :wire ?
                             :stranded_core : region.source.tag,
    legend_labels = Dict(:stranded_core => "Stranded core"),
    legend_position = :right,
    legend_attributes = (; nbanks = 2),
    colorbar_position = :right,
    colorbar_attributes = (; vertical = false)
)
cable_preview.figure #hide
````

`legend_group` accepts either the shown callback over a `PlacedRegion` or a
tag-to-group dictionary. `legend_labels` maps the resulting presentation group
to text. This keeps display grouping independent of physical tags. For detailed
annotation, mutate `cable_preview.axes[1]` and add ordinary Makie plots. Each
object in `cable_preview.colorbars` is a native `Colorbar`.

### Cable-design collection

A collection preview repeats the same detached geometry path for every design
and assigns the caller-requested layout. The material ranges are aggregated
once, so every panel uses comparable colors and one shared set of scales.
Insulating regions carry sparse diagonal marks over their existing material
colors. Pass `display_dielectric_pattern=false` to omit these marks in a
design, collection, or system preview. Semicon and conductor fills retain
their own material colors.

````@example plotting
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
````

Use any sufficient `(rows, columns)` layout. The default omits layer legends;
request selected local legends with `panel_legends`, using the same logical
grid-position rules as line dashboards.

### Cable-system cross-section

A system preview resolves every placed region into the system frame, adds
reference geometry such as the earth interface, and derives limits from the
physical placement or `zoom_factor`. Earth properties use the same atomic
color-scheme rules as cable materials, with a separate logarithmic
resistivity palette: slate at 0.1 Ω·m, taupe at 100 Ω·m, and ochre at 10⁴ Ω·m.
Horizontal earth fills follow pan, zoom and figure resizing while interfaces
stay at their physical depths. A semi-infinite basement covers the remainder
of the view; a finite final layer retains its declared bottom.

````@example plotting
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
````

Modify physical contents before previewing; modify visual annotations after
previewing. `zoom_factor` changes only the initial view, and the reset control
returns to those computed limits.
The light blue sky fades from the upper axis limit to transparency at `z=0`,
stretching with the view and disappearing in entirely underground views.
It is a surface cue with no material-property meaning; use
`display_surface_gradient=false` to disable it independently of earth colors.
Vertical strata are currently not rendered in the system preview.

## Material colors and native colorbars

### One reusable material scheme

A color scheme is the atom. `material_property_ranges` obtains values from a
design, a design collection, or the material defaults. `materialcolors`
turns exactly one property and range into a Makie-compatible named tuple:
`label`, `colormap`, `limits`, and `ticks`. `materialscale!` places exactly
that one scheme wherever the caller chooses.

````@example plotting
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
````

The scheme contains no placement. Put that `Colorbar` in any `GridPosition`,
combine it with a heatmap, or reuse the same scheme in another figure. Define a
new property by defining another palette producer next to `materialcolors`;
physical range collection remains a DataModel concern.
Relative permeability uses the existing logarithmic range from 1 to 300.
Unity leaves the base color unchanged; indigo tint becomes visible before
the progression toward magenta at high permeability. The permeability scale
shows that same tint over a neutral reference color.
Dielectric marks are native Makie pattern tiles. Cairo embeds the small
bitmap tiles in SVG/PDF output while preserving the surrounding geometry.

### High-level material-scale reference

`show_material_scale` is a preview option that combines the three
default property schemes. Its result still exposes three independent native
colorbars; it is not the reusable unit of the color API.

````@example plotting
material_scale = show_material_scale(
    backend = :cairo,
    display_plot = false,
    controls = false,
    size = (850, 360),
    figure_title = "Reusable material-property schemes",
    colorbar_attributes = (; vertical = false)
)
material_scale.figure #hide
````

For one property, use the preceding `materialscale!(position, scheme)` pattern.
For a different high-level combination, compose schemes in the caller's own
Makie layout.

## Monte Carlo result recipes

Monte Carlo methods first publish the requested marginal from retained
samples and/or `HistogramDensity`. Cable-constant requests accept `R`, `L`,
`C`, `G`, or `(selector, assembly)`. Matrix-valued line results require an
exact request such as `@observe R[1, 1, 3]`. Native Makie function identity
chooses the visual primitive.

### Sample histogram

`Makie.hist` uses retained samples and calls native `hist!`. `bins` and
`normalization` retain their Makie meanings; unit options are consumed while
publishing the marginal.

````@example plotting
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
````

Pass ordinary `hist!` attributes such as `color`, `strokewidth`, or
`transparency` through the remaining keywords. Use `(R, assembly)` when a
`CableConstants` result contains more than one assembly.

### Retained-model probability density

`Makie.stairs` consumes the retained histogram model and calls native
`stairs!` with post-step edges. It does not regenerate a stochastic result.

````@example plotting
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
````

With `bins=n`, a different binning is derived from retained samples; the
stored model is unchanged and no stochastic calculation is repeated. If
samples were not retained, only the stored model's bin count is available.
With `bins=nothing` (the default), the retained model is reused, or a model
is derived with automatic binning if only samples were retained. Constant
samples always produce one finite-width bin. Use native `stairs!` keywords
for appearance.

### Empirical cumulative distribution

`Makie.ecdfplot` consumes retained samples and delegates the empirical curve
to Makie's `ecdfplot!` recipe.

````@example plotting
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
````

The returned axis can be combined with confidence bands or additional native
curves; the addon only owns the initial empirical series and its legend group.

### Retained-model cumulative distribution

Raw `Makie.lines` first acquires the UQ-owned CDF product. Its observed method
draws retained coordinates with `lines!`; the renderer does not estimate a CDF.

````@example plotting
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
````

Increase visual resolution by extending the owner-side model grid settings; add
purely visual reference curves directly to `model_cdf.axes[1]`.

### Sample/model Q-Q plot

`Makie.qqplot` requests both retained products, lets the UQ owner calculate
matching quantile pairs, draws native scatter points, and optionally adds the
identity reference line.

````@example plotting
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
    legend_attributes = (; halign=:left,valign=:top,backgroundcolor=(:white,0.92)),
    legend_overflow = :show_all,
    color = :purple,
    markersize = 10
)
quantile_plot.figure #hide
````

Set `qqline=:none` to remove the reference. Other scatter attributes are
forwarded to `scatter!`; data and units remain observation concerns.

## Callable controls and retained assembly points

A widget callable receives the final live handle once per figure. Its native
callback uses the same status and frame-preserving shell operations as built-ins.
`controls=false` skips standard and custom widgets. Closing a display window
leaves the retained handle editable, exportable, and redisplayable.

````@example plotting
function add_reset_y!(p)
    addwidget!((plot,cell) -> Button(cell;label="Reset Y"),p,:reset_y;
        event=button -> button.clicks,
        callback=(plot,_) -> resetview!(plot;x=false,y=true),
        success="Y view reset")
    nothing
end
widget_plot=Makie.plot(parameters; ydata=((R,1,1,:),),layout=(1,1),
    backend=:cairo,display_plot=false,widgets=(add_reset_y!,));
widget_plot.figure #hide
````

CableConstants supply categorical assembly coordinates. The first-seen union
below is core, sheath, screen; neither observation invents the other's point.
Their intervals retain their original uncertainty meaning. Categorical points
are primary data, so decorative-marker sampling does not remove them.

````@example plotting
assembly_a=CableConstants([:core,:sheath],[1.,2.].*1e-4,[2.,3.].*1e-7,
    [3.,4.].*1e-10,[1.,2.].*1e-9,50.);
assembly_b=CableConstants([:core,:screen],[1.1,2.2].*1e-4,[2.1,3.1].*1e-7,
    [3.1,4.1].*1e-10,[1.1,2.1].*1e-9,50.);
assembly_plot=Makie.plot((first=assembly_a,second=assembly_b);ydata=(R,),
    backend=:cairo,display_plot=false,controls=false);
assembly_plot.figure #hide
````

## Composing figures with `plotwindow`

### A caller-owned 2×2 plot

`plotwindow` is the escape hatch when no high-level recipe is appropriate. It
creates the figure and controls, then passes its content `GridLayout` to the
callback. The callback uses normal Makie constructors; afterward `plotwindow`
discovers the native axes and attaches shared numeric scale, reset and export
controls. Explicit `axis=(...)` overrides apply to callback-created axes;
otherwise their native construction settings are retained.

````@example plotting
custom_dashboard = LineCableModels.plotwindow(
    title = "Caller-owned diagnostics",
    figure_title = "Four caller-owned Makie axes",
    size = (900, 650),
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
````

Nothing prevents nested layouts, `Axis3`, `Colorbar`, `Slider`, custom Makie
recipes, or arbitrary plot primitives in the callback. If a native composition
needs no LineCableModels controls or export settings, use `Figure` directly.

## What the addons do

### Scientific axes, scale controls, and limits

Scientific recipes supply quantity and unit labels. PlotBuilder formats
every linear numeric axis, including previews, statistical plots, and `plotwindow`
axes: the displayed limits determine one engineering power-of-ten multiplier
(powers of three) in the axis label. Ticks show plain decimal mantissas, never
another exponent. Tick density follows each data rectangle and the rendered
label size. For example, a linear frequency view can show `0, 2, 4, 6, 8, 10`
with `Frequency [Hz] ×10⁶`; the retained frequencies are unchanged.
Zooming, panning, changing limits, and resetting keep the ticks and multiplier
synchronized. Native custom tick formatters or explicit tick labels override
automatic formatting; setting the formatter back to `Makie.automatic` restores
it. Native tick positions persist across scale changes. Caller tick functions
and custom locator objects retain native formatting; resetting the tick
attribute to `Makie.automatic` restores the shared locator.
Date/category axes and other native transforms retain Makie's own presentation.
These defaults use the native tick locators in Makie 0.24.11 or newer.
Logarithmic views spanning less than two decades show decimal values at
logarithmic positions, with one engineering multiplier when needed. Broader
views show integer powers of ten without another multiplier. Both x and y
use the same policy, fitting actual label spacing after the coordinate transform.
The x/y toggles validate current visible data and uncertainty bounds before
changing the page. Native numeric `plotwindow` axes share these controls.
Automatic near-constant positive log ranges use modest multiplicative padding
(`c/1.05` to `c*1.05`, enlarged for uncertainty), not whole-decade bounds.

Native Axis keywords (`xticks`, `limits`, `ytickformat`, etc.) and native
series attributes (`linewidth`, `color`, etc.) can be passed at construction.
Explicit `axis=(...)`, `figure=(...)` and per-series `series_attributes`
override shared defaults. Subsequent native mutations remain authoritative.
For names shared by Axis and a plot, the unqualified form targets Axis.
An explicit native `figure.size` overrides `fig_size`; portrait sizes stay portrait.

UQ benchmark plots with ordinary `ydata=(R,L,G,C)` overlay mean ±1 standard
deviation from retained owner moments. No uncertain result is constructed to draw error bars. Explicit `(statistics,L,std)`
requests retain statistic-only plots; request these two products in separate
calls. `uncertain(result, configuration)` returns the stored uncertainty-bearing
core without reconstruction. MC constructs this marginal representation during
aggregation, without inferring joint correlations.
Every eligible numeric route retains log controls for zero or negative support;
adaptive logarithmic panels use a sign-preserving pseudo-log transform with
`log1p`/`expm1` evaluation to retain tiny signed values near zero.

Limits are calculated from finite visible data, including measurement error
bounds. Constant and near-constant series receive at least ±5% padding around
a nonzero baseline, enlarged for visible uncertainty. An exactly zero series
without uncertainty uses a neutral nonzero range. Near-constant means that the
endpoints agree within `sqrt(eps(Float64))` relatively in view coordinates.
This is only an automatic-view rule: no sample is rounded. Small physical
values with meaningful relative variation still receive a tightly fitted view.
Legend visibility changes
trigger another limit pass, so hiding a dominant curve exposes the remaining
data instead of leaving a stale range.
Explicit native limits, including one-sided limits, remain authoritative.
Reset refits automatic bounds and restores explicit bounds. Changing x/y scale
refits that automatic dimension while preserving the other dimension's current
view. Incompatible log limits are rejected before the page changes. Use native
`autolimits!(axis)` explicitly when manual bounds should be discarded.

### Legends and docks

`legend_position` selects only placement. `:inside` overlays the union of the
figure's axis viewports or the single axis viewport for a panel legend;
native `halign`/`valign` select corner, center, or fractional alignment. Side-grid slots
remain outside the plot area. `legend_attributes` is merged into the native
`Legend` constructor, so orientation, bank count, padding, background,
alignment, and other Makie options remain available. `legend_overflow` is the
one addon settings: `:ellipsis` fits entries to the current bounding box and
restores them when space returns; `:show_all` always retains all entries.
`figurelegend!` and `panellegend!` move retained native objects. Label changes
retain their source bindings, so placement, title, and semantic labels can be
changed after construction. Clicking/toggling a grouped Makie legend entry
continues to affect every plot handle in that group.

High-level recipes may put a legend and colorbars in the same dock. PlotBuilder
creates a nested `GridLayout` and gives each block its own cell. Resizing moves
the entire canvas/dock group; it does not anchor a block to the raw window
while stretching only the outer figure.

### Colorbars

`colorbar_position` and `colorbar_attributes` mirror the legend placement
rules. A preview can request several schemes, but
`_addon_colorbar!` always consumes one scheme and creates one native
`Colorbar`. The reusable public atom remains
`materialcolors(property, range)` plus `materialscale!(position, scheme)`.

### Observables and ownership

Two distinct concepts are intentionally present. `@observe` belongs to the
scientific observation grammar: it selects a physically meaningful view of a
result vault, including matrix and frequency coordinates. Makie's
`Observable` type drives live UI state such as controls, scales, limits,
visibility, and layout bounds. The adapter consumes the former and wires the
latter; it does not replace either one with an `AxisBehavior` aggregate.
Keeping the `UIPlot` alive keeps the native figure and its subscriptions alive.

Native mutation is therefore the normal extension mechanism:

````@example plotting
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
````

### Responsive layout

Automatic diagonal and preview flow pages may reflow locally when resized.
Their page membership and identities stay fixed. Explicit layouts and matrix
topology stay fixed. Physical aspect belongs to each panel, so circular designs,
wide systems, and 1×4 preview collections remain correctly scaled. A later guide,
title, or widget change refits only that window around its current data frames.

### Current-state SVG export

Load CairoMakie before creating a plot to include its Save button. For GL
interactivity and SVG export, import both CairoMakie and GLMakie, then select
`backend=:gl`. Installed but unloaded CairoMakie does not enable the button.
Loading CairoMakie later enables direct export of an existing plot; recreate
the plot to add its button. Export never loads packages. The toolbar reports
file errors in the status row; direct `export_svg` calls throw them to the caller.

[`export_svg`](@ref) saves the current live figure through CairoMakie. For a
publication export it temporarily hides the toolbar and status row, switches
the figure's font roles to Makie's LaTeX font theme, uses a white background,
and then restores every changed observable. Saving does not activate a backend.
Caller-added plots, visibility, scales, limits, and annotations are saved
because export does not reconstruct an earlier specification.
Interactive zoom and pan are retained in both the SVG and the live window.

````@example plotting
export_directory = mktempdir();
export_svg(
    owned_plot;
    path = joinpath(export_directory, "caller_owned_resistance.svg"),
    theme = :publication,
    open_file = false
); #hide
nothing #hide
````

## Adding or changing a managed recipe

A new high-level recipe should preserve the same ownership boundary:

1. Expose numerical observations or physical geometry through the owning
   module's existing public protocol. Do not add plot preparation, labels,
   colors, layout, or Makie types to scientific owners.
2. Add request normalization and the narrow public dispatch method in
   `LineCableModelsMakieExt`, then call
   native Makie constructors or primitives there.
3. Reuse only the addon services the recipe needs: `_addon_shell`,
   `_addon_axis!`, `_addon_finish!`, or `plotwindow`. Do not create a second
   plot specification or an optional adapter hierarchy.
4. Return `UIPlot` with the actual native objects and leave further mutation to
   the caller.
5. Add the real call to this literate gallery, a Cairo rendering assertion,
   and an interactive GL inspection fixture when resizing or widgets matter.

A purely visual variation usually needs no new managed recipe: pass a native
attribute, mutate the returned block, add a Makie primitive, or start from
`plotwindow`. A new recipe is justified when LineCableModels owns meaningful
retained coordinate presentation, physical geometry, or a reusable piece
of scientific interaction.

## Implementation

Scientific objects supply observations, geometry, material properties, and
units. `LineCableModelsMakieExt` converts these values into Makie plots.
Matrix coordinates identify subplots; result containers identify overlaid
series. PlotBuilder adds labels, controls, and SVG export.

| Implementation | Responsibility |
|:--|:--|
| `src/plotbuilder/` | Optional entry points and `UIPlot` |
| `ext/LineCableModelsMakieExt/recipes/line_data.jl` and `comparison_data.jl` | Select observations for plotting |
| `src/datamodel/preview/geometry.jl` and `materials.jl` | Geometry and material ranges |
| `ext/LineCableModelsMakieExt/recipes/preview_data.jl` | Prepare geometry for drawing |
| `ext/LineCableModelsMakieExt/material_colors.jl` | Material palettes |
| `ext/LineCableModelsMakieExt/shell.jl` | Figure layout, axes, limits, and controls |
| `ext/LineCableModelsMakieExt/recipes/*_render.jl` | Draw declaration geometry |
| `ext/LineCableModelsMakieExt/layout.jl` and `guides.jl` | Capacity, native frames, and guide composition |
| `ext/LineCableModelsMakieExt/controls.jl` | Native widget ownership and common actions |
| `ext/LineCableModelsMakieExt/export_presentation.jl` | Temporary native export presentation and restoration |
| `ext/LineCableModelsMakieExt/montecarlo.jl` | Statistical plots |
| `ext/LineCableModelsMakieExt/native_export.jl` | SVG export |
