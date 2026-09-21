"""
    plot(observed::ObservedResult, selection=nothing; ydata=nothing, kwargs...)
    plot(observed::AbstractVector{<:ObservedResult}, selection=nothing; ydata=nothing, kwargs...)
    plot(completed_result, selection=nothing; ydata=nothing, kwargs...)

Present retained scientific quantities with a loaded Makie backend. Completed
numerical results are conveniences: they construct `ObservedResult` objects and
call the same public observed-input method. `PlotBuilder.plot` and
`LineCableModels.plot` are the same function.

`selection` and `ydata` are alternative spellings of the same request; supplying
both is an error. Requests retain `@observe` syntax and original coordinates.
An existing observation supplies values, units, descriptions, uncertainty, and
scientific groups. Rendering does not reacquire a result or run a comparison.

# Selection and acquisition

Raw conveniences accept line, series, shunt, cable-constant, parametric, and UQ
results, ordinary supported tuples/vectors, named collections, and report
artifacts. Standalone series/shunt inputs also accept a frequency vector before
the selection. Raw references become separate atomic observations. Report
artifacts forward their observed candidates and their observed reference.

Raw-only acquisition keywords are `clip`, `atol`, and `frequencies`. Single
primary requests use the observation owner's pair completion; statistical
products do not trigger it. Observed inputs reject new clipping decisions,
thresholds, or replacement sample coordinates.

`units`, `length_unit`, `quantity_units`, and `frequency_unit` also apply to
retained inputs, including reports. The existing `ObservedResult(existing; ...)`
operation re-expresses compatible units once before drawing; omitted options
preserve recorded units, masks, errors, timings, and uncertainty dependencies.
`freq_unit` is a spelling of `frequency_unit`; supplying both is an error.
For example, `plot(report; ydata=(R,), length_unit=:base)` displays its retained
candidate and reference curves per meter without rebuilding the report.

`problem` and `formulations` select original recorded identities. `band` selects
saved comparison samples through the observation owner, retaining each trace's
own coordinates and reference association. It does not calculate new errors.

# Layout and native presentation

- `layout=nothing` resolves the nominal panel capacity from selected products:
  the largest declared matrix extent, or a near-square flow arrangement.
  An explicit `(rows, columns)` supplies a positive nominal capacity.
- Each quantity/statistical meaning has separate figure families. A full 3×3
  matrix at `layout=(2,2)` has four pages with extents `(2,2)`, `(2,1)`, `(1,2)`,
  `(1,1)` per quantity. `layout=(1,1)` produces nine pages per quantity.
  Explicit diagonal products paginate compactly with original `(i,i)` identities.
- `fig_size` describes the complete nominal capacity. `figure=(size=...,)` wins.
  Residual pages fit their occupied tracks with the same initial data frames;
  there is no forced landscape orientation or permanent sibling synchronization.
- Assembly products use categorical points and uncertainty intervals, with a
  first-seen union of recorded assembly names. Scalars/vectors without physical
  coordinates use honest element-index views. Full modal matrices retain every
  selected coefficient, including zero and residual off-diagonals.
- `series_labels`, `reference`, and `series_attributes` control trace identity and
  native appearance. Attributes accept one NamedTuple or an aligned tuple/vector.
  Candidate slots are assigned before filtering; a separate reference does not
  shift them. References default to black solid curves and hollow circles.
- `errorbar_sampling` defaults to `:staggered` for multiple displayed series and
  `:all` for one. Full uncertainty support still controls limits. Explicit curve
  markers use every original sample; categorical/scalar points are never thinned.
- `xscale`, `yscale`, `xlabel`, `ylabel`, native limits/ticks/formatters, `axis=(;)`,
  and `figure=(;)` configure native objects. Constructor groups beat shared
  attributes; per-series overrides beat shared series settings. Later native
  edits retain authority. Unknown attributes fail with a diagnostic.
- Frequency X defaults to adaptive `:log10`; other numeric dimensions default
  to linear. Adaptive log uses stable signed log when visible support or bounds
  contain zero/negative values. The native `log10` function remains strict.
  Categorical axes retain native conversion and have no X-log toggle.
- Titles use `title`, `title_prefix`, `figure_title`, `title_attributes`, and
  `panel_titles`. Positional panel titles bind before pagination; dictionary and
  function selections retain original panel identities.
- Legends use `legend_position`, `legend_title`, `legend_attributes`,
  `legend_overflow`, and `panel_legends`. Multiple/explicitly labelled result
  series default to a bottom legend with `:show_all`; one unlabelled series has
  none. Measured wrapping preserves labels and native interaction targets.
- `colorbar_position`, `colorbar_group_attributes`, `colorbar_attributes`, and
  `guide_gap=8` place existing scale content; they do not fabricate scales.
  Native `halign`/`valign` supply symbolic or fractional alignment.
- `backend`, `display_plot=true`, `controls=true`, `widgets=()`,
  `export_theme=:default`, and `open_export=true` control display and UI behavior.
  Widget callables receive the final live handle once per figure. Hiding controls
  does not disable [`axisscale!`](@ref) or [`resetview!`](@ref).

Numerical axes share size-aware ticks, engineering multipliers, and relative
near-constant padding. Scale changes preflight the complete page and preserve
orthogonal views and configured bounds. Observation owns clipping and uncertainty
meaning; the shell never clips, rounds, or recalculates scientific eligibility.
Guide/title/widget changes refit the affected outer window around its current
frames. Closing a display window leaves its retained handle reusable.

# Returns

One [`UIPlot`](@ref), or an ordinary `Vector{UIPlot}` for multiple figures. Its
figure, axes, guides, controls, and status remain live native objects.

# Examples

```julia
using LineCableModels
using LineCableModels: plot
using CairoMakie
frequency = [1.0, 10.0, 100.0]
impedance = reshape(complex.([1.0, 2.0, 3.0], [2.0, 3.0, 4.0]), 1, 1, :)
raw = LineParameters(impedance, impedance .* 1e-6, frequency)
r_request = @observe R[:, :, :]
pair = (@observe(R[:, :, :]), @observe(L[:, :, :]))
a = plot(raw; ydata=(r_request,), clip=true, length_unit=:kilo)
o = ObservedResult(raw, (r_request,); complete_pairs=true,
    clip=true, length_unit=:kilo)
b = plot(o; ydata=(r_request,), length_unit=:base)
c = plot(raw; ydata=pair, layout=(1,1))
d = plot(ObservedResult(raw, pair); ydata=pair, layout=(1,1))
```

`c` and `d` each contain separate R and L figures.
"""
function plot end

function plot(args...; kwargs...)
    throw(ArgumentError(
        "Plotting is optional. Load CairoMakie, GLMakie, or WGLMakie before calling plot.",
    ))
end

"""
    preview(source; kwargs...)

Preview a cable design, a collection of cable designs, or a cable system with
a loaded Makie backend.

# Keywords

- `display_dielectric_pattern=true`: Fill insulating regions with sparse
  diagonal marks over their material color. Applies to all three preview routes.
- `earth_model=nothing`: Static earth model for a system preview. Horizontal
  strata retain their physical depths while their visible coverage follows the
  axis view. Vertical strata are not rendered.
- `display_surface_gradient=true`: For a system with horizontal earth, add a
  light blue sky strongest at the upper axis limit, fading toward the white
  or transparent background at the surface `z=0` \\[m\\]. The fade stretches
  with the view and is hidden when the view lies entirely underground. This
  decoration does not encode a material property.
- `zoom_factor=nothing`: Initial system-view span multiplier. The reset
  control restores that initial view.

# Returns

- One [`UIPlot`](@ref), or an ordinary vector when a collection exceeds `layout`
  capacity. Collection panels retain original integer indices; material ranges
  are shared across pages. `size` supplies the nominal figure dimensions.

# Notes

Material colors retain nominal physical properties. Magnetic tint progresses
from indigo to magenta on a logarithmic relative-permeability range; earth uses
its own logarithmic resistivity palette. Dielectric marks use native Makie
pattern tiles, including in SVG/PDF exports.
"""
function preview end

function preview(args...; kwargs...)
    throw(ArgumentError(
        "Plotting is optional. Load CairoMakie, GLMakie, or WGLMakie before calling preview.",
    ))
end

"""
    show_material_scale(; kwargs...)

Display the three independently defined material color schemes as a compact
reference figure. Use [`materialscale!`](@ref) to place any one scheme in a
caller-owned Makie layout.
"""
function show_material_scale end

function show_material_scale(args...; kwargs...)
    throw(ArgumentError(
        "Plotting is optional. Load CairoMakie, GLMakie, or WGLMakie before calling show_material_scale.",
    ))
end

"""
    export_svg(plot::UIPlot; path=nothing, theme=nothing, open_file=nothing)

Save the current live Makie figure in `plot` as SVG through CairoMakie and
return the absolute output path. `theme` may be `:default` or `:publication`.
The SVG retains the current zoom and pan without resetting the interactive view.

Load CairoMakie explicitly before exporting. For an interactive GLMakie window
with SVG export, import both backends and select `backend=:gl` when plotting.
The SVG button is created only when CairoMakie is already loaded. Loading it
later enables this function on existing plots; recreate a plot to add its button.
Export does not activate a backend and restores the live figure state. The toolbar
reports file errors in the window's status row; direct calls throw the
corresponding exception. An unloaded CairoMakie renderer raises `ArgumentError`
before filesystem or figure changes.
"""
function export_svg end

"""
    figurelegend!(plot::UIPlot; position, title, overflow, legend_labels, kwargs...)

Update the figure legend from the shell's native series groups. Omitted options
preserve current state; `position=nothing` detaches and hides the guide. Native
`halign`, `valign`, margins, and style attributes remain editable. Removal and
restoration retain native styles and series visibility.
"""
function figurelegend! end

"""
    panellegend!(plot::UIPlot, panel; kwargs...)

Create or replace a native Makie legend scoped to one logical plot panel.
`panel` may be the stable panel identity returned by a recipe or its compatible
grid position.
"""
function panellegend! end

"""
    figuretitle!(plot::UIPlot, title; kwargs...)

Create, replace, or remove the figure-wide native Makie title. Pass `nothing`
to remove it.
"""
function figuretitle! end

"""
    paneltitle!(plot::UIPlot, panel, title)

Set the native axis title for one logical plot panel. Pass `nothing` to clear
it.
"""
function paneltitle! end

"""
    plotwindow(callback; title, figure_title=nothing, size=(800, 400), kwargs...)

Build the standard Makie shell, pass its caller-owned content `GridLayout` to
`callback`, and return a [`UIPlot`](@ref). The callback uses ordinary Makie and
is not constrained by a renderer-independent plot specification.
Numeric axes share scale, reset and export controls. Native `axis=(...)`
attributes explicitly override callback-created axes; otherwise their native
construction settings are retained. Shared series and `figure=(...)` attributes
follow the same rules as [`plot`](@ref).
"""
function plotwindow end

"""
    materialcolors(property, [range]; alpha=1.0)

Construct one reusable material color scheme. Palette selection is separate
from [`materialscale!`](@ref), which only renders a supplied scheme.
"""
function materialcolors end

"""
    materialscale!(position, scheme; kwargs...)

Place one native Makie color scale at `position`. `scheme` supplies one label,
colormap, limits, and tick definition. The primitive never chooses or combines
material properties.
"""
function materialscale! end

"""
    axisscale!(p::UIPlot, dimension::Symbol, scale; panel=nothing)

Set one displayed coordinate scale on all eligible numeric axes, or on the
original identity selected by `panel`. `:linear` selects the identity scale;
`:log10` adapts to signed support, while the native `log10` function requires
strictly positive support. `:pseudolog10` explicitly selects stable signed log.
Preflight covers the complete selection before mutation and preserves the
orthogonal view and configured limits. Return `p`.
"""
function axisscale! end

"""
    resetview!(p::UIPlot; panel=nothing, x::Bool=true, y::Bool=true)

Refit selected automatic view dimensions to current visible support, including
full enabled uncertainty intervals. Preserve configured full or partial limits.
`panel=nothing` selects every panel; another value selects its original identity.
The Boolean keywords `x=true` and `y=true` select dimensions. Return `p`.
"""
function resetview! end

"""
    addwidget!(builder, p::UIPlot, key::Symbol; event=nothing, callback=nothing, success=nothing)

Construct one custom native control with `builder(p, slot)` and register it
under a unique nonstandard Symbol `key`. Supply both `event` and `callback`, or
neither. `event(control)` returns the event observable; actual notifications call
`callback(p, value)` once and may set `success` on `p.status`. Return the native
control. Figures constructed with `controls=false` reject widget additions.
"""
function addwidget! end

"""
    removewidget!(p::UIPlot, key::Symbol)

Remove a custom widget, its native subtree, and its owned event subscriptions.
Release and repack its toolbar slot without changing other controls. Standard
and unknown keys are rejected. Return `p`.
"""
function removewidget! end

"""
    figurecolorbars!(p::UIPlot; position, group_attributes, native_bar_attributes...)

Update the placement and native attributes of the figure's retained color scales.
Omitted arguments preserve current state; `position=nothing` removes displayed
scales. `group_attributes` controls group presentation. Return the current native
colorbar collection. This operation does not fabricate scales for line results.
"""
function figurecolorbars! end
