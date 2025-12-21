

# -----------------------------------------------------------------------------
# Spec-level traits (configuration surface)
# -----------------------------------------------------------------------------

"""
Kind of plotting primitive this spec corresponds to.

Typical values: :line, :heatmap, :hist, :bar, :surface, ...
"""
plot_kind(::Type{S}) where {S <: AbstractPlotSpec} = :unknown

"""
Axes that can be toggled to log-scale at the UI level.

Returns a tuple of axis dims, e.g.:

	enable_logscale(::Type{MySpec}) = (:x,)       # only x can log
	enable_logscale(::Type{OtherSpec}) = (:x,:y)  # x and y
"""
enable_logscale(::Type{S}) where {S <: AbstractPlotSpec} = ()

"""
Domain/container type this spec expects to dispatch on.

Example:

	dispatch_on(::Type{MyRPlotSpec}) = LineParameters
"""
dispatch_on(::Type{S}) where {S <: AbstractPlotSpec} = Any

"""
Default figure size for this spec type, in pixels.

This removes hard-coded consts in PlotUIComponents and pushes that
decision into the grammar layer.
"""
default_figsize(::Type{S}) where {S <: AbstractPlotSpec} = (800, 400)

"""
	axis_quantity(::Type{S}, ::Val{dim}) where {S<:AbstractPlotSpec, dim}

Return the default semantic quantity for axis `dim` in spec `S`, when it
does not depend on which data source is selected.

	axis_quantity(::Type{S}, ::Val{dim}, ::Val{datakey})

Higher-ranked variant: given a data selector `datakey` (e.g. :f, :R, :L, :Z),
return the semantic quantity for axis `dim`.

The `datakey` is a symbol describing *where* data comes from in the container;
it is not necessarily equal to the quantity name used in `QuantityTag{Q}`.
"""
axis_quantity(::Type{S}, dim::Symbol) where {S <: AbstractPlotSpec} =
	QuantityTag{:unknown}()

axis_quantity(::Type{S}, ::Val{dim}) where {S <: AbstractPlotSpec, dim} =
	axis_quantity(S, dim)

axis_quantity(
	::Type{S},
	::Val{dim},
	::Val{datakey},
) where {S <: AbstractPlotSpec, dim, datakey} =
	axis_quantity(S, dim)

"""
	index_keys(::Type{S}) where {S<:AbstractPlotSpec}

Semantic index parameters this spec uses to address elements of its underlying
tensors (e.g. (:i, :j) for matrix-like data, (:i, :j, :k) for 3D, etc.).

By default, no index keys are assumed. Specs that work on per-element data
over frequencies should typically override this to `(:i, :j)`.
"""
index_keys(::Type{S}) where {S <: AbstractPlotSpec} = ()

"""
	ranged_keys(::Type{S}) where {S<:AbstractPlotSpec}

Index keys among `index_keys(S)` that may also be specified as ranges.

A ranged index may be:

  • `Int`            → single position (e.g. `k = 5`)
  • `AbstractUnitRange{<:Int}`  → slice (e.g. `k = 1:10`, `k = 3:2:99`)
  • `:`              → full range

Default: empty tuple (no ranged indices).
"""
ranged_keys(::Type{S}) where {S <: AbstractPlotSpec} = ()

"""
	geom_axes(::Type{S}) where {S<:AbstractPlotSpec}

Geometric axes used by this spec, in order.
Default is 2D (:x, :y). If your spec is 3D, override to return (:x, :y, :z).
"""
geom_axes(::Type{S}) where {S <: AbstractPlotSpec} = (:x, :y)


# --------------------------------------------------------------------------
# Title / legend grammar traits
# --------------------------------------------------------------------------

"""
	default_title(::Type{S}, nt) where {S<:AbstractPlotSpec}

Return the default plot title for this spec, given the resolved input `nt`.
`nt` is the output of `resolve_input(S, ...)`, so its structure is spec-defined.
"""
default_title(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec} = ""

"""
	legend_labels(::Type{S}, nt) where {S<:AbstractPlotSpec}

Return the legend entry labels for this spec, given the resolved input `nt`.
Length of the returned vector must match the number of primitives produced
by `build_series(S, nt)`.
"""
legend_labels(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec} = String[]


# -----------------------------------------------------------------------------
# Quantity-level unit and label hooks (using UnitHandler)
# -----------------------------------------------------------------------------

"""
Spec-level hook: unit for specific spec + axis dim.

By default, delegates to `display_unit(quantity)`, since plotting is a
display concern. Specs can override for special cases if needed or to
honour user overrides.
"""
axis_unit(::Type{S}, q::QuantityTag, dim::Symbol) where {S <: AbstractPlotSpec} =
	display_unit(q)


# --------------------------------------------------------------------------
# Input / backend grammar traits
# --------------------------------------------------------------------------

"""
	input_kwargs(::Type{S}) where {S<:AbstractPlotSpec}

Plot-level *semantic* kwargs understood by this spec.

These describe what is plotted or how the data is selected/sliced
(e.g. :quantity, :stat, :i, :j, :k, :values_expr, ...).
"""
input_kwargs(::Type{S}) where {S <: AbstractPlotSpec} = ()

"""
	renderer_kwargs(::Type{S}) where {S<:AbstractPlotSpec}

Figure kwargs that are simply forwarded to the renderer that will be processed by the backend (Makie today, whatever tomorrow).
"""
renderer_kwargs(::Type{S}) where {S <: AbstractPlotSpec} = ()

"""
Root container inside `obj` for the data of axis `dim`.

Default: `nothing` → use `obj` itself.

For example, if `obj.stats` is a NamedTuple of tensors and y-axis data
comes from there, define:

	data_container(::Type{MySpec}, ::Val{:y}) = :stats
"""
data_container(::Type{S}, ::Val{dim}) where {S <: AbstractPlotSpec, dim} = nothing

# Child key for axis `dim` inside the container returned by `data_container`.
# child_key = getfield(nt, select_field(S, Val(dim)))  # e.g. :mean
# data      = obj.container.datakey[i,j,k].(child_key)
select_field(::Type{S}, ::Val{dim}) where {S <: AbstractPlotSpec, dim} = nothing

"""
	input_defaults(::Type{S}, obj) where {S<:AbstractPlotSpec}

Defaults for semantic kwargs declared in `input_kwargs(S)`.

May depend on the dispatched object `obj` (e.g. pick default `:quantity`
from `obj` contents).
"""
input_defaults(::Type{S}, obj) where {S <: AbstractPlotSpec} = NamedTuple()

"""
	renderer_defaults(::Type{S}, obj) where {S<:AbstractPlotSpec}

Defaults for figure kwargs declared in `renderer_kwargs(S)`.

This is where a spec declares its default color/linestyle/whatever,
possibly depending on `obj` (e.g. per-phase colors).
"""
renderer_defaults(::Type{S}, obj) where {S <: AbstractPlotSpec} = NamedTuple()

"""
How to group dataseries into figures.
Options: :auto -> let the machinery decide;
		 :single -> one dataseries in one plot area (panel), same axis;
		 :overlay_ij -> one plot area (panel), overlay all (i,j) on the same axis - target/leaf resolved to one field;
		 :overlay_fields -> one plot area (panel), overlay all fields on the same axis - data container resolved to one pair (i,j).
		 :per_ij_overlay_fields -> multiple plot areas (panels), one per (i,j), overlay all fields.
"""
grouping_mode(::Type{S}) where {S <: AbstractPlotSpec} = :auto

"""
How to render figures into Makie windows.
Options: :windows -> one panel per window;
		 :grid -> all panels in a single window, arranged in a grid;
		 :tabs -> TBD: all panels in a single window, arranged in tabs.
"""
figure_layout(::Type{S}) where {S <: AbstractPlotSpec} = :windows  # default

# --------------------------------------------------------------------------
# Complex quantity / "as" traits (default: disabled)
# --------------------------------------------------------------------------

has_complex_qty(
	::Type{S},
	::Val{dim},
	::Val{datakey},
) where {S <: AbstractPlotSpec, dim, datakey} =
	false

complex_as(
	::Type{S},
	::Val{dim},
	::Val{datakey},
) where {S <: AbstractPlotSpec, dim, datakey} =
	(:re, :im, :abs, :angle)

complex_as_default(
	::Type{S},
	::Val{dim},
	::Val{datakey},
) where {S <: AbstractPlotSpec, dim, datakey} =
	:re

# View-aware axis_quantity: fallback keeps existing grammar intact
axis_quantity(
	::Type{S},
	::Val{dim},
	::Val{datakey},
	::Val{as},
) where {S <: AbstractPlotSpec, dim, datakey, as} =
	axis_quantity(S, Val(dim), Val(datakey))

# Z: re/im correspond to R/X
axis_quantity(
	::Type{S},
	::Val{dim},
	::Val{:Z},
	::Val{:re},
) where {S <: AbstractPlotSpec, dim} = QuantityTag{:resistance}()
axis_quantity(
	::Type{S},
	::Val{dim},
	::Val{:Z},
	::Val{:im},
) where {S <: AbstractPlotSpec, dim} = QuantityTag{:reactance}()

# Y: re/im correspond to G/B
axis_quantity(
	::Type{S},
	::Val{dim},
	::Val{:Y},
	::Val{:re},
) where {S <: AbstractPlotSpec, dim} = QuantityTag{:conductance}()
axis_quantity(
	::Type{S},
	::Val{dim},
	::Val{:Y},
	::Val{:im},
) where {S <: AbstractPlotSpec, dim} = QuantityTag{:susceptance}()

axis_quantity(
	::Type{S},
	::Val{dim},
	::Val{:Z},
	::Val{:abs},
) where {S <: AbstractPlotSpec, dim} = QuantityTag{(:impedance, :abs)}()
axis_quantity(
	::Type{S},
	::Val{dim},
	::Val{:Z},
	::Val{:angle},
) where {S <: AbstractPlotSpec, dim} = QuantityTag{(:impedance, :angle)}()

axis_quantity(
	::Type{S},
	::Val{dim},
	::Val{:Y},
	::Val{:abs},
) where {S <: AbstractPlotSpec, dim} = QuantityTag{(:admittance, :abs)}()
axis_quantity(
	::Type{S},
	::Val{dim},
	::Val{:Y},
	::Val{:angle},
) where {S <: AbstractPlotSpec, dim} = QuantityTag{(:admittance, :angle)}()
