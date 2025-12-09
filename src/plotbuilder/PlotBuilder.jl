module PlotBuilder

using Base: @kwdef

import ..UnitHandler: Units, QuantityTag, get_label, get_symbol, display_unit, scale_factor

export PlotAxis, PlotRenderer,
	plot_kind, enable_logscale, dispatch_on,
	axis_quantity, axis_unit, axis_label,
	default_figsize,
	parse_kwargs, resolve_input, build_payloads, make_renderer

abstract type AbstractPlotSpec end
include("plotspecs.jl")

# Submodule `BackendHandler`
include("backendhandler/BackendHandler.jl")
using .BackendHandler

# Submodule `PlotUIComponents`
include("plotuicomponents/PlotUIComponents.jl")
using .PlotUIComponents


# -----------------------------------------------------------------------------
# Axis semantics
# -----------------------------------------------------------------------------

"""
PlotAxis is the fully-decided axis description used by payloads.

- `dim`      : :x, :y, or :z
- `quantity` : QuantityTag{:freq}, QuantityTag{:R}, etc.
- `units`    : Units (from UnitHandler)
- `label`    : full label text, including unit symbol
- `scale`    : :linear or :log10
"""
struct PlotAxis
	dim::Symbol
	quantity::QuantityTag
	units::Units
	label::String
	scale::Symbol
end

"""
PlotRenderer ties together:

- a plot spec type `S <: AbstractPlotSpec` (the recipe / grammar),
- a concrete payload (NamedTuple) with all numeric data and Makie kwargs.

The render pipeline should only see `PlotRenderer` values and never touch
domain objects or validation logic.
"""
struct PlotRenderer{S <: AbstractPlotSpec, P <: NamedTuple}
	spec_type::Type{S}
	payload::P
end


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

"""
	build_axes(::Type{S}, nt::NamedTuple) where {S<:AbstractPlotSpec}

Build axes for spec `S` using quantity tags stored in `nt` as
fields `x_quantity`, `y_quantity`, `z_quantity` (when applicable).

Returns:
	(xaxis = PlotAxis or nothing,
	 yaxis = PlotAxis or nothing,
	 zaxis = PlotAxis or nothing)
"""
function build_axes(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec}
	dims = geom_axes(S)

	xaxis = nothing
	yaxis = nothing
	zaxis = nothing

	for dim in dims
		qfield = Symbol(dim, :_quantity)   # :x_quantity, :y_quantity, :z_quantity
		q = getproperty(nt, qfield)        # expected to be a QuantityTag

		u  = axis_unit(S, q, dim)
		ax = build_axis(S, dim, q, u)

		if dim === :x
			xaxis = ax
		elseif dim === :y
			yaxis = ax
		elseif dim === :z
			zaxis = ax
		else
			Base.error("Unsupported axis dim $(dim) in geom_axes for $(S)")
		end
	end

	return (xaxis = xaxis, yaxis = yaxis, zaxis = zaxis)
end

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
by `build_payloads(S, nt)`.
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

"""
Build the final axis label from quantity + units.

Default pattern: "Label [symbol]" if units are non-empty.
Uses `get_label(q::QuantityTag)` and `get_label(u::Units)` from UnitHandler.
"""
function axis_label(q::QuantityTag, u::Units)
	base = get_label(q)        # human-readable quantity name
	usym = get_label(u)        # unit symbol string, e.g. "Ω/km"
	return isempty(usym) ? base : string(base, " [", usym, "]")
end


"""
Internal: decide :linear or :log10 for axis dim `dim` in spec `S`,
based on `enable_logscale(S)`.
"""
function _axis_scale(::Type{S}, dim::Symbol) where {S <: AbstractPlotSpec}
	dims = enable_logscale(S)
	return dim in dims ? :log10 : :linear
end

"""
Build a PlotAxis for spec `S`, axis dim `dim`, quantity `q`, units `u`.
"""
function build_axis(::Type{S},
	dim::Symbol,
	q::QuantityTag,
	u::Units) where {S <: AbstractPlotSpec}

	lab = axis_label(q, u)
	sc  = _axis_scale(S, dim)
	return PlotAxis(dim, q, u, lab, sc)
end

# High-level builder: dim only (static semantics)
function build_axis(::Type{S}, ::Val{dim}) where {S <: AbstractPlotSpec, dim}
	q = axis_quantity(S, Val(dim))
	u = axis_unit(S, q, dim)
	return build_axis(S, dim, q, u)
end

# High-level builder: dim + semantic code (e.g. :R/:L/:C/:G)
function build_axis(
	::Type{S},
	::Val{dim},
	::Val{qty},
) where {S <: AbstractPlotSpec, dim, qty}
	q = axis_quantity(S, Val(dim), Val(qty))
	u = axis_unit(S, q, dim)
	return build_axis(S, dim, q, u)
end

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
	backend_kwargs(::Type{S}) where {S<:AbstractPlotSpec}

Backend / style kwargs that are simply forwarded to the plotting backend
(Makie today, whatever tomorrow).

Typical examples: :color, :linewidth, :marker, :markersize, :linestyle,
:colormap, :alpha, etc.
"""
backend_kwargs(::Type{S}) where {S <: AbstractPlotSpec} = ()

"""
Root container inside `obj` for the data of axis `dim`.

Default: `nothing` → use `obj` itself.

For example, if `obj.stats` is a NamedTuple of tensors and y-axis data
comes from there, define:

	data_container(::Type{MySpec}, ::Val{:y}) = :stats
"""
data_container(::Type{S}, ::Val{dim}) where {S <: AbstractPlotSpec, dim} = nothing

# Child key for axis `dim` inside the container returned by `data_container`.
# child_key = getfield(nt, axis_key(S, Val(dim)))  # e.g. :mean
# data      = obj.container.datakey[i,j,k].(child_key)
axis_key(::Type{S}, ::Val{dim}) where {S <: AbstractPlotSpec, dim} = nothing

"""
	input_defaults(::Type{S}, obj) where {S<:AbstractPlotSpec}

Defaults for semantic kwargs declared in `input_kwargs(S)`.

May depend on the dispatched object `obj` (e.g. pick default `:quantity`
from `obj` contents).
"""
input_defaults(::Type{S}, obj) where {S <: AbstractPlotSpec} = NamedTuple()

"""
	backend_defaults(::Type{S}, obj) where {S<:AbstractPlotSpec}

Defaults for backend kwargs declared in `backend_kwargs(S)`.

This is where a spec declares its default color/linestyle/whatever,
possibly depending on `obj` (e.g. per-phase colors).
"""
backend_defaults(::Type{S}, obj) where {S <: AbstractPlotSpec} = NamedTuple()

"""
How to group dataseries into figures.
Options: :none -> each dataseries spawns its own figure;
		 :overlay -> all dataseries in one figure, same axis;
		 :grid -> all dataseries in one figure, arranged in a grid.
"""
grouping_mode(::Type{S}) where {S <: AbstractPlotSpec} = :none

# -----------------------------------------------------------------------------

# split_kwargs – purely “what did the user say?”
function split_kwargs(
	::Type{S},
	kwargs::NamedTuple,
	ik::Tuple,
	bk::Tuple,
	idx::Tuple,
	dims::Tuple,
) where {S <: AbstractPlotSpec}

	# Axis selector keys: (:x, :y, :z)
	axis_keys = dims

	semantic_keys = (ik..., idx..., axis_keys...)
	allowed       = (semantic_keys..., bk...)

	spec_pairs    = Tuple(filter(((k, _),) -> k in semantic_keys, pairs(kwargs)))
	backend_pairs = Tuple(filter(((k, _),) -> k in bk, pairs(kwargs)))

	for k in keys(kwargs)
		k in allowed || @warn "Unknown plot keyword for $(S): :$(k)"
	end

	spec    = NamedTuple(spec_pairs)
	backend = NamedTuple(backend_pairs)

	return spec, backend
end

# merge_defaults – “how does this spec fill in the blanks?”
function merge_defaults(
	::Type{S},
	obj,
	spec::NamedTuple,
	backend::NamedTuple,
) where {S <: AbstractPlotSpec}

	idefault = input_defaults(S, obj)
	bdefault = backend_defaults(S, obj)

	spec_merged = merge(idefault, spec)
	backend_merged = merge(bdefault, backend)

	return spec_merged, backend_merged
end

# normalize_indices – enforce Int vs range-capable
function normalize_indices(
	::Type{S},
	spec::NamedTuple,
	idx::Tuple,
	rk::Tuple,
) where {S <: AbstractPlotSpec}

	# No index keys → nothing to normalize
	isempty(idx) && return spec

	idx_vals = Tuple(
		begin
			is_ranged = k in rk
			# Default: `:` for ranged (all), 1 for scalar indices
			default = is_ranged ? Colon() : 1

			v = get(spec, k, default)

			if is_ranged
				# Accept Int, AbstractUnitRange{<:Int}, or colon
				(v isa Int ||
				 v isa AbstractUnitRange{<:Int} ||
				 v isa Colon) ||
					Base.error(
						"Index $(k) must be Int, Int range, or `:`; got $(typeof(v))",
					)
			else
				v isa Int || Base.error("Index $(k) must be Int, got $(typeof(v))")
			end

			v
		end for k in idx
	)

	idx_nt = NamedTuple{idx}(idx_vals)
	return merge(spec, idx_nt)
end

# normalize selectors of datasources – ensure sources for xdata/ydata/... exist and are Symbols
function verify_selectors(
	::Type{S},
	spec::NamedTuple,
	dims::Tuple,
) where {S <: AbstractPlotSpec}

	for d in dims
		val = get(spec, d, nothing)
		val === nothing &&
			Base.error("Missing axis selector $(d) for spec $(S) after defaults")
		val isa Symbol ||
			Base.error(
				"Axis selector $(d) must be a Symbol, got $(typeof(val)) for spec $(S)",
			)
	end

	return
end

function container_array(
	::Type{S},
	obj,
	dim::Symbol,
	datakey::Symbol,
) where {S <: AbstractPlotSpec}

	container = data_container(S, Val(dim))

	if container === nothing
		hasproperty(obj, datakey) ||
			Base.error(
				"For spec $(S), axis $(dim) expects obj.$(datakey), " *
				"but $(typeof(obj)) has no such field.",
			)
		return getproperty(obj, datakey)
	else
		container isa Symbol ||
			Base.error(
				"data_container(::Type{$(S)}, Val($(dim))) must be Symbol or nothing; got $(typeof(container))",
			)

		hasproperty(obj, container) ||
			Base.error(
				"data_container(::Type{$(S)}, Val($(dim))) = :$(container), " *
				"but $(typeof(obj)) has no field :$(container).",
			)

		parent = getproperty(obj, container)

		if parent isa AbstractDict
			haskey(parent, datakey) ||
				Base.error(
					"Container field :$(container) for $(S) has no key :$(datakey) for axis $(dim).",
				)
			return parent[datakey]
		elseif parent isa NamedTuple && haskey(parent, datakey)
			return parent[datakey]
		elseif hasproperty(parent, datakey)
			return getproperty(parent, datakey)
		else
			try
				return parent[datakey]
			catch
				Base.error(
					"Container field :$(container) of type $(typeof(parent)) " *
					"does not provide data for key :$(datakey) for axis $(dim) in $(S).",
				)
			end
		end
	end
end

# normalize_shapes – centralized structural sanity
# container resolution,
# i/j bounds,
# sample length alignment.
# adjusted to respect the axis-level data_container(::Type{S}, ::Val{dim}) contract:
function verify_shapes(
	::Type{S},
	obj,
	spec::NamedTuple,
	dims::Tuple,
	idx_keys::Tuple,
) where {S <: AbstractPlotSpec}

	# Axis → datakey mapping (:x → :f, :y → :R, etc.)
	datakeys = Dict{Symbol, Symbol}()
	for d in dims
		datakeys[d] = getfield(spec, d)  # verified by verify_selectors
	end

	has_i = :i in idx_keys
	has_j = :j in idx_keys
	has_k = :k in idx_keys

	i_val = has_i ? spec.i : nothing
	j_val = has_j ? spec.j : nothing
	k_val = has_k ? spec.k : Colon()  # normalized to Int / range / Colon by normalize_indices

	local function _check_index(name::Symbol, v, n::Int)
		v isa Int || Base.error("Index $(name) must be Int, got $(typeof(v)) for spec $(S)")
		(1 <= v <= n) ||
			error("Index $(name) = $(v) out of bounds 1:$(n) for spec $(S)")
		v
	end

	lengths = Dict{Symbol, Int}()

	for d in dims
		datakey = datakeys[d]
		arr = container_array(S, obj, d, datakey)

		nd = ndims(arr)
		nd == 0 &&
			Base.error("Axis $(d) data for $(S) is scalar; expected an array.")

		# Check i/j bounds using first/second dims when present
		if has_i && nd >= 1
			_check_index(:i, i_val, size(arr, 1))
		end
		if has_j && nd >= 2
			_check_index(:j, j_val, size(arr, 2))
		end

		# Effective sample length along last dimension, accounting for k
		n_samp = size(arr, nd)
		n_samp == 0 && Base.error(
			"Axis $(d) data for $(S) has zero samples; no data to plot.",
		)
		len = if has_k
			kv = k_val
			if kv isa Int
				_check_index(:k, kv, n_samp)
				1
			elseif kv isa AbstractUnitRange{<:Int}
				first(kv) >= 1 && last(kv) <= n_samp ||
					error(
						"Range k = $(kv) out of bounds 1:$(n_samp) for spec $(S) on axis $(d).",
					)
				length(kv)
			elseif kv isa Colon
				n_samp
			else
				Base.error(
					"Index :k must be Int, Int range, or `:` after normalization; " *
					"got $(typeof(kv)) for spec $(S).",
				)
			end
		else
			nd == 1 ? length(arr) : n_samp
		end

		lengths[d] = len

		# guard axis_key semantics
		kfield = axis_key(S, Val(d))
		if kfield !== nothing
			# Decide how to interpret axis_key:
			# - if kfield is a field in spec → spec[kfield] must be a Symbol (indirect mode)
			# - else                  → kfield itself is the child key (direct mode)
			sym = if kfield in keys(spec)
				v = spec[kfield]
				v isa Symbol || Base.error(
					"axis_key($(S), Val($(d))) = :$(kfield) but spec.$(kfield) " *
					"is not a Symbol; got $(typeof(v)).",
				)
				v
			else
				kfield::Symbol
			end

			# Check data element type and key existence
			first_el = first(arr)
			first_el isa NamedTuple || Base.error(
				"Data for axis $(d) in $(S) must be NamedTuple when axis_key is used; " *
				"got $(typeof(first_el)).",
			)
			haskey(first_el, sym) || Base.error(
				"NamedTuple data for axis $(d) in $(S) has no key $(sym).",
			)
		end
	end

	vals = collect(values(lengths))
	isempty(vals) && return

	ref = first(vals)
	for (d, len) in lengths
		len == ref || Base.error(
			"Mismatched sample lengths for spec $(S): axis $(d) has length $(len), " *
			"expected $(ref). Containers must align along their sample dimension.",
		)
	end

	return
end




@inline function trait_to_tuple(::Type{S}, raw, name) where {S <: AbstractPlotSpec}
	raw === () && return ()
	raw isa Tuple && return raw
	@warn "Trait $(name) for $(S) should be a Tuple; got $(typeof(raw)). Coercing to 1-tuple."
	return (raw,)
end

"""
	parse_kwargs(::Type{S}, obj, kwargs::NamedTuple) where {S<:AbstractPlotSpec}

Grammar-level normalization phase.

Responsibilities:

  1. Decide which kwargs matter for this spec (`input_kwargs`, `backend_kwargs`,
	 `index_keys`, `geom_axes`) and partition user kwargs into semantic vs
	 backend.
  2. Merge user kwargs with `input_defaults(S, obj)` and
	 `backend_defaults(S, obj)`.
  3. Normalize indices (`index_keys(S)` / `ranged_keys(S)`) to the allowed
	 types and fill in defaults.
  4. Ensure axis data keys (`x`, `y`, ...) exist and are `Symbol`s.
  5. Run grammar-level structural checks:
	 - data sources exist under `obj` according to `data_container`,
	 - indices are in bounds,
	 - all active axes have compatible sample lengths.

Returns a canonical NamedTuple:

	(; obj = obj, spec = spec_nt, backend = backend_nt)

to be consumed by `resolve_input`.
"""
function parse_kwargs(::Type{S}, obj, kwargs::NamedTuple) where {S <: AbstractPlotSpec}
	# Raw traits
	ik_raw   = input_kwargs(S)
	bk_raw   = backend_kwargs(S)
	idx_raw  = index_keys(S)
	dims_raw = geom_axes(S)
	rk_raw   = ranged_keys(S)

	# Coerce to tuples with warnings if someone was lazy
	ik   = trait_to_tuple(S, ik_raw, "input_kwargs")
	bk   = trait_to_tuple(S, bk_raw, "backend_kwargs")
	idx  = trait_to_tuple(S, idx_raw, "index_keys")
	dims = trait_to_tuple(S, dims_raw, "geom_axes")
	rk   = trait_to_tuple(S, rk_raw, "ranged_keys")

	# Basic trait sanity: they should all be Symbols, and axes only from :x,:y,:z
	for (name, tup) in (("input_kwargs", ik), ("backend_kwargs", bk),
		("index_keys", idx), ("ranged_keys", rk))
		all(k -> k isa Symbol, tup) ||
			@warn "$(name)(::Type{$(S)}) should be a Tuple of Symbols, got $(tup)."
	end

	for d in dims
		d in (:x, :y, :z) ||
			Base.error(
				"geom_axes(::Type{$(S)}) returned unsupported axis $(d). " *
				"Valid axes are :x, :y, :z.",
			)
	end

	# 1) Split user kwargs into semantic vs backend
	spec_inputs, backend_inputs = split_kwargs(S, kwargs, ik, bk, idx, dims)

	# 2) Merge with defaults
	spec_nt, backend_nt = merge_defaults(S, obj, spec_inputs, backend_inputs)

	# 3) Normalize indices (i,j,k,...) according to index/ranged traits
	spec_nt = normalize_indices(S, spec_nt, idx, rk)

	# 4) Ensure axis data keys exist and are Symbols
	verify_selectors(S, spec_nt, dims)

	# 5) Grammar-level structural sanity: containers, bounds, lengths
	verify_shapes(S, obj, spec_nt, dims, idx)

	return (; obj = obj, spec = spec_nt, backend = backend_nt)
end

# Convenience varargs wrapper
parse_kwargs(::Type{S}, obj; kwargs...) where {S <: AbstractPlotSpec} =
	parse_kwargs(S, obj, (; kwargs...))

"""
Resolve raw inputs into a normalized NamedTuple understood by `build_payloads`.

This is where a spec implements its own mini-grammar:

- parse `values_expr` / `ijk`,
- pick matrix indices/slices,
- decide which quantities (R/L/C/G etc.) and which kind (:hist, :heatmap, ...).

Default is identity; spec types are expected to override.
"""
function resolve_input(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec}
	obj     = nt.obj
	spec    = nt.spec
	backend = nt.backend

	dims = geom_axes(S)
	dims = dims isa Tuple ? dims : (dims,)

	xsel = :x in dims ? spec.x : nothing
	ysel = :y in dims ? spec.y : nothing
	zsel = :z in dims && haskey(spec, :z) ? spec.z : nothing

	xq = nothing
	yq = nothing
	zq = nothing

	if :x in dims
		xq = axis_quantity(S, Val(:x), Val(xsel))
	end
	if :y in dims
		yq = axis_quantity(S, Val(:y), Val(ysel))
	end
	if :z in dims && zsel !== nothing
		zq = axis_quantity(S, Val(:z), Val(zsel))
	end

	out = spec

	if :x in dims
		out = merge(out, (; x = xsel, x_quantity = xq))
	end
	if :y in dims
		out = merge(out, (; y = ysel, y_quantity = yq))
	end
	if :z in dims && zsel !== nothing
		out = merge(out, (; z = zsel, z_quantity = zq))
	end

	return merge(out, (; obj = obj, backend = backend))
end

# --------------------------------------------------------------------------
# Axis-level transform hook
# --------------------------------------------------------------------------

"""
	axis_transform(::Type{S}, ::Val{dim}, ::Val{datakey}, nt, axis::PlotAxis, data) where {S,dim,datakey}

Per-spec hook to post-process the sliced axis data *before* unit scaling.

- `dim`      : :x, :y, or :z
- `datakey`  : axis selector symbol (e.g. :f, :R, :Z, ...)
- `nt`       : resolved input NamedTuple from `resolve_input`
- `axis`     : PlotAxis for this dimension (quantity + units + label + scale)
- `data`     : 1D numeric/complex array returned by `axis_slice`

Default is identity; specs override this to apply `abs`, `angle`, imperial
conversions, etc.
"""
axis_transform(
	::Type{S},
	::Val{dim},
	::Val{datakey},
	nt::NamedTuple,
	axis::PlotAxis,
	data,
) where {S <: AbstractPlotSpec, dim, datakey} = data

"""
	axis_slice(::Type{S}, nt, axis::PlotAxis, ::Val{dim}) where {S<:AbstractPlotSpec}

Return a 1D slice for axis `dim` using the grammar:

  * Use `data_container(S, Val(dim))` and the axis selector `nt.<dim>`
	(e.g. `nt.x`, `nt.y`) to locate the raw storage in `nt.obj`.
  * Apply indices `i, j` if present in `nt`, assuming the sample dimension
	is the last array dimension.
  * Optionally unwrap child fields using `axis_key(S, Val(dim))` if it is
	non-`nothing` and elements are NamedTuples.

No unit scaling and no numeric check happen here; those are handled by
`axis_transform` and `build_payloads`.
"""
function axis_slice(
	::Type{S},
	nt::NamedTuple,
	axis::PlotAxis,
	::Val{dim},
) where {S <: AbstractPlotSpec, dim}

	obj = nt.obj

	# Axis selector: what the user (or defaults) chose for this axis, e.g. :f, :R, ...
	selector = getfield(nt, dim)::Symbol

	# --- Fetch raw array via centralized container logic ---
	raw_arr = container_array(S, obj, dim, selector)

	# --- Apply indices (i,j,k) → 1D slice along sample dimension ---
	arr = raw_arr
	nd  = ndims(arr)

	has_i = haskey(nt, :i)
	has_j = haskey(nt, :j)
	has_k = haskey(nt, :k)

	# First slice in i,j where applicable
	if has_i && has_j && nd >= 3
		arr = view(arr, nt.i, nt.j, :)
	elseif has_i && nd >= 2 && !has_j
		arr = view(arr, nt.i, :)
	elseif has_j && nd >= 2 && !has_i
		arr = view(arr, :, nt.j)
	end

	# Then slice in k along last dimension (sample dim)
	if has_k
		k = nt.k
		nd2 = ndims(arr)

		if nd2 == 0
			Base.error(
				"Axis $(dim) for $(S) has scalar data after i/j slicing; cannot apply k index.",
			)
		end

		if nd2 == 1
			if k isa Int
				arr = view(arr, k:k)
			elseif k isa AbstractUnitRange{<:Int} || k isa Colon
				arr = view(arr, k)
			else
				Base.error(
					"Index :k must be Int, Int range, or `:` after normalization; " *
					"got $(typeof(k)) for spec $(S) on axis $(dim).",
				)
			end
		else
			# nd2 ≥ 2, index last dimension
			lastdim = nd2
			if k isa Int
				inds = ntuple(d -> d == lastdim ? (k:k) : Colon(), lastdim)
			elseif k isa AbstractUnitRange{<:Int} || k isa Colon
				inds = ntuple(d -> d == lastdim ? k : Colon(), lastdim)
			else
				Base.error(
					"Index :k must be Int, Int range, or `:` after normalization; " *
					"got $(typeof(k)) for spec $(S) on axis $(dim).",
				)
			end
			arr = view(arr, inds...)
		end
	end

	ndims(arr) == 1 ||
		Base.error(
			"Axis $(dim) for $(S) expected to resolve to a 1D slice after indexing; " *
			"got $(ndims(arr))-dimensional array.",
		)

	vec_arr = arr

	# --- NamedTuple unwrapping via axis_key ---
	kfield = axis_key(S, Val(dim))

	if kfield === nothing
		return collect(vec_arr)
	else
		# Same dual semantics:
		# - if nt has a field kfield → that field must be a Symbol (indirect mode)
		# - else                    → kfield itself is the child key (direct mode)
		ksym = if kfield in propertynames(nt)
			v = getfield(nt, kfield)
			v isa Symbol || Base.error(
				"Field $(kfield) in resolved input for $(S) on axis $(dim) " *
				"must be a Symbol; got $(typeof(v)).",
			)
			v
		else
			kfield::Symbol
		end

		first_el = first(vec_arr)
		first_el isa NamedTuple ||
			Base.error(
				"axis_key($(S), Val($(dim))) expects NamedTuple elements; " *
				"got $(typeof(first_el)).",
			)

		haskey(first_el, ksym) || Base.error(
			"NamedTuple elements on axis $(dim) for $(S) have no key $(ksym). " *
			"Available keys: $(collect(keys(first_el))).",
		)

		return [el[ksym] for el in vec_arr]
	end

end


# Process one axis if present
@inline function axis_data(
	::Type{S},
	dim::Symbol,
	nt::NamedTuple,
	axis::Union{PlotAxis, Nothing},
) where {S <: AbstractPlotSpec}
	axis === nothing && return nothing

	# axis selector: nt.x / nt.y / nt.z
	selector = getfield(nt, dim)::Symbol

	# 1) slice + axis_key unwrapping (no scaling)
	raw_vec = axis_slice(S, nt, axis, Val(dim))

	# 2) spec-level transform
	transformed = axis_transform(S, Val(dim), Val(selector), nt, axis, raw_vec)

	# 3) numeric check
	transformed isa AbstractArray ||
		Base.error(
			"Axis $(dim) for $(S) did not resolve to an array; got $(typeof(transformed)).",
		)

	eltype(transformed) <: Number ||
		Base.error(
			"Axis $(dim) for $(S) did not resolve to numeric data; got eltype $(eltype(transformed)).",
		)

	# 4) unit scaling
	sf = scale_factor(axis.quantity, axis.units)
	return sf .* transformed
end

"""
Final construction step: from normalized NamedTuple to concrete payloads.

Return a `Vector{<:NamedTuple}` where each NamedTuple is the full Makie
payload for a *single* plot primitive, with keys as close as possible
to Makie’s own vocabulary, e.g.:

	(
		xdata   = ::AbstractVector,
		ydata   = ::AbstractVector,
		zdata   = ::Union{AbstractArray,Nothing},
		xaxis   = ::PlotAxis,
		yaxis   = ::PlotAxis,
		zaxis   = ::Union{PlotAxis,Nothing},
		title   = ::String,
		legend  = ::Vector{String},
		kwargs  = ::NamedTuple,   # color, linestyle, markersize, ...
	)

Specs may override this when they need multiple primitives, grids, etc.
"""
# --------------------------------------------------------------------------
# Generic payload builder
# --------------------------------------------------------------------------

function build_payloads(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec}
	dims    = geom_axes(S)
	dims    = dims isa Tuple ? dims : (dims,)
	backend = nt.backend

	axes  = build_axes(S, nt)
	xaxis = axes.xaxis
	yaxis = axes.yaxis
	zaxis = axes.zaxis

	xdata = nothing
	ydata = nothing
	zdata = nothing

	:x in dims && (xdata = axis_data(S, :x, nt, xaxis))
	:y in dims && (ydata = axis_data(S, :y, nt, yaxis))
	:z in dims && (zdata = axis_data(S, :z, nt, zaxis))

	title  = default_title(S, nt)
	legend = legend_labels(S, nt)

	payload = (
		xdata  = xdata,
		ydata  = ydata,
		zdata  = zdata,
		xaxis  = xaxis,
		yaxis  = yaxis,
		zaxis  = zaxis,
		title  = title,
		legend = legend,
		kwargs = backend,
	)

	return [payload]
end

"""
High-level API: from domain object + kwargs to a vector of PlotRenderer.

Checks that the object type is compatible with `dispatch_on(S)` and then
runs:

	parse_kwargs → resolve_input → build_payloads

`build_payloads` returns vector of payload NamedTuples; `make_renderer` wraps them into
`PlotRenderer` values that the UI layer will later assemble into actual
windows/layouts.
"""
function make_renderer(::Type{S}, obj; kwargs...) where {S <: AbstractPlotSpec}
	Tdispatch = dispatch_on(S)
	obj isa Tdispatch ||
		Base.error("Spec $(S) cannot dispatch on $(typeof(obj)); expected $(Tdispatch)")

	raw      = parse_kwargs(S, obj; kwargs...)
	norm     = resolve_input(S, raw)
	payloads = build_payloads(S, norm)  # ::Vector{<:NamedTuple}

	renderers = Vector{PlotRenderer}(undef, length(payloads))
	@inbounds for (k, p) in pairs(payloads)
		renderers[k] = PlotRenderer(S, p)
	end
	return renderers
end

end # module PlotBuilder
