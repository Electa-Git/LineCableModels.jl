module PlotBuilder

using Base: @kwdef

import ..UnitHandler: Units, QuantityTag, get_label, get_symbol, display_unit

export PlotAxis, PlotRenderer,
	plot_kind, enable_logscale, dispatch_on,
	axis_quantity, axis_unit, axis_label,
	default_figsize,
	parse_kwargs, resolve_input, build_payloads, make

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
			error("Unsupported axis dim $(dim) in geom_axes for $(S)")
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

	# Axis data keys: (:x, :y, :z) → (:xdata, :ydata, :zdata)
	data_keys = Tuple(Symbol(d, :data) for d in dims)

	semantic_keys = (ik..., idx..., data_keys...)
	allowed       = (semantic_keys..., bk...)

	# Split kwargs into semantic vs backend
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
					error(
						"Index $(k) must be Int, Int range, or `:`; got $(typeof(v))",
					)
			else
				v isa Int || error("Index $(k) must be Int, got $(typeof(v))")
			end

			v
		end for k in idx
	)

	idx_nt = NamedTuple{idx}(idx_vals)
	return merge(spec, idx_nt)
end

# normalize_datasources – ensure xdata/ydata/... exist and are Symbols
function verify_datasources(
	::Type{S},
	spec::NamedTuple,
	dims::Tuple,
) where {S <: AbstractPlotSpec}

	for d in dims
		key = Symbol(d, :data) # :xdata, :ydata, :zdata
		val = get(spec, key, nothing)
		val === nothing &&
			Base.error("Missing data source $(key) for spec $(S) after defaults")
		val isa Symbol ||
			Base.error("$(key) must be a Symbol, got $(typeof(val)) for spec $(S)")
	end

	return
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
		key = Symbol(d, :data) # :xdata, :ydata, :zdata
		val = get(spec, key, nothing)
		datakeys[d] = val
	end

	# Extract indices if present
	has_i = :i in idx_keys
	has_j = :j in idx_keys
	i_val = has_i ? spec.i : nothing
	j_val = has_j ? spec.j : nothing

	# Helper: basic index bounds check (Int only for now; ranges will need
	# additional logic if/when you allow them for i/j).
	local function _check_index(name::Symbol, v, n::Int)
		v isa Int || Base.error("Index $(name) must be Int, got $(typeof(v)) for spec $(S)")
		(1 <= v <= n) ||
			error("Index $(name) = $(v) out of bounds 1:$(n) for spec $(S)")
		v
	end

	# Helper: resolve the raw storage array for a given axis + data source key.
	# This is where `data_container(::Type{S}, ::Val{dim})` is obeyed.
	local function _data_array(dim::Symbol, datakey::Symbol)
		container = data_container(S, Val(dim))

		if container === nothing
			# Direct: obj.R, obj.f, ...
			return getproperty(obj, datakey)
		else
			# Indirect: obj.stats[:R], obj.some_nt[:foo], ...
			parent = getproperty(obj, container)
			return parent[datakey]
		end
	end

	# Compute sample length per axis (last dimension or vector length)
	lengths = Dict{Symbol, Int}()

	for d in dims
		datakey = datakeys[d]
		arr = _data_array(d, datakey)

		nd = ndims(arr)
		nd == 0 && error("Axis $(d) data for $(S) is scalar; expected an array.")

		# Check i/j bounds when there are at least that many dims
		if has_i && nd >= 1
			_check_index(:i, i_val, size(arr, 1))
		end
		if has_j && nd >= 2
			_check_index(:j, j_val, size(arr, 2))
		end

		# Sample length = length of vector (1D) or size along last dimension (≥2D)
		len = nd == 1 ? length(arr) : size(arr, nd)
		lengths[d] = len
	end

	# Ensure all axes share the same sample length
	vals = collect(values(lengths))
	isempty(vals) && return

	ref = first(vals)
	for (d, len) in lengths
		len == ref || error(
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
  4. Ensure axis data keys (`xdata`, `ydata`, ...) exist and are `Symbol`s.
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
			error(
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
	verify_datasources(S, spec_nt, dims)

	# 5) Grammar-level structural sanity: containers, bounds, lengths
	verify_shapes(S, obj, spec_nt, dims, idx)

	return (; obj = obj, spec = spec_nt, backend = backend_nt)
end

# Convenience varargs wrapper
parse_kwargs(::Type{S}, obj; kwargs...) where {S <: AbstractPlotSpec} =
	parse_kwargs(S, obj, (; kwargs...))


# # -----------------------------------------------------------------------------
# # build_payloads pipeline (spec → payloads)
# # -----------------------------------------------------------------------------
# """
# 	parse_kwargs(::Type{S}, obj, kwargs::NamedTuple) where {S<:AbstractPlotSpec}
# 	parse_kwargs(::Type{S}, obj; kwargs...) where {S<:AbstractPlotSpec}

# Normalize user keyword arguments for a plot specification `S` into a canonical
# NamedTuple with three fields: `obj`, `spec`, and `style`.

# This function implements the trait-driven plot grammar. It splits the incoming
# `kwargs` into:

#   • semantic control parameters declared by `input_kwargs(S)`, merged with
# 	`input_defaults(S, obj)` into `spec`;
#   • backend or style parameters declared by `backend_kwargs(S)`, merged with
# 	`backend_defaults(S, obj)` into `style`.

# Any keyword not listed in `input_kwargs(S)` or `backend_kwargs(S)` is ignored
# for the merge and may trigger a warning. No other filtering, mutation, or
# validation is performed here; downstream methods are expected to consume the
# canonical form only.

# The returned NamedTuple is the sole input shape that `resolve_input` and
# `build_payloads` are expected to handle for specification `S`.

# # Arguments

#   • `S::Type{<:AbstractPlotSpec}`: plot specification type.
#   • `obj`: domain object dispatched by `S` \\(for example, a Monte Carlo result
# 	container\\).
#   • `kwargs::NamedTuple`: user keyword arguments collected from the plotting
# 	call.

# # Returns

#   • `NamedTuple`: a normalized tuple

# 	  (; obj = obj, spec = spec_nt, style = style_nt)

# 	where `spec_nt` and `style_nt` are NamedTuples obtained by partitioning
# 	and merging `kwargs` with the trait-provided defaults.

# # See also

# `input_kwargs`, `input_defaults`, `backend_kwargs`, `backend_defaults`,
# `resolve_input`, `build_payloads`.
# """
# function parse_kwargs(::Type{S}, obj, kwargs::NamedTuple) where {S <: AbstractPlotSpec}
# 	# Raw trait values
# 	ik_raw   = input_kwargs(S)
# 	bk_raw   = backend_kwargs(S)
# 	idx_raw  = index_keys(S)
# 	dims_raw = geom_axes(S)
# 	rk_raw   = ranged_keys(S)

# 	# Coerce to tuples, with warnings if dev was lazy.
# 	ik =
# 		ik_raw === () ? () :
# 		ik_raw isa Tuple ? ik_raw :
# 		begin
# 			@warn "input_kwargs(::Type{$(S)}) should return a Tuple of Symbols. " *
# 				  "Wrapping single value `$(ik_raw)` into a 1-tuple. Use (:key,) instead of :key."
# 			(ik_raw,)
# 		end

# 	bk =
# 		bk_raw === () ? () :
# 		bk_raw isa Tuple ? bk_raw :
# 		begin
# 			@warn "backend_kwargs(::Type{$(S)}) should return a Tuple of Symbols. " *
# 				  "Wrapping single value `$(bk_raw)` into a 1-tuple. Use (:key,) instead of :key."
# 			(bk_raw,)
# 		end

# 	idx =
# 		idx_raw === () ? () :
# 		idx_raw isa Tuple ? idx_raw :
# 		begin
# 			@warn "index_keys(::Type{$(S)}) should return a Tuple of Symbols. " *
# 				  "Wrapping single value `$(idx_raw)` into a 1-tuple. Use (:i,) or (:i,:j)."
# 			(idx_raw,)
# 		end

# 	dims =
# 		dims_raw === () ? () :
# 		dims_raw isa Tuple ? dims_raw :
# 		begin
# 			@warn "geom_axes(::Type{$(S)}) should return a Tuple of axis symbols, e.g. (:x, :y). " *
# 				  "Wrapping single value `$(dims_raw)` into a 1-tuple."
# 			(dims_raw,)
# 		end

# 	rk =
# 		rk_raw === () ? () :
# 		rk_raw isa Tuple ? rk_raw :
# 		begin
# 			@warn "ranged_keys(::Type{$(S)}) should return a Tuple of Symbols. " *
# 				  "Wrapping single value `$(rk_raw)` into a 1-tuple. Use (:k,) instead of :k."
# 			(rk_raw,)
# 		end

# 	# Optionally sanity-check symbol-ness
# 	if !all(d -> d isa Symbol, dims)
# 		@warn "geom_axes(::Type{$(S)}) returned non-Symbol entries $(dims). Expected dims like (:x, :y, :z)."
# 	end
# 	for (name, tup) in (("input_kwargs", ik), ("backend_kwargs", bk),
# 		("index_keys", idx), ("ranged_keys", rk))
# 		all(k -> k isa Symbol, tup) ||
# 			@warn "$(name)(::Type{$(S)}) should be a Tuple of Symbols, got $(tup)."
# 	end

# 	# For each dim (:x,:y,:z) define corresponding data key (:xdata,:ydata,:zdata)
# 	data_keys = Tuple(Symbol(d, :data) for d in dims)

# 	semantic_keys = (ik..., idx..., data_keys...)

# 	allowed  = (semantic_keys..., bk...)
# 	idefault = input_defaults(S, obj)
# 	bdefault = backend_defaults(S, obj)

# 	spec_pairs    = Tuple(filter(((k, _),)->k in semantic_keys, pairs(kwargs)))
# 	backend_pairs = Tuple(filter(((k, _),)->k in bk, pairs(kwargs)))

# 	for k in keys(kwargs)
# 		k in allowed || @warn "Unknown plot keyword for $(S): :$(k)"
# 	end

# 	spec_user    = NamedTuple(spec_pairs)
# 	backend_user = NamedTuple(backend_pairs)

# 	spec0 = merge(idefault, spec_user)
# 	style = merge(bdefault, backend_user)

# 	# Enforce index semantics: Int vs range-capable
# 	if isempty(idx)
# 		spec = spec0
# 	else
# 		idx_vals = Tuple(
# 			begin
# 				is_ranged = k in rk
# 				default = is_ranged ? Colon() : 1   # `:` means "all" for ranged axes

# 				v = get(spec0, k, default)

# 				if is_ranged
# 					# Accept Int, AbstractUnitRange{<:Int}, or colon
# 					(v isa Int ||
# 					 v isa AbstractUnitRange{<:Int} ||
# 					 v === Colon()) ||
# 						error(
# 							"Index $(k) must be Int, Int range, or `:`; got $(typeof(v))",
# 						)
# 				else
# 					v isa Int || error("Index $(k) must be Int, got $(typeof(v))")
# 				end

# 				v
# 			end for k in idx
# 		)

# 		idx_nt = NamedTuple{idx}(idx_vals)
# 		spec   = merge(spec0, idx_nt)
# 	end

# 	# Sanity check: enforce basic grammar invariants after normalization.
# 	#
# 	# Invariants:
# 	#   - For each geometric axis dim ∈ geom_axes(S),
# 	#       there must be a data source key :xdata/:ydata/:zdata in `spec`.
# 	#   - For each axis, the resolved storage must be an array:
# 	#       * 1D: treated as a sample vector.
# 	#       * ≥2D: last dimension is the sample axis.
# 	#   - All active axes must have the same sample length.
# 	#   - Index keys (:i,:j,...) must be in bounds for any storage that uses
# 	#     the corresponding dimension (we assume dim1→i, dim2→j as per your
# 	#     cable grammar).
# 	#
# 	# No content/physics is inspected; this only checks shapes and bounds.
# 	dims = geom_axes(S)
# 	dims = dims isa Tuple ? dims : (dims,)

# 	idx_keys = index_keys(S)
# 	idx_keys = idx_keys isa Tuple ? idx_keys : (idx_keys,)

# 	# Axis → datakey mapping (:x → :f, :y → :R, etc.)
# 	datakeys = Dict{Symbol, Symbol}()

# 	for d in dims
# 		key = Symbol(d, :data) # :xdata, :ydata, :zdata
# 		val = get(spec, key, nothing)
# 		val === nothing && error("Missing data source $(key) for spec $(S)")
# 		val isa Symbol || error("$(key) must be a Symbol, got $(typeof(val))")
# 		datakeys[d] = val
# 	end

# 	# Extract indices if present
# 	has_i = :i in idx_keys
# 	has_j = :j in idx_keys
# 	i_val = has_i ? spec.i : nothing
# 	j_val = has_j ? spec.j : nothing

# 	# Helper: basic index bounds check (Int only for now; ranges handled later if/when you add :k)
# 	_check_index(name::Symbol, v, n::Int) = begin
# 		v isa Int || error("Index $(name) must be Int, got $(typeof(v))")
# 		(1 <= v <= n) ||
# 			error("Index $(name) = $(v) out of bounds 1:$(n) for spec $(S)")
# 		v
# 	end

# 	# Helper: resolve the raw storage array for a given data source key.
# 	# Obeys the `data_container` trait.
# 	function _data_array(::Type{S}, obj, datakey::Symbol) where {S <: AbstractPlotSpec}
# 		container = data_container(S, Val(datakey))

# 		if container === nothing
# 			# Direct: obj.R, obj.f, ...
# 			return getproperty(obj, datakey)
# 		else
# 			# Indirect: obj.stats[:R], obj.some_nt[:foo], ...
# 			parent = getproperty(obj, container)
# 			return parent[datakey]
# 		end
# 	end

# 	# Compute sample length per axis (last dimension or vector length)
# 	lengths = Dict{Symbol, Int}()

# 	for d in dims
# 		datakey = datakeys[d]
# 		arr = _data_array(S, obj, datakey)

# 		nd = ndims(arr)
# 		nd == 0 && error("Axis $(d) data for $(S) is scalar; expected an array.")

# 		# Check i/j bounds when there are at least that many dims
# 		if has_i && nd >= 1
# 			_check_index(:i, i_val, size(arr, 1))
# 		end
# 		if has_j && nd >= 2
# 			_check_index(:j, j_val, size(arr, 2))
# 		end

# 		# Sample length = length of vector (1D) or size along last dimension (≥2D)
# 		len = nd == 1 ? length(arr) : size(arr, nd)

# 		lengths[d] = len
# 	end

# 	# Ensure all axes share the same sample length
# 	vals = collect(values(lengths))
# 	isempty(vals) && return

# 	ref = first(vals)
# 	for (d, len) in lengths
# 		len == ref || error(
# 			"Mismatched sample lengths for spec $(S): axis $(d) has length $(len), " *
# 			"expected $(ref). Containers must align along their sample dimension.",
# 		)
# 	end

# 	return (; obj = obj, spec = spec, style = style)
# end


# parse_kwargs(::Type{S}, obj; kwargs...) where {S <: AbstractPlotSpec} =
# 	parse_kwargs(S, obj, (; kwargs...))


# """
# Resolve raw inputs into a normalized NamedTuple understood by `build_payloads`.

# This is where a spec implements its own mini-grammar:

# - parse `values_expr` / `ijk`,
# - pick matrix indices/slices,
# - decide which quantities (R/L/C/G etc.) and which kind (:hist, :heatmap, ...).

# Default is identity; spec types are expected to override.
# """
# function resolve_input(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec}
# 	obj   = nt.obj
# 	spec  = nt.spec
# 	style = nt.style

# 	# Geometry + index traits
# 	dims = geom_axes(S)
# 	dims = dims isa Tuple ? dims : (dims,)

# 	idx_keys = index_keys(S)
# 	idx_keys = idx_keys isa Tuple ? idx_keys : (idx_keys,)

# 	# ------------------------------------------------------------------
# 	# Axis data sources per dim (xdata, ydata, zdata)
# 	# ------------------------------------------------------------------
# 	xdata_src = nothing
# 	ydata_src = nothing
# 	zdata_src = nothing

# 	for d in dims
# 		key = Symbol(d, :data)        # :xdata, :ydata, :zdata
# 		val = get(spec, key, nothing)
# 		val === nothing && error("Missing required data source $(key) for spec $(S)")
# 		val isa Symbol || error("$(key) must be a Symbol, got $(typeof(val))")

# 		if d === :x
# 			xdata_src = val
# 		elseif d === :y
# 			ydata_src = val
# 		elseif d === :z
# 			zdata_src = val
# 		else
# 			error("Unsupported axis dim $(d) in geom_axes for $(S)")
# 		end
# 	end

# 	# ------------------------------------------------------------------
# 	# Axis quantity semantics (only for existing axes)
# 	# ------------------------------------------------------------------
# 	xq = nothing
# 	yq = nothing
# 	zq = nothing

# 	if :x in dims
# 		xq = axis_quantity(S, Val(:x), Val(xdata_src))
# 	end
# 	if :y in dims
# 		yq = axis_quantity(S, Val(:y), Val(ydata_src))
# 	end
# 	if :z in dims
# 		zq = axis_quantity(S, Val(:z), Val(zdata_src))
# 	end

# 	# ------------------------------------------------------------------
# 	# Flatten spec into top-level NT + attach axis data & quantities
# 	# ------------------------------------------------------------------
# 	out = spec

# 	if :x in dims
# 		out = merge(out, (; xdata = xdata_src, x_quantity = xq))
# 	end
# 	if :y in dims
# 		out = merge(out, (; ydata = ydata_src, y_quantity = yq))
# 	end
# 	if :z in dims
# 		out = merge(out, (; zdata = zdata_src, z_quantity = zq))
# 	end

# 	# Attach object and style
# 	out = merge(out, (; obj = obj, style = style))

# 	return out
# end

"""
Resolve raw inputs into a normalized NamedTuple understood by `build_payloads`.

This is where a spec implements its own mini-grammar:

- parse `values_expr` / `ijk`,
- pick matrix indices/slices,
- decide which quantities (R/L/C/G etc.) and which kind (:hist, :heatmap, ...).

Default is identity; spec types are expected to override.
"""
function resolve_input(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec}
	obj = nt.obj
	spec = nt.spec
	backend_nt = nt.backend

	# Geometry trait
	dims = geom_axes(S)
	dims = dims isa Tuple ? dims : (dims,)

	# ------------------------------------------------------------------
	# Axis data sources (symbols) — invariants guaranteed by parse_kwargs
	# ------------------------------------------------------------------
	xdata_src = :x in dims ? spec.xdata : nothing
	ydata_src = :y in dims ? spec.ydata : nothing
	zdata_src = :z in dims && haskey(spec, :zdata) ? spec.zdata : nothing

	# ------------------------------------------------------------------
	# Axis quantity semantics (only for existing axes)
	# ------------------------------------------------------------------
	xq = nothing
	yq = nothing
	zq = nothing

	if :x in dims
		xq = axis_quantity(S, Val(:x), Val(xdata_src))
	end
	if :y in dims
		yq = axis_quantity(S, Val(:y), Val(ydata_src))
	end
	if :z in dims && zdata_src !== nothing
		zq = axis_quantity(S, Val(:z), Val(zdata_src))
	end

	# ------------------------------------------------------------------
	# Flatten spec into top-level NT + attach axis data & quantities
	# ------------------------------------------------------------------
	out = spec

	if :x in dims
		out = merge(out, (; xdata = xdata_src, x_quantity = xq))
	end
	if :y in dims
		out = merge(out, (; ydata = ydata_src, y_quantity = yq))
	end
	if :z in dims && zdata_src !== nothing
		out = merge(out, (; zdata = zdata_src, z_quantity = zq))
	end

	# Attach object and style
	return merge(out, (; obj = obj, backend = backend_nt))
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

Every concrete spec type MUST override this.
"""
function build_payloads(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec}
	error("`build_payloads(::Type{$(S)}, nt)` not implemented for this plot spec")
end

"""
High-level API: from domain object + kwargs to a vector of PlotRenderer.

Checks that the object type is compatible with `dispatch_on(S)` and then
runs:

	parse_kwargs → resolve_input → build_payloads

`build_payloads` returns vector of payload NamedTuples; `make` wraps them into
`PlotRenderer` values that the UI layer will later assemble into actual
windows/layouts.
"""
function make(::Type{S}, obj; kwargs...) where {S <: AbstractPlotSpec}
	Tdispatch = dispatch_on(S)
	obj isa Tdispatch ||
		error("Spec $(S) cannot dispatch on $(typeof(obj)); expected $(Tdispatch)")

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
