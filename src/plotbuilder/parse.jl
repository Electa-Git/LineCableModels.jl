
# split_kwargs – purely “what did the user say?”
function split_kwargs(
	::Type{S},
	kwargs::NamedTuple,
	input_keys::Tuple,
	renderer_keys::Tuple,
	idx::Tuple,
	dims::Tuple,
) where {S <: AbstractPlotSpec}

	# AxisSpec selector keys: (:x, :y, :z)
	select_fields = dims

	semantic_keys = (input_keys..., idx..., select_fields...)
	allowed       = (semantic_keys..., renderer_keys...)

	spec_pairs = Tuple(filter(((k, _),) -> k in semantic_keys, pairs(kwargs)))
	renderer_pairs = Tuple(filter(((k, _),) -> k in renderer_keys, pairs(kwargs)))
	for k in keys(kwargs)
		k in allowed || @warn "Unknown plot keyword for $(S): :$(k)"
	end

	spec = NamedTuple(spec_pairs)
	renderer = NamedTuple(renderer_pairs)

	return spec, renderer
end

# merge_defaults – “how does this spec fill in the blanks?”
function merge_defaults(
	::Type{S},
	obj,
	spec::NamedTuple,
	renderer::NamedTuple,
) where {S <: AbstractPlotSpec}

	idefault = input_defaults(S, obj)
	bdefault = renderer_defaults(S, obj)

	spec_merged = merge(idefault, spec)
	renderer_merged = merge(bdefault, renderer)

	return spec_merged, renderer_merged
end

# normalize_indices – enforce Int vs range-capable
function normalize_indices(
	::Type{S},
	spec::NamedTuple,
	idx_keys::Tuple{Vararg{Symbol}},
	ranged_keys::Tuple{Vararg{Symbol}},
) where {S <: AbstractPlotSpec}

	# No index keys → nothing to normalize
	isempty(idx_keys) && return spec

	out = spec

	for k in idx_keys
		is_ranged = k in ranged_keys

		if is_ranged
			# Sample-like index (typically :k, optionally :l)
			# Default to full range when not provided at all.
			v = get(out, k, Colon())

			(v isa Int ||
			 v isa AbstractUnitRange{<:Int} ||
			 v isa Colon) ||
				Base.error(
					"Index $(k) for spec $(S) must be Int, AbstractUnitRange{<:Int} or `:`; " *
					"got $(typeof(v)).",
				)

			out = merge(out, NamedTuple{(k,)}((v,)))
		else
			# Selector indices (:i, :j, ...) – only normalized if explicitly present.
			# No defaults are invented for these.
			if haskey(out, k)
				v = getfield(out, k)
				v isa Int || Base.error(
					"Index $(k) for spec $(S) must be Int when provided; got $(typeof(v)).",
				)
				out = merge(out, NamedTuple{(k,)}((v,)))
			end
		end
	end

	return out
end


# sanity check selectors of datasources – ensure sources for xdata/ydata/... exist and are Symbols
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
				"AxisSpec selector $(d) must be a Symbol, got $(typeof(val)) for spec $(S)",
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

	# AxisSpec → datakey mapping (:x → :f, :y → :R, etc.)
	datakeys = Dict{Symbol, Symbol}()
	for d in dims
		datakeys[d] = getfield(spec, d)  # verified by verify_selectors
	end

	# Index presence semantics:
	# - i/j are selector indices: present iff field exists in spec NT
	# - k is sample-like only if it is in ranged_keys(S)
	rk = ranged_keys(S)

	has_i = (:i in idx_keys) && haskey(spec, :i)
	has_j = (:j in idx_keys) && haskey(spec, :j)
	has_k = (:k in rk)       # sample-like dimension iff ranged_keys(S) contains :k

	i_val = has_i ? spec.i : nothing
	j_val = has_j ? spec.j : nothing
	k_val = has_k ? spec.k : Colon()  # normalized in normalize_indices

	local function _check_index(name::Symbol, v, n::Int)
		v isa Int || Base.error(
			"Index $(name) must be Int, got $(typeof(v)) for spec $(S)",
		)
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
			Base.error("AxisSpec $(d) data for $(S) is scalar; expected an array.")

		# Enforce the same storage contract as axis_slice,
		# except that :x may be a global 1D vector.
		if has_i && has_j && !(d === :x && nd == 1) && nd < 3
			Base.error(
				"Invalid axis storage for $(d): spec uses indices :i and :j, " *
				"but container_array($(S), $(d)) returned an array with $(nd) dimension(s). " *
				"When both :i and :j are active, the underlying array must be at least 3D " *
				"(Ni, Nj, Nk...).",
			)
		end

		# Check i/j bounds using first/second dims when present.
		# Skip for global 1D :x vectors.
		if !(d === :x && nd == 1)
			if has_i && nd >= 1
				_check_index(:i, i_val, size(arr, 1))
			end
			if has_j && nd >= 2
				_check_index(:j, j_val, size(arr, 2))
			end
		end

		# Determine sample length along k or last dimension
		n_samp = if nd == 1
			length(arr)
		else
			size(arr, nd)
		end

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
			# No ranged k for this spec → sample length is the 1D length (nd == 1)
			# or the last dimension if nd ≥ 2; caller already ensured alignment.
			nd == 1 ? length(arr) : n_samp
		end

		lengths[d] = len

		# guard select_field semantics
		kfield = select_field(S, Val(d))
		if kfield !== nothing
			# select_field is interpreted strictly as a spec field name.
			# If that field is provided in the spec NT, then it must be a Symbol
			# and the data must be NamedTuple with that key.
			# If not provided, we only enforce that elements are NamedTuple;
			# the grammar decides how to use the keys later.
			isempty(arr) && continue

			first_el = first(arr)

			if kfield in keys(spec)
				v = spec[kfield]
				v isa Symbol || Base.error(
					"select_field($(S), Val($(d))) = :$(kfield) but spec.$(kfield) " *
					"is not a Symbol; got $(typeof(v)).",
				)
				sym = v

				first_el isa NamedTuple || Base.error(
					"Data for axis $(d) in $(S) must be NamedTuple when a leaf " *
					"field is selected via select_field; got $(typeof(first_el)).",
				)
				haskey(first_el, sym) || Base.error(
					"NamedTuple data for axis $(d) in $(S) has no key $(sym).",
				)
			else
				# select_field is defined but no concrete field has been bound yet.
				# Enforce that the data are NamedTuple; actual key usage is left
				# to the generic make_series/make_views logic.
				first_el isa NamedTuple || Base.error(
					"Data for axis $(d) in $(S) must be NamedTuple when " *
					"select_field($(S), Val($(d))) is defined; got $(typeof(first_el)).",
				)
			end
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

  1. Decide which kwargs matter for this spec (`input_kwargs`, `renderer_kwargs`,
	 `index_keys`, `geom_axes`) and partition user kwargs into semantic vs
	 renderer.
  2. Merge user kwargs with `input_defaults(S, obj)` and
	 `renderer_defaults(S, obj)`.
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
	bk_raw   = renderer_kwargs(S)
	idx_raw  = index_keys(S)
	dims_raw = geom_axes(S)
	rk_raw   = ranged_keys(S)

	# Coerce to tuples with warnings if someone was lazy
	ik   = trait_to_tuple(S, ik_raw, "input_kwargs")
	bk   = trait_to_tuple(S, bk_raw, "renderer_kwargs")
	idx  = trait_to_tuple(S, idx_raw, "index_keys")
	dims = trait_to_tuple(S, dims_raw, "geom_axes")
	rk   = trait_to_tuple(S, rk_raw, "ranged_keys")

	# Basic trait sanity: they should all be Symbols, and axes only from :x,:y,:z
	for (name, tup) in (("input_kwargs", ik), ("renderer_kwargs", bk),
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

	# Additional trait sanity for ranged_keys:
	# - ranged_keys ⊆ index_keys
	# - only :k and :l are allowed to be ranged (sample-like dims)
	if !isempty(rk)
		# ranged_keys ⊆ index_keys
		for key in rk
			key in idx || Base.error(
				"ranged_keys(::Type{$(S)}) includes $(key), which is not in " *
				"index_keys(::Type{$(S)}) = $(idx).",
			)
		end

		# Only :k and :l allowed as rangeable indices (sample dimensions)
		for key in rk
			(key === :k || key === :l) || Base.error(
				"ranged_keys(::Type{$(S)}) may only contain :k and/or :l. " *
				"Got $(key). Allowing :i or :j here would break the axis " *
				"semantics (sample dimension must remain unique).",
			)
		end
	end

	# 1) Split user kwargs into semantic vs backend
	spec_inputs, renderer_inputs = split_kwargs(S, kwargs, ik, bk, idx, dims)

	# 2) Merge with defaults
	spec_nt, renderer_nt = merge_defaults(S, obj, spec_inputs, renderer_inputs)

	# 3) Normalize indices (i,j,k,...) according to index/ranged traits
	spec_nt = normalize_indices(S, spec_nt, idx, rk)

	# 4) Ensure axis data keys exist and are Symbols
	verify_selectors(S, spec_nt, dims)

	# 5) Grammar-level structural sanity: containers, bounds, lengths
	verify_shapes(S, obj, spec_nt, dims, idx)

	return (; obj = obj, spec = spec_nt, renderer = renderer_nt)
end

# Convenience varargs wrapper
parse_kwargs(::Type{S}, obj; kwargs...) where {S <: AbstractPlotSpec} =
	parse_kwargs(S, obj, (; kwargs...))

"""
Resolve raw inputs into a normalized NamedTuple understood by `make_pages`.

This is where a spec implements its own mini-grammar:

- parse `values_expr` / `ijk`,
- pick matrix indices/slices,
- decide which quantities (R/L/C/G etc.) and which kind (:hist, :heatmap, ...).

Default is identity; spec types are expected to override.
"""
function resolve_input(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec}
	obj = nt.obj
	spec = nt.spec
	renderer_nt = nt.renderer

	dims = geom_axes(S)
	dims = dims isa Tuple ? dims : (dims,)

	xsel = :x in dims ? spec.x : nothing
	ysel = :y in dims ? spec.y : nothing
	zsel = :z in dims && haskey(spec, :z) ? spec.z : nothing

	# raw user knob (global, for now)
	has_as = haskey(spec, :as)
	raw_as = has_as ? spec.as : nothing

	xq = yq = zq = nothing
	xas = yas = zas = nothing

	any_complex = false

	if :x in dims
		if has_complex_qty(S, Val(:x), Val(xsel))
			any_complex = true
			xas = raw_as === nothing ? complex_as_default(S, Val(:x), Val(xsel)) : raw_as

			allowed = complex_as(S, Val(:x), Val(xsel))
			xas in allowed ||
				Base.error("Invalid as=$(xas) for x=$(xsel). Allowed: $(allowed).")

			xq = axis_quantity(S, Val(:x), Val(xsel), Val(xas))
		else
			xq = axis_quantity(S, Val(:x), Val(xsel))
		end
	end

	if :y in dims
		if has_complex_qty(S, Val(:y), Val(ysel))
			any_complex = true
			yas = raw_as === nothing ? complex_as_default(S, Val(:y), Val(ysel)) : raw_as
			allowed = complex_as(S, Val(:y), Val(ysel))
			yas in allowed ||
				Base.error("Invalid as=$(yas) for y=$(ysel). Allowed: $(allowed).")

			yq = axis_quantity(S, Val(:y), Val(ysel), Val(yas))
		else
			yq = axis_quantity(S, Val(:y), Val(ysel))
		end
	end

	if :z in dims && zsel !== nothing
		if has_complex_qty(S, Val(:z), Val(zsel))
			any_complex = true
			zas = raw_as === nothing ? complex_as_default(S, Val(:z), Val(zsel)) : raw_as

			allowed = complex_as(S, Val(:z), Val(zsel))
			zas in allowed ||
				Base.error("Invalid as=$(zas) for z=$(zsel). Allowed: $(allowed).")

			zq = axis_quantity(S, Val(:z), Val(zsel), Val(zas))
		else
			zq = axis_quantity(S, Val(:z), Val(zsel))
		end
	end

	# Tight API: if user asked for as= but nothing is complex, that's nonsense.
	if has_as && !any_complex
		Base.error(
			"Keyword as= is only valid for complex selectors with trait has_complex_qty == true.",
		)
	end

	out = spec

	if :x in dims
		out = merge(out, (; x = xsel, x_quantity = xq))
		xas === nothing || (out = merge(out, (; x_as = xas)))
	end
	if :y in dims
		out = merge(out, (; y = ysel, y_quantity = yq))
		yas === nothing || (out = merge(out, (; y_as = yas)))
	end
	if :z in dims && zsel !== nothing
		out = merge(out, (; z = zsel, z_quantity = zq))
		zas === nothing || (out = merge(out, (; z_as = zas)))
	end

	return merge(out, (; obj = obj, renderer = renderer_nt))
end
