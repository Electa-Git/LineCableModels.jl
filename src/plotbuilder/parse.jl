
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
Resolve raw inputs into a normalized NamedTuple understood by `build_figures`.

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