
# --------------------------------------------------------------------------
# Axis-level transform hook
# --------------------------------------------------------------------------

"""
	axis_transform(::Type{S}, ::Val{dim}, ::Val{datakey}, nt, axis::Axis, data) where {S,dim,datakey}

Per-spec hook to post-process the sliced axis data *before* unit scaling.

- `dim`      : :x, :y, or :z
- `datakey`  : axis selector symbol (e.g. :f, :R, :Z, ...)
- `nt`       : resolved input NamedTuple from `resolve_input`
- `axis`     : Axis for this dimension (quantity + units + label + scale)
- `data`     : 1D numeric/complex array returned by `axis_slice`

Default is identity; specs override this to apply `abs`, `angle`, imperial
conversions, etc.
"""
axis_transform(
	::Type{S},
	::Val{dim},
	::Val{datakey},
	nt::NamedTuple,
	axis::Axis,
	data,
) where {S <: AbstractPlotSpec, dim, datakey} = data

"""
	axis_slice(::Type{S}, nt, axis::Axis, ::Val{dim}) where {S<:AbstractPlotSpec}

Return a 1D slice for axis `dim` using the grammar:

  * Use `data_container(S, Val(dim))` and the axis selector `nt.<dim>`
	(e.g. `nt.x`, `nt.y`) to locate the raw storage in `nt.obj`.
  * Apply indices `i, j` if present in `nt`, assuming the sample dimension
	is the last array dimension.
  * Optionally unwrap child fields using `select_field(S, Val(dim))` if it is
	non-`nothing` and elements are NamedTuples.

No unit scaling and no numeric check happen here; those are handled by
`axis_transform` and `build_series`.
"""
function axis_slice(
	::Type{S},
	nt::NamedTuple,
	axis::Axis,
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
	@show nd
	@show arr
	has_i = haskey(nt, :i)
	has_j = haskey(nt, :j)
	has_k = haskey(nt, :k)

	# First slice in i,j where applicable
	if has_i && has_j
		if nd < 3
			Base.error(
				"Invalid axis storage for $(dim): spec uses indices :i and :j, " *
				"but container_array($(S), $(dim)) returned an array with $(nd) dimension(s). " *
				"When both :i and :j are active, the underlying array must be at least 3D " *
				"(Ni, Nj, Nk...). Check index_keys($(S)) and container_array($(S), $(dim)).",
			)
		end
		# canonical case: Ni×Nj×Nk...
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

	# --- NamedTuple unwrapping via select_field ---
	kfield = select_field(S, Val(dim))

	if kfield === nothing
		return collect(vec_arr)
	else
		# select_field is interpreted strictly as a spec field name that, when present
		# in the resolved input `nt`, holds the Symbol of the NamedTuple field to
		# extract. If the field is not present, we *do not* guess: we simply
		# return the NamedTuple vector and let higher-level grammar decide what
		# to do (overlay all fields, facet, etc.).
		if haskey(nt, kfield)
			v = getfield(nt, kfield)
			v isa Symbol || Base.error(
				"Field $(kfield) in resolved input for $(S) on axis $(dim) " *
				"must be a Symbol; got $(typeof(v)).",
			)
			ksym = v

			first_el = first(vec_arr)
			first_el isa NamedTuple ||
				Base.error(
					"select_field($(S), Val($(dim))) expects NamedTuple elements; " *
					"got $(typeof(first_el)).",
				)

			haskey(first_el, ksym) || Base.error(
				"NamedTuple elements on axis $(dim) for $(S) have no key $(ksym). " *
				"Available keys: $(collect(keys(first_el))).",
			)

			return [el[ksym] for el in vec_arr]
		else
			# No leaf field bound yet; just enforce NamedTuple contract and return as-is.
			first_el = first(vec_arr)
			first_el isa NamedTuple ||
				Base.error(
					"select_field($(S), Val($(dim))) is defined but resolved input has no " *
					"field $(kfield); data elements on axis $(dim) for $(S) must be " *
					"NamedTuple; got $(typeof(first_el)).",
				)
			return collect(vec_arr)
		end
	end
end


# Process one axis if present
@inline function axis_data(
	::Type{S},
	dim::Symbol,
	nt::NamedTuple,
	axis::Union{Axis, Nothing},
) where {S <: AbstractPlotSpec}
	axis === nothing && return nothing

	# axis selector: nt.x / nt.y / nt.z
	selector = getfield(nt, dim)::Symbol

	# 1) slice + select_field unwrapping (no scaling)
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
	build_series(::Type{S}, nt, axes) where {S<:AbstractPlotSpec}

Builds the vector of Dataseries for the given spec and resolved input `nt`.

Defaults to a single Dataseries corresponding to the primary primitive of
this spec, using the axis data computed by `axis_data`. Specs that need
multiple traces (overlays, histogram + CDF, etc.) should override this
method and typically still call `axis_data` under the hood.
"""
function build_series(
	::Type{S},
	nt::NamedTuple,
	axes::NamedTuple,
) where {S <: AbstractPlotSpec}
	dims = geom_axes(S)
	dims = dims isa Tuple ? dims : (dims,)

	xaxis = axes.xaxis
	yaxis = axes.yaxis
	zaxis = axes.zaxis

	xdata = :x in dims ? axis_data(S, :x, nt, xaxis) : nothing
	ydata = :y in dims ? axis_data(S, :y, nt, yaxis) : nothing
	zdata = :z in dims ? axis_data(S, :z, nt, zaxis) : nothing

	kind   = plot_kind(S)
	labels = legend_labels(S, nt)
	label  = isempty(labels) ? nothing : first(labels)

	series = Dataseries(
		kind,
		xdata,
		ydata,
		zdata,
		label,
	)

	return Dataseries[series]
end