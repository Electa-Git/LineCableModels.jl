"""
	build_figures(::Type{S}, nt, panels) where {S<:AbstractPlotSpec}

Packs Panel values into Figure payloads.

Default behavior:
- a single Figure is created for this spec call,
- `layout` is chosen based on `figure_layout(S)` trait.

Specs that want multiple OS windows or more complex layout policies may
override this method.
"""
function build_figures(
	::Type{S},
	nt::NamedTuple,
	panels::Vector{Panel},
) where {S <: AbstractPlotSpec}
	layout = figure_layout(S)                # :windows, :grid, ...

	# Do not materialize plots if length == 1.
	# Upstream grouping/materialization stays intact, but the final Figure[] is empty.
	if !isempty(panels)
		p1 = first(panels)
		if !isempty(p1.series)
			ds1 = first(p1.series)

			data1 =
				ds1.xdata === nothing ? (ds1.ydata === nothing ? ds1.zdata : ds1.ydata) :
				ds1.xdata
			if data1 !== nothing
				n = length(data1)
				if n <= 1
					@warn "Skipping plot for spec $(S): sample length is $(n) (need ≥ 2)."
					return Figure[]
				end
			end
		end
	end

	fig_kwargs = nt.renderer
	figsize = default_figsize(S)
	fig_title = default_title(S, nt)

	fig = Figure(fig_title, figsize, layout, panels, fig_kwargs)
	return Figure[fig]
end

# Determine if, for spec S and resolved input nt, the leaf seen by axis_slice /
# build_series behaves as "scalar" (numeric/complex) or as a NamedTuple that
# should be exploded at the figure / panel level.
function is_leaf(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec}
	dims_raw = geom_axes(S)
	dims = dims_raw isa Tuple ? dims_raw : (dims_raw,)

	for dim in dims
		kfield = select_field(S, Val(dim))
		kfield === nothing && continue

		if haskey(nt, kfield)
			# A concrete field Symbol is pinned in the resolved input; axis_slice
			# will unwrap the NamedTuple and return numeric data.
			return true
		else
			# Spec uses select_field on this axis but no field was chosen yet;
			# axis_slice will return a vector of NamedTuples.
			return false
		end
	end

	# No axis uses select_field at all → grammar never unwraps fields.
	return true
end

# Decide high-level figure mode given spec S and resolved input nt.
function resolve_group_mode(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec}

	trait_mode = grouping_mode(S)  #  Symbol
	if trait_mode !== :auto
		return trait_mode
	end

	idx_raw = index_keys(S)
	idx = idx_raw isa Tuple ? idx_raw : (idx_raw,)

	# Matrix-like selector indices (pure selectors, never ranged)
	has_i = :i in idx
	has_j = :j in idx
	has_matrix = has_i || has_j

	i_defined = has_i && haskey(nt, :i)
	j_defined = has_j && haskey(nt, :j)


	if is_leaf(S, nt)
		if !has_matrix
			return :single
		end

		all_pinned = (!has_i || i_defined) && (!has_j || j_defined)
		if all_pinned
			return :single
		else
			return :overlay_ij
		end
	else
		if !has_matrix
			return :overlay_fields
		end

		all_pinned = (!has_i || i_defined) && (!has_j || j_defined)
		if all_pinned
			return :overlay_fields
		else
			return :per_ij_overlay_fields
		end
	end
end

# Infer matrix dimensions (Ni, Nj) from any axis container that carries the
# i/j selector structure. Falls back to (1,1) if no matrix indices are used.
function matrix_size(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec}
	idx_raw = index_keys(S)
	idx = idx_raw isa Tuple ? idx_raw : (idx_raw,)

	has_i = :i in idx
	has_j = :j in idx

	Ni = 1
	Nj = 1

	(!has_i && !has_j) && return Ni, Nj

	dims_raw = geom_axes(S)
	dims = dims_raw isa Tuple ? dims_raw : (dims_raw,)
	obj = nt.obj

	for dim in dims
		# skip if this axis has no selector in nt (e.g. z unused)
		!haskey(nt, dim) && continue

		datakey = getfield(nt, dim)
		datakey isa Symbol || continue

		arr = container_array(S, obj, dim, datakey)

		arr isa AbstractArray || continue

		if has_i && size(arr, 1) > 1
			Ni = size(arr, 1)
		end
		if has_j && ndims(arr) >= 2 && size(arr, 2) > 1
			Nj = size(arr, 2)
		end

		if (!has_i || Ni > 1) && (!has_j || Nj > 1)
			break
		end
	end

	return Ni, Nj
end

# For specs where target data is a NamedTuple (select_field used but no field
# chosen yet), infer which axis carries the NamedTuple leaf and what its
# field keys are, using axis_slice to get a vector of NamedTuples.
function get_target_fields(
	::Type{S},
	nt::NamedTuple,
	axes::NamedTuple,
) where {S <: AbstractPlotSpec}
	dims_raw = geom_axes(S)
	dims = dims_raw isa Tuple ? dims_raw : (dims_raw,)

	# Find the first axis that uses select_field
	dim_field = nothing
	kfield    = nothing
	for dim in dims
		kf = select_field(S, Val(dim))
		kf === nothing && continue
		dim_field = dim
		kfield    = kf
		break
	end

	if dim_field === nothing || kfield === nothing
		Base.error("get_target_fields: no axis uses select_field for spec $(S).")
	end

	# axis descriptor
	axis =
		dim_field === :x ? axes.xaxis :
		dim_field === :y ? axes.yaxis :
		axes.zaxis

	axis === nothing &&
		Base.error(
			"get_target_fields: axis $(dim_field) has no Axis for spec $(S).",
		)

	# axis_slice will return a vector of NamedTuples in the "NamedTuple leaf" modes
	vec = axis_slice(S, nt, axis, Val(dim_field))
	isempty(vec) &&
		Base.error(
			"get_target_fields: empty data along axis $(dim_field) for spec $(S); cannot infer NamedTuple fields.",
		)

	leaf = first(vec)
	leaf isa NamedTuple ||
		Base.error(
			"get_target_fields: expected NamedTuple leaf for spec $(S) axis $(dim_field), got $(typeof(leaf)).",
		)

	field_keys = collect(keys(leaf))
	return dim_field, kfield, field_keys
end

"""
	is_modal(obj)

Return `true` iff `domain(obj)` is a modal-like tag.
"""
@inline is_modal(x) = (D = domain(x); D !== nothing && D <: ModalDomain)

# --------------------------------------------------------------------------
# Index pair iterator: PhaseDomain vs ModalDomain
# --------------------------------------------------------------------------

"""
	index_pairs(::Type{S}, nt) where {S<:AbstractPlotSpec}

Return the list of \\((i,j)\\) index pairs that this spec should materialize
for the given resolved input `nt`.

Semantics:

- Phase-like domain (default, or `domain(nt.obj) === nothing`):
	* If :i is in `index_keys(S)`:
		- If `nt` pins `i`, use only that value.
		- Otherwise, use `1:Ni` where `Ni` is inferred from `matrix_size(S, nt)`.
	* If :j is in `index_keys(S)`:
		- Same, using `nj` / `Nj`.

  Result: full rectangular coverage over the active ranges.

- ModalDomain AND both :i and :j are in `index_keys(S)` AND neither is pinned in `nt`:
	* Let `(Ni, Nj) = matrix_size(S, nt)` and `N = min(Ni, Nj)`.
	* Return only diagonal pairs: `(1,1), (2,2), ..., (N,N)`.

- If the user pins `i` and/or `j`, user intent takes precedence regardless
  of domain tag.
"""
function index_pairs(
	::Type{S},
	nt::NamedTuple,
) where {S <: AbstractPlotSpec}
	idx_raw = index_keys(S)
	idx = idx_raw isa Tuple ? idx_raw : (idx_raw,)

	has_i = :i in idx
	has_j = :j in idx

	Ni, Nj = matrix_size(S, nt)

	i_defined = has_i && haskey(nt, :i)
	j_defined = has_j && haskey(nt, :j)

	obj = nt.obj
	modal = is_modal(obj)

	# Modal diagonal semantics:
	# - both :i and :j are matrix selectors,
	# - neither is pinned by the user or defaults,
	# - domain is ModalDomain (or subtype).
	if modal && has_i && has_j && !i_defined && !j_defined
		N = min(Ni, Nj)
		pairs = Vector{Tuple{Int, Int}}(undef, N)
		@inbounds for k in 1:N
			pairs[k] = (k, k)
		end
		return pairs
	end

	# Phase-like / generic semantics (or user-pinned case in modal)
	I_range =
		if has_i
			i_defined ? (nt.i:nt.i) : (1:Ni)
		else
			1:1
		end

	J_range =
		if has_j
			j_defined ? (nt.j:nt.j) : (1:Nj)
		else
			1:1
		end

	pairs = Tuple{Int, Int}[]
	@inbounds for i in I_range
		for j in J_range
			push!(pairs, (i, j))
		end
	end

	return pairs
end


"""
	build_figures(::Type{S}, nt) where {S<:AbstractPlotSpec}

Top-level figure builder for spec `S`.

Decides a high-level figure mode based on:
- whether the leaf behaves as scalar or NamedTuple (via select_field traits),
- whether index selectors :i and :j are present in `index_keys(S)`,
- whether :i and :j are defined or free in the resolved input `nt`.

It then delegates to `build_figures(::Type{S}, ::Val{mode}, nt, axes)` where
`mode` is one of:
- :single                 → single atom, no i/j or field expansion
- :overlay_ij              → scalar leaf, some of :i/:j free → overlay all (i,j)
- :overlay_fields          → NamedTuple leaf, fixed (i,j) → overlay all fields
- :per_ij_overlay_fields → NamedTuple leaf, some of :i/:j free → one Panel
							  per (i,j), overlaying all fields in each panel.
"""
function build_figures(
	::Type{S},
	nt::NamedTuple,
) where {S <: AbstractPlotSpec}
	axes = build_axes(S, nt)
	mode = resolve_group_mode(S, nt)
	return build_figures(S, Val(mode), nt, axes)
end

function build_figures(
	::Type{S},
	::Val{:single},
	nt::NamedTuple,
	axes::NamedTuple,
) where {S <: AbstractPlotSpec}
	series = build_series(S, nt, axes)
	panels = build_panels(S, nt, axes, series)
	return build_figures(S, nt, panels)
end

function build_figures(
	::Type{S},
	::Val{:overlay_ij},
	nt::NamedTuple,
	axes::NamedTuple,
) where {S <: AbstractPlotSpec}
	idx_raw = index_keys(S)
	idx = idx_raw isa Tuple ? idx_raw : (idx_raw,)

	has_i = :i in idx
	has_j = :j in idx

	all_series = Dataseries[]

	for (i, j) in index_pairs(S, nt)
		nt_ij = nt
		if has_i
			nt_ij = merge(nt_ij, (; i = i))
		end
		if has_j
			nt_ij = merge(nt_ij, (; j = j))
		end

		series_ij = build_series(S, nt_ij, axes)
		append!(all_series, series_ij)
	end

	panels = build_panels(S, nt, axes, all_series)
	return build_figures(S, nt, panels)
end

function build_figures(
	::Type{S},
	::Val{:overlay_fields},
	nt::NamedTuple,
	axes::NamedTuple,
) where {S <: AbstractPlotSpec}
	_, kfield, field_keys = get_target_fields(S, nt, axes)

	all_series = Dataseries[]

	for fk in field_keys
		nt_fk = merge(nt, (; kfield => fk))
		series_fk = build_series(S, nt_fk, axes)
		append!(all_series, series_fk)
	end

	panels = build_panels(S, nt, axes, all_series)
	return build_figures(S, nt, panels)
end

function build_figures(
	::Type{S},
	::Val{:per_ij_overlay_fields},
	nt::NamedTuple,
	axes::NamedTuple,
) where {S <: AbstractPlotSpec}
	idx_raw = index_keys(S)
	idx = idx_raw isa Tuple ? idx_raw : (idx_raw,)

	has_i = :i in idx
	has_j = :j in idx

	_, kfield, field_keys = get_target_fields(S, nt, axes)

	areas = Panel[]

	for (i, j) in index_pairs(S, nt)
		nt_ij = nt
		if has_i
			nt_ij = merge(nt_ij, (; i = i))
		end
		if has_j
			nt_ij = merge(nt_ij, (; j = j))
		end

		series_ij = Dataseries[]

		for fk in field_keys
			nt_ij_fk = merge(nt_ij, (; kfield => fk))
			series_fk = build_series(S, nt_ij_fk, axes)
			append!(series_ij, series_fk)
		end

		# Populate area key with whatever matrix indices this spec actually uses.
		key =
			if has_i && has_j
				(; i = i, j = j)
			elseif has_i
				(; i = i)
			elseif has_j
				(; j = j)
			else
				(;)
			end

		title = default_title(S, nt_ij)

		area = Panel(
			axes.xaxis,
			axes.yaxis,
			axes.zaxis,
			title,
			series_ij,
			key,
		)

		push!(areas, area)
	end

	return build_figures(S, nt, areas)
end


