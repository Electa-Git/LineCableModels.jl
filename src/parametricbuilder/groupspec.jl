# ─────────────────────────────────────────────────────────────────────────────
# Grouped formations: PositionGroupSpec
#
# These are *lazy* group specs. They do NOT carry concrete coordinates; they
# carry:
#   - an arrangement symbol (:trifoil, :hflat, :vflat, …)
#   - the number of legs n
#   - an anchor (x0,y0)
#   - a spacing spec (values,pct) via the same grammar as everything else
#   - per-leg connection maps (Dict{String,Int})
#
# They are materialized to concrete (x,y,conn) tuples *after* the CableDesign is
# known, so we can enforce a min spacing of 2 * outer_radius.
# ─────────────────────────────────────────────────────────────────────────────
struct PositionGroupSpec <: AbstractPositionSpec
	arrangement::Symbol                 # :trifoil, :hflat, :vflat, …
	n::Int                              # number of cables in the group
	anchor::Tuple{Float64, Float64}     # (x0,y0)
	d::Tuple{Any, Any}              # (valuespec, pctspec)
	conn::Vector{Dict{String, Int}}    # per-leg connection maps
end

# ─────────────────────────────────────────────────────────────────────────────
# Public sugar constructors
# ─────────────────────────────────────────────────────────────────────────────

"""
	trifoil(; x0 = 0.0, y0, d, phases)

Lazily describes a 3-cable trifoil formation. The anchor `(x0,y0)` is passed to
`trifoil_formation(x0,y0,d)` when the group is materialized.

The spacing `d` follows the usual `(valuespec, pctspec)` grammar; it will be
expanded lazily and clamped at runtime to avoid overlaps.
"""
function trifoil(; x0::Real = 0.0, y0::Real, d, phases)
	conn = make_phase_maps(phases, 3)
	return PositionGroupSpec(
		:trifoil,
		3,
		(float(x0), float(y0)),
		_spec(d),
		conn,
	)
end

"""
	hflat(; x0 = 0.0, y0 = 0.0, d, n = 3, phases)

Horizontal flat formation: first cable at `(x0, y0)`, remaining `n-1` cables at
`(x0 + k*d, y0)` for `k = 1, …, n-1`.

`d` accepts the `(valuespec, pctspec)` grammar.
"""
function hflat(; x0::Real = 0.0, y0::Real = 0.0, d, n::Integer = 3, phases)
	n < 1 && error("hflat requires n ≥ 1")
	conn = make_phase_maps(phases, n)
	return PositionGroupSpec(
		:hflat,
		n,
		(float(x0), float(y0)),
		_spec(d),
		conn,
	)
end

"""
	vflat(; x0 = 0.0, y0 = 0.0, d, n = 3, phases)

Vertical flat formation: first cable at `(x0, y0)`, remaining `n-1` cables at
`(x0, y0 - k*d)` for `k = 1, …, n-1`.

`d` accepts the `(valuespec, pctspec)` grammar.
"""
function vflat(; x0::Real = 0.0, y0::Real = 0.0, d, n::Integer = 3, phases)
	n < 1 && error("vflat requires n ≥ 1")
	conn = make_phase_maps(phases, n)
	return PositionGroupSpec(
		:vflat,
		n,
		(float(x0), float(y0)),
		_spec(d),
		conn,
	)
end

# ─────────────────────────────────────────────────────────────────────────────
# Group materialization (PositionGroupSpec)
# ─────────────────────────────────────────────────────────────────────────────

# Expand spacing spec and clamp out overlapping choices, based on radius
function _get_valid_spacings(g::PositionGroupSpec, rout)
	min_spacing = to_nominal(rout)

	# Full grid of values × pct → Measurement or plain Real
	raw = collect(_make_range(g.d[1]; pct = g.d[2]))

	# Nothing at all? auto-min with same pct grammar.
	if isempty(raw)
		return collect(_make_range(min_spacing; pct = g.d[2]))
	end

	valid     = Any[]
	discarded = 0

	# Filter by geometry, but KEEP the original object (Measurement or Real)
	for s in raw
		ds = to_nominal(s)
		if ds >= min_spacing
			push!(valid, s)          # don't strip uncertainty
		else
			discarded += 1
		end
	end

	# CASE 1: all invalid → pure AUTO: min_spacing with all pcts
	if isempty(valid)
		@debug "Spacing spec produced only overlapping layouts; clamping to minimum with % uncertainty grid." min_spacing=min_spacing
		return collect(_make_range(min_spacing; pct = g.d[2]))
	end

	# CASE 2: some valid, some discarded → inject ONE batch at min_spacing,
	# but only for spacing+uncertainty combos that are not already present.
	if discarded > 0
		@debug "Dropped $discarded spacing samples below minimum center-to-center distance; including one batch at the minimum feasible spacing." min_spacing=min_spacing

		autos_all = collect(_make_range(min_spacing; pct = g.d[2]))

		# Use a Set to avoid injecting exact duplicates (same Measurement).
		valid_set = Set(valid)
		autos = Any[]
		for a in autos_all
			if !(a in valid_set)
				push!(autos, a)
			end
		end

		# Prepend autos so min_spacing layouts come first, but WITHOUT
		# multiplying cardinality by cloning identical points.
		valid = vcat(autos, valid)
	end

	return valid
end


"""
	_materialize(g::PositionGroupSpec, des::CableDesign)

Lazily expands a grouped formation into concrete `(x, y, conn)` blocks after the
external radius is known (`des` is the fully materialized design).

Returns a generator of `Vector{Tuple{Float64,Float64,Dict{String,Int}}}`, one
vector per valid spacing choice.
"""
function _materialize(g::PositionGroupSpec, des::CableDesign)
	r = get_outer_radius(des)
	spacings = _get_valid_spacings(g, r)

	# Build one concrete layout (vector of (x,y,conn)) for a given spacing d
	function _make_layout(g::PositionGroupSpec, d::Real)
		x0, y0 = g.anchor

		coords =
			g.arrangement == :trifoil ?
			begin
				g.n == 3 || error("trifoil formation expects n = 3, got $(g.n)")
				x0p, y0p, dp = promote(x0, y0, d)
				xa, ya, xb, yb, xc, yc = DataModel.trifoil_formation(x0p, y0p, dp)
				[(xa, ya), (xb, yb), (xc, yc)]
			end :
			g.arrangement == :hflat ? begin
				[(x0 + d * (i - 1), y0) for i in 1:g.n]
			end :
			g.arrangement == :vflat ? begin
				[(x0, y0 - d * (i - 1)) for i in 1:g.n]
			end :
			error("Unknown position group arrangement $(g.arrangement)")

		return [(x, y, g.conn[i]) for (i, (x, y)) in enumerate(coords)]
	end

	return (_make_layout(g, d) for d in spacings)
end


"""
	_expand_position(position_defs, des)

Top-level helper that yields flattened `Vector{(x,y,conn)}` for every allowed
combination of positions/groups.
"""
function _expand_position(position_defs::Vector{AbstractPositionSpec}, des::CableDesign)
	spaces = Vector{Any}(undef, length(position_defs))

	for (i, p) in pairs(position_defs)
		if p isa PositionGroupSpec
			# group: generator of Vector{(x,y,conn)}
			spaces[i] = _materialize(p, des)
		elseif p isa PositionSpec
			# single: wrap each (x,y,conn) into a 1-element vector so the
			# outer logic can always `vcat` vectors.
			spaces[i] = (
				[(x, y, p.conn)]
				for x in _axis(p.x0, p.dx),
				y in _axis(p.y0, p.dy)
			)
		else
			error("Unsupported position spec type: $(typeof(p))")
		end
	end

	return (reduce(vcat, combo) for combo in product(spaces...))
end
