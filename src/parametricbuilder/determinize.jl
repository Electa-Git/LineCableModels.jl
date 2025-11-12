# ─────────────────────────────────────────────────────────────────────────────
# Deterministic collapse (Monte Carlo on deterministic ranges only)
# Policy: transform (valuespec, pctspec) → (merged_valuespec, nothing)
# ─────────────────────────────────────────────────────────────────────────────

# Percent helpers
@inline _pct(u) = float(u) / 100
@inline _expand_nom(nom::Number, u::Number) = (nom*(1 - _pct(u)), nom*(1 + _pct(u)))
@inline _expand_bounds(lo::Number, hi::Number, u1::Number, u2::Number) =
	(lo*(1 - _pct(u1)), hi*(1 + _pct(u2)))

# Deterministic collapse with pct interpreted as percent (not absolute)
@inline function _det_pair(spec, pct)
	pct === nothing && return (spec, nothing)

	# helper: largest percent magnitude in the tuple
	_umax(u1, u2) = max(abs(float(u1)), abs(float(u2)))

	# A) spec = (lo,hi,N1), pct = (u1,u2,N2)
	if (spec isa Tuple && length(spec)==3 && all(x->x isa Number, spec)) &&
	   (pct isa Tuple && length(pct) == 3 && all(x->x isa Number, pct))
		lo, hi, N1 = float(spec[1]), float(spec[2]), Int(spec[3])
		u1, u2, N2 = float(pct[1]), float(pct[2]), Int(pct[3])
		u = _umax(u1, u2)
		lo_det = lo * (1 - _pct(u))
		hi_det = hi * (1 + _pct(u))
		return ((lo_det, hi_det, N1 * N2), nothing)
	end

	# B) spec = (lo,hi,N1), pct = u
	if (spec isa Tuple && length(spec)==3 && all(x->x isa Number, spec)) && (pct isa Number)
		lo, hi, N1 = float(spec[1]), float(spec[2]), Int(spec[3])
		u = abs(float(pct))
		lo_det = lo * (1 - _pct(u))
		hi_det = hi * (1 + _pct(u))
		return ((lo_det, hi_det, N1), nothing)
	end

	# C) spec = nom, pct = (u1,u2,N2)
	if (spec isa Number) && (pct isa Tuple && length(pct)==3 && all(x->x isa Number, pct))
		nom = float(spec)
		u1, u2, N2 = float(pct[1]), float(pct[2]), Int(pct[3])
		u = _umax(u1, u2)
		lo_det = nom * (1 - _pct(u))
		hi_det = nom * (1 + _pct(u))
		return ((lo_det, hi_det, max(N2, 2)), nothing)
	end

	# D) spec = nom, pct = u
	if (spec isa Number) && (pct isa Number)
		nom = float(spec);
		u = abs(float(pct))
		lo_det = nom * (1 - _pct(u))
		hi_det = nom * (1 + _pct(u))
		return ((lo_det, hi_det, 2), nothing)
	end

	# E) fallback
	return (spec, nothing)
end

# Normalizer: accept a field already in (spec,pct) or as a scalar → return (spec’, nothing)
@inline _det_field(x) = (x isa Tuple && length(x)==2) ? _det_pair(x[1], x[2]) : (x, nothing)

# ---- MaterialSpec ----
function determinize(ms::MaterialSpec)
	MaterialSpec(
		rho   = _det_field(ms.rho),
		eps_r = _det_field(ms.eps_r),
		mu_r  = _det_field(ms.mu_r),
		T0    = _det_field(ms.T0),
		alpha = _det_field(ms.alpha),
	)
end

# ---- PartSpec (dim, args, material) ----
function determinize(ps::PartSpec)
	dim_det = _det_field(ps.dim)
	# each arg can be scalar or (spec,pct)
	args_det = map(a -> (a isa Tuple && length(a)==2) ? _det_field(a) : a, ps.args) |> Tuple
	mat_det  = determinize(ps.material)
	return PartSpec(
		ps.component,
		ps.part_type,
		ps.n_layers;
		dim = dim_det,
		args = args_det,
		material = mat_det,
	)
end

# ---- CableBuilderSpec (vector/nested parts) ----
function determinize(cbs::CableBuilderSpec)
	parts_det = PartSpec[determinize(p) for p in cbs.parts]
	return CableBuilderSpec(cbs.cable_id, parts_det, cbs.nominal)
end

# ─────────────────────────────────────────────────────────────────────────────
# Deterministic collapse for SystemBuilderSpec (non-materializing)
# ─────────────────────────────────────────────────────────────────────────────

@inline _det_axis(a) = (a isa Tuple && length(a)==2) ? _det_pair(a[1], a[2]) : a
# determinize EarthSpec
function determinize(e::EarthSpec)
	EarthSpec(
		rho   = _det_field(e.rho),
		eps_r = _det_field(e.eps_r),
		mu_r  = _det_field(e.mu_r),
		t     = _det_field(e.t),
	)
end

# determinize _Pos (keep anchors; just collapse dx/dy specs)
function determinize(p::_Pos)
	_Pos(p.x0, p.y0, _det_axis(p.dx), _det_axis(p.dy), p.conn)
end

# determinize SystemBuilderSpec
function determinize(s::SystemBuilderSpec)
	SystemBuilderSpec(
		s.system_id,
		determinize(s.builder),
		[determinize(p) for p in s.positions];
		length      = _det_field(s.length),
		temperature = _det_field(s.temperature),
		earth       = determinize(s.earth),
		f           = s.frequencies,
	)
end