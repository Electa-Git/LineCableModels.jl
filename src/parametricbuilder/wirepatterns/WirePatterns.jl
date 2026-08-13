module WirePatterns

# ────────────────────────────────────────────────────────────────────────────
# Public API
# ────────────────────────────────────────────────────────────────────────────

# export ScreenPattern, HexaPattern
export make_stranded, make_screened

# ────────────────────────────────────────────────────────────────────────────
# Types
# ────────────────────────────────────────────────────────────────────────────

"""
	struct HexaPattern

Result for a single design choice.

Fields:
- `layers::Int`               — number of concentric layers (1 = center only).
- `wires::Int`                — total number of wires, N(L) = 1 + 3L(L-1).
- `wire_diameter_m::Float64`  — strand diameter [m].
- `total_area_m2::Float64`    — summed metallic area [m²].
- `awg::String`               — AWG label from the table (informative).
"""
struct HexaPattern
	layers::Int
	wires::Int
	wire_diameter_m::Float64
	total_area_m2::Float64
	awg::String
end

"""
	struct ScreenPattern

Screen wires design.

Fields:
- `wires::Int`               — number of wires on the wire array (N).
- `wire_diameter_m::Float64` — strand diameter [m].
- `lay_diameter_m::Float64`  — laying diameter Dm [m].
- `radius_m::Float64`        — wire array centerline radius = (Dm + d)/2 [m].
- `total_area_m2::Float64`   — N * π/4 * d^2 [m²].
- `coverage_pct::Float64`    — 100 * N*d / (π*Dm*sinα) [%].
- `awg::String`              — AWG label from the table (informative).
"""
struct ScreenPattern
	wires::Int
	wire_diameter_m::Float64
	lay_diameter_m::Float64
	radius_m::Float64
	total_area_m2::Float64
	coverage_pct::Float64
	awg::String
end

# ────────────────────────────────────────────────────────────────────────────
# Utils
# ────────────────────────────────────────────────────────────────────────────

_wire_area(dw::Real) = (pi/4) * (dw^2)  # area of one wire

# ---- AWG exact formulas (solid wire) ----
const _AWG_BASE  = 92.0
const _D0_MM     = 0.127            # 0.005 in in mm
const _AREA0_MM2 = 0.012668         # (π/4)*0.127^2
const _LN_BASE   = log(_AWG_BASE)

awg_to_d_mm(n::Real) = _D0_MM * (_AWG_BASE ^ ((36 - n)/39))
awg_to_area_mm2(n::Real) = _AREA0_MM2 * (_AWG_BASE ^ ((36 - n)/19.5))

d_mm_to_awg(d_mm::Real) = 36 - 39 * (log(d_mm/_D0_MM) / _LN_BASE)
area_mm2_to_awg(A_mm2::Real) = 36 - 19.5 * (log(A_mm2/_AREA0_MM2) / _LN_BASE)

function awg_label(n::Integer)
	n == -3 && return "0000 (4/0)"
	n == -2 && return "000 (3/0)"
	n == -1 && return "00 (2/0)"
	n == 0 && return "0 (1/0)"
	return string(n)
end

"Generate (label, diameter_m) for AWG n in [nmin, nmax]."
function awg_sizes(nmin::Integer = -3, nmax::Integer = 40)
	out = Tuple{String, Float64}[]
	@inbounds for n in nmin:nmax
		d_m = awg_to_d_mm(n) / 1000.0
		push!(out, (awg_label(n), d_m))
	end
	return out
end

"Apply a compaction/fill factor to solid area to approximate stranded metallic CSA."
stranded_area_mm2(n::Real; fill_factor::Real = 0.94) = fill_factor * awg_to_area_mm2(n)

# ────────────────────────────────────────────────────────────────────────────
# Hexagonal strand patterns
# ────────────────────────────────────────────────────────────────────────────

# ---- wire-count constraints per target area (mm²) ----
const _WIRE_RULES = Tuple{Int, Int, Union{Int, Nothing}}[
	(10, 6, 7),
	(16, 6, 7),
	(25, 6, 7),
	(35, 6, 7),
	(50, 6, 19),
	(70, 12, 19),
	(95, 15, 19),
	(120, 15, 37),
	(150, 15, 37),
	(185, 30, 37),
	(240, 30, 37),
	(300, 30, 61),
	(400, 53, 61),
	(500, 53, 61),
	(630, 53, 91),
	(800, 53, 91),
	(1000, 53, 91),
]

"""
	make_stranded(target_area_m2::Real; nmin::Integer=-3, nmax::Integer=40)

Compute hexagonal-pattern strand layouts that approximate or meet the target metallic cross-section, imposing allowed total-wire ranges by target area.

Inputs:
- `target_mm2`     — target metallic area [mm²].
- `nmin`,`nmax`    — AWG range to consider (default 4/0 … 40).

Returns:
- `best_match` — within allowed N(L), minimize |A − target|.
- `min_layers` — within allowed N(L) and A ≥ target, minimize layers (tie: smallest excess, then smaller diameter).
				 Fallback: within allowed, pick largest A < target (tie: smaller L, then smaller diameter).
- `min_diam`   — within allowed N(L) and A ≥ target, minimize diameter, then layers, then excess.
				 Fallback: within allowed, pick smallest diameter with largest A < target (then smallest L).
"""
function make_stranded(target_mm2::Real; nmin::Integer = -3, nmax::Integer = 40)
	@assert target_mm2 > 0 "Target cross-section must be positive."
	@assert nmin <= nmax "nmin must be ≤ nmax."

	target_area_m2 = target_mm2 * 1e-6  # m²
	# ---- hex geometry ----
	_hex_N(L::Int) = 1 + 3L*(L - 1)         # total wires after L layers
	_to_choice((dw, L, N, A, awg)) = HexaPattern(L, N, dw, A, awg)

	# Return (minN, maxN::Union{Int,Nothing}) for target in mm²
	function _allowed_wires(target_mm2::Real)
		for (thr, minN, maxN) in _WIRE_RULES
			if target_mm2 <= thr
				return (minN, maxN)
			end
		end
		return (53, nothing)  # > 1000 mm² -> min 53, no maximum
	end

	@inline function _allowed_N(N::Int, minN::Int, maxN::Union{Int, Nothing})
		maxN === nothing ? (N >= minN) : (N >= minN && N <= maxN)
	end

	# Allowed wire-count range from target (mm²)
	minN, maxN = _allowed_wires(target_mm2)

	# AWG sizes (label, d_m)
	sizes = awg_sizes(nmin, nmax)
	@assert !isempty(sizes) "AWG range produced no sizes."

	# Build allowed candidates: (dw, L, N, A, awg)
	candidates = Vector{Tuple{Float64, Int, Int, Float64, String}}()
	for (awg, dw) in sizes
		a1 = _wire_area(dw)
		@inbounds for L in 1:300
			N = _hex_N(L)
			if _allowed_N(N, minN, maxN)
				A = N * a1
				push!(candidates, (dw, L, N, A, awg))
			end
			if maxN !== nothing && N > maxN
				break
			end
		end
	end
	@assert !isempty(candidates) "No allowed candidates under the imposed wire-count span."

	# ---- best_match: minimize |A - target| (tie: smaller dw, then smaller L) ----
	rank_keys = [(abs(A - target_area_m2), dw, L) for (dw, L, N, A, _) in candidates]
	best_match = _to_choice(candidates[argmin(rank_keys)])
	if best_match.wires > 271
		@warn "Best match stranded pattern with a very high wire count ($(best_match.wires) wires). Consider revising the target cross-section or choosing a different wire configuration."
	end

	# Split feasible/infeasible for next selectors
	feas   = filter(((dw, L, N, A, awg),)->A >= target_area_m2, candidates)
	infeas = filter(((dw, L, N, A, awg),)->A < target_area_m2, candidates)

	# ---- min_layers ----
	if !isempty(feas)
		# minimal layers, then minimal excess, then smaller diameter
		keys_L = [(L, A - target_area_m2, dw) for (dw, L, N, A, _) in feas]
		min_layers = _to_choice(feas[argmin(keys_L)])
	else
		# fallback: closest from below (largest A), then minimal L, then smaller dw
		keys_fb = [(-A, L, dw) for (dw, L, N, A, _) in infeas]
		min_layers = _to_choice(infeas[argmin(keys_fb)])
	end

	# ---- min_diam ----
	if !isempty(feas)
		# smallest diameter; for it, minimal layers; then smallest excess
		sort!(feas, by = x -> (x[1], x[2], x[4] - target_area_m2))  # (dw asc, L asc, excess asc)
		min_diam = _to_choice(first(feas))
	else
		# fallback: smallest diameter with best undershoot; then minimal layers
		sort!(infeas, by = x -> (x[1], -(x[4]), x[2]))             # (dw asc, A desc, L asc)
		min_diam = _to_choice(first(infeas))
	end

	return (; best_match, min_layers, min_diam)
end

# ────────────────────────────────────────────────────────────────────────────
# Screen (single wire array) patterns
# ────────────────────────────────────────────────────────────────────────────

"""
	make_screened(A_req_m2::Real, Dm_m::Real;
				alpha_deg::Real=15.0, coverage_min_pct::Real=85.0,
				gap_frac::Real=0.0, min_wires::Int=3, extra_span::Int=8,
				nmin::Integer=-3, nmax::Integer=40)

Compute screen wire layouts that approximate or meet the target metallic cross-section, imposing:

  1) CSA:   N * (π/4) * d^2 ≥ A_req_m2
  2) Cover: (N*d)/(π*Dm*sinα) * 100 ≥ coverage_min_pct

while enforcing no-overlap wire array geometry (with optional clearance `gap_frac`).

Arguments:
- `A_req_mm2`         — required metallic cross-section [mm²].
- `Dm_mm`             — laying diameter (screen centerline) [mm].
- `alpha_deg`         — lay angle α in degrees (default 20°).
- `coverage_min_pct`  — required geometric coverage (default 85%).
- `gap_frac`          — extra clearance fraction in the no-overlap check (default 0).
- `min_wires`         — lower bound on N to avoid degenerate “non-wire-array" cases (default 3; set 6 for stronger symmetry).
- `extra_span`        — consider up to this many extra wires above the minimal requirement for better best_match search.
- `nmin`,`nmax`       — AWG range to consider (default 4/0 … 40).

Returns:
- `min_wires` — minimal N (≥ min_wires) that satisfies both CSA & coverage & geometry;
				tie-break: smaller d, then smaller excess area.
- `min_diam`  — smallest d that can satisfy both constraints; for it, minimal feasible N;
				tie-break: smaller excess area.
- `best_match`— among all feasible combos, area closest to A_req_m2; tie: smaller N, then smaller d.
"""
function make_screened(A_req_mm2::Real, Dm_mm::Real;
	alpha_deg::Real = 15.0, coverage_min_pct::Real = 85.0,
	gap_frac::Real = 0.0, min_wires::Int = 6, extra_span::Int = 8,
	nmin::Integer = -3, nmax::Integer = 40,
	coverage_max_pct::Real = 100.0,               # NEW: cap coverage to a single layer
	max_overshoot_pct::Real = 10.0,                # NEW: optional cap on A overshoot (∞ to disable)
	custom_diameters_mm::AbstractVector{<:Real} = Float64[])

	@assert 0.0 < coverage_min_pct <= 100.0
	@assert coverage_max_pct >= coverage_min_pct
	@assert max_overshoot_pct ≥ 0
	@assert A_req_mm2 > 0
	@assert Dm_mm > 0

	A_req_m2 = A_req_mm2 * 1e-6  # m²
	Dm_m     = Dm_mm * 1e-3      # m

	# --- helpers ---
	function _max_wires_single_layer(Dm::Real, d::Real; gap_frac::Real = 0.0)
		s = d*(1 + gap_frac) / (Dm + d)
		if !(0.0 < s < 1.0)
			;
			return 0;
		end
		return max(0, floor(Int, pi / asin(s)))
	end
	_to_choice((N, d, Dm, A, cov, awg)) = ScreenPattern(N, d, Dm, 0.5*(Dm + d), A, cov, awg)

	α  = deg2rad(alpha_deg)
	sα = sin(α);
	@assert sα > 0

	# AWG sizes + optional customs
	sizes = awg_sizes(nmin, nmax)
	for d in custom_diameters_mm
		push!(sizes, ("custom($(round(d; digits=3)) mm)", Float64(d)))
	end
	@assert !isempty(sizes)

	# Build candidates that satisfy BOTH constraints + geometry + coverage upper bound
	candidates = Tuple{Int, Float64, Float64, Float64, Float64, String}[]  # (N,d,Dm,A,cov,awg)

	for (awg, d) in sizes
		a1    = _wire_area(d)
		N_csa = ceil(Int, A_req_m2 / a1)
		N_cov = ceil(Int, (coverage_min_pct/100.0) * (pi*Dm_m*sα) / d)
		N_min = max(min_wires, N_csa, N_cov)
		N_max = _max_wires_single_layer(Dm_m, d; gap_frac = gap_frac)
		if N_max <= 0 || N_min > N_max
			continue
		end

		# Try N from N_min upward but reject coverage > coverage_max_pct and overshoot > max_overshoot_pct
		upper = min(N_min + extra_span, N_max)
		@inbounds for N in N_min:upper
			A   = N * a1
			cov = 100.0 * (N*d) / (pi*Dm_m*sα)
			if cov > coverage_max_pct
				break  # for fixed d, cov grows linearly with N; larger N will also violate
			end
			if isfinite(max_overshoot_pct)
				if A > A_req_m2 * (1 + max_overshoot_pct/100)
					continue
				end
			end
			push!(candidates, (N, d, Dm_m, A, cov, awg))
		end
	end

	@assert !isempty(candidates) "No feasible screen with given CSA, Dm, α, coverage bounds, and geometry."

	# --- selectors (tweaked) ---

	# 1) min_wires: minimize N; tie → minimize |A−Areq|; then smaller d
	keys_minN = [(N, abs(A - A_req_m2), d) for (N, d, _, A, _, _) in candidates]
	min_wires = _to_choice(candidates[argmin(keys_minN)])

	# 2) min_diam: smallest d; for it, minimal |A−Areq|; then minimal N
	sort!(candidates, by = x -> (x[2], abs(x[4] - A_req_m2), x[1]))  # (d asc, |ΔA| asc, N asc)
	min_diam = _to_choice(first(candidates))

	# 3) best_match: closest area to A_req; tie → smaller N, then smaller d
	keys_best = [(abs(A - A_req_m2), N, d) for (N, d, _, A, _, _) in candidates]
	best_match = _to_choice(candidates[argmin(keys_best)])

	return (; min_wires, min_diam, best_match)
end


end # module
