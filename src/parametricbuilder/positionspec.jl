struct PositionSpec <: AbstractPositionSpec
	x0::Real
	y0::Real
	dx::Any
	dy::Any
	conn::Dict{String, Int}
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase mapping helpers
# Accepts:
#   :core => 1
#   ("core", 1)
#   (:core, 1)
#   [ :core => 1, :sheath => 0 ]
# etc.
# ─────────────────────────────────────────────────────────────────────────────
const _PhaseMapInputs = Union{
	Tuple{Symbol, Any},
	Tuple{String, Any},
	Pair{Symbol, Any},
	Pair{String, Any},
}

# normalize phases input to a splattable tuple of _PhaseMapInputs
_normalize_phase_map(p::_PhaseMapInputs) = (p,)
_normalize_phase_map(p::Tuple) = p
_normalize_phase_map(v::AbstractVector) = Tuple(v)
_normalize_phase_map(::Nothing) = ()

"""
	make_phase_maps(phases, n::Int)

Unified helper to process phase DSL inputs.
- If `n=1`, returns a vector with one Dict (used by `at`).
- If `n>1`, distributes values:
  - Scalars (e.g. `1`) are broadcast to all `n` legs.
  - Tuples/Vectors (e.g. `(1,2,3)`) are distributed to respective legs.
"""
function make_phase_maps(phases, n::Int)
	items = _normalize_phase_map(phases)
	out = [Dict{String, Int}() for _ in 1:n]

	for item in items
		# Extract key/value
		key_raw, val_raw = item isa Pair ? (first(item), last(item)) : (item[1], item[2])
		key = string(key_raw)

		# Distribute
		if val_raw isa Integer
			# Scalar broadcast
			v = Int(val_raw)
			for i in 1:n
				out[i][key] = v
			end
		elseif (val_raw isa Tuple || val_raw isa AbstractVector)
			# Vector distribution
			if length(val_raw) != n
				error(
					"Dimension mismatch for phase '$key': expected $n elements, got $(length(val_raw))",
				)
			end
			for i in 1:n
				out[i][key] = Int(val_raw[i])
			end
		else
			error(
				"Invalid phase value for '$key': expected Integer or collection of length $n, got $(typeof(val_raw))",
			)
		end
	end
	return out
end

function at(; x, y, dx = 0.0, dy = 0.0, phases = nothing)
	# n=1 for single position
	maps = make_phase_maps(phases, 1)
	return PositionSpec(x, y, dx, dy, maps[1])
end
