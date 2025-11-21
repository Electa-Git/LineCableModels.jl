struct PositionSpec
	x0::Real
	y0::Real
	dx::Any
	dy::Any
	conn::Dict{String, Int}
end

# ─────────────────────────────────────────────────────────────────────────────
# Left-hand mapping syntax
# Accepts (:core,1) etc., with positions:
#   at(x=..., y=..., dx=..., dy=..., phases = (:core=>1, :sheath=>0, :jacket=>0))
# ─────────────────────────────────────────────────────────────────────────────
const _PhaseMapItem = Union{
	Tuple{Symbol, Int}, Tuple{String, Int},
	Pair{Symbol, Int}, Pair{String, Int},
}

_parse_phase_map(items::_PhaseMapItem...) = Dict{String, Int}(
	(it isa Pair ? (string(first(it)) => Int(last(it)))
		: (string(it[1]) => Int(it[2])))
	for it in items
)

# normalize phases input to a splattable tuple of _PhaseMapItem
_iter_phases(p) = (p,)                                  # single Pair or Tuple(:sym,Int)
_iter_phases(p::Tuple) = p                              # tuple of items
_iter_phases(v::AbstractVector{<:_PhaseMapItem}) = Tuple(v)  # vector of items

function at(; x, y, dx = 0.0, dy = 0.0, phases = nothing)
	items = phases === nothing ? () : _iter_phases(phases)
	return PositionSpec(x, y, dx, dy, _parse_phase_map(items...))
end