
# ─────────────────────────────────────────────────────────────────────────────
# Left-hand mapping syntax
# Accepts (:core,1) etc., with positions:
#   at(x=..., y=..., dx=..., dy=..., phases = (:core=>1, :sheath=>0, :jacket=>0))
# ─────────────────────────────────────────────────────────────────────────────
const _MapItem = Union{
	Tuple{Symbol, Int}, Tuple{String, Int},
	Pair{Symbol, Int}, Pair{String, Int},
}

_mapdict(items::_MapItem...) = Dict{String, Int}(
	(it isa Pair ? (string(first(it)) => Int(last(it)))
		: (string(it[1]) => Int(it[2])))
	for it in items
)

struct _Pos
	x0::Number
	y0::Number
	dx::Any
	dy::Any
	conn::Dict{String, Int}
end

# normalize phases input to a splattable tuple of _MapItem
_iter_phases(p) = (p,)                                  # single Pair or Tuple(:sym,Int)
_iter_phases(p::Tuple) = p                              # tuple of items
_iter_phases(v::AbstractVector{<:_MapItem}) = Tuple(v)  # vector of items

function at(; x, y, dx = 0.0, dy = 0.0, phases = nothing)
	items = phases === nothing ? () : _iter_phases(phases)
	return _Pos(x, y, dx, dy, _mapdict(items...))
end

# ─────────────────────────────────────────────────────────────────────────────
# Earth and system specs
# ─────────────────────────────────────────────────────────────────────────────
struct EarthSpec
	rho::Any;
	eps_r::Any;
	mu_r::Any;
	t::Any
end
EarthSpec(; rho, eps_r = 1.0, mu_r = 1.0, t = Inf) =
	EarthSpec(_spec(rho), _spec(eps_r), _spec(mu_r), _spec(t))

Earth(; rho, eps_r = 1.0, mu_r = 1.0, t = Inf) =
	EarthSpec(_spec(rho), _spec(eps_r), _spec(mu_r), _spec(t))

struct SystemBuilderSpec
	system_id::String
	builder::CableBuilderSpec
	positions::Vector{_Pos}
	length::Any         # (valuespec, pctspec) or scalar
	temperature::Any    # (valuespec, pctspec) or scalar
	earth::EarthSpec
	frequencies::Vector{Float64}
end

function SystemBuilderSpec(id::AbstractString, cbs::CableBuilderSpec,
	positions::Vector{_Pos};
	length = 1000.0, temperature = 20.0, earth::EarthSpec, f::AbstractVector{<:Real})
	return SystemBuilderSpec(
		String(id),
		cbs,
		positions,
		_spec(length),
		_spec(temperature),
		earth,
		collect(float.(f)),
	)
end

SystemBuilder(id::AbstractString, cbs::CableBuilderSpec,
	positions::Vector{_Pos};
	length = 1000.0, temperature = 20.0, earth::EarthSpec, f::AbstractVector{<:Real}) = SystemBuilderSpec(id, cbs, positions; length, temperature, earth, f)

# ─────────────────────────────────────────────────────────────────────────────
# Internals: expand range/% grammar via ParametricBuilder helpers
# ─────────────────────────────────────────────────────────────────────────────
@inline _expand_pair(specpair) = _make_range(specpair[1]; pct = specpair[2])

# (nothing, pct) on dx/dy ⇒ attach % to the anchor itself (no displacement sweep)
@inline function _axis(anchor::Number, dspec)
	spec, pct = _spec(dspec)
	if spec === nothing
		return _make_range(anchor; pct = pct)                # uncertain anchor
	else
		return (anchor .+ v for v in _make_range(spec; pct = pct))  # displaced anchor
	end
end

_expand_position(p::_Pos) =
	((x, y, p.conn) for x in _axis(p.x0, p.dx), y in _axis(p.y0, p.dy))

_expand_earth(e::EarthSpec) = (
	(ρ, ε, μ, t)
	for ρ in _expand_pair(e.rho), #_make_range(e.rho[1]; pct = e.rho[2]),
	ε in _expand_pair(e.eps_r), #_make_range(e.epsr[1]; pct = e.epsr[2]),
	μ in _expand_pair(e.mu_r), #_make_range(e.mur[1]; pct = e.mur[2]),
	t in _expand_pair(e.t) #_make_range(e.t[1]; pct = e.t[2])
)

# ─────────────────────────────────────────────────────────────────────────────
# Main iterator: yields fully-formed LineParametersProblem objects
# Overlaps are *not* emitted (skipped with warning  by catching the geometry error).
# Designs are identical per system realization (no cross-mixing).
# ─────────────────────────────────────────────────────────────────────────────
function iterate_problems(spec::SystemBuilderSpec)
	return Channel{LineParametersProblem}(32) do ch
		produced = 0
		try
			for des in spec.builder
				for L in _expand_pair(spec.length)
					pos_spaces = map(_expand_position, spec.positions)
					for choice in product(pos_spaces...)
						try
							x1, y1, c1 = choice[1]
							sys = DataModel.LineCableSystem(spec.system_id, L,
								DataModel.CablePosition(des, x1, y1, c1))
							for k in Iterators.drop(eachindex(choice), 1)
								xk, yk, ck = choice[k]
								sys = add!(sys, des, xk, yk, ck)
							end
							for T in _expand_pair(spec.temperature)
								for (ρ, ε, μ, t) in _expand_earth(spec.earth)
									em = EarthModel(spec.frequencies, ρ, ε, μ; t = t)
									prob = LineParametersProblem(sys;
										temperature = T, earth_props = em,
										frequencies = spec.frequencies)
									put!(ch, prob)
									produced += 1
								end
							end
						catch e
							if occursin("overlap", sprint(showerror, e))
								@warn sprint(showerror, e)
								@warn "Skipping..."
								continue
							else
								rethrow()
							end
						end
					end
				end
			end
		catch e
			@error "iterate SystemBuilderSpec failed" exception=(e, catch_backtrace())
		finally
			@debug "iterate SystemBuilderSpec finished" produced=produced upper_bound=cardinality(
				spec,
			)
		end
	end
end
