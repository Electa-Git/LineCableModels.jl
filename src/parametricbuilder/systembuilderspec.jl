# Positions in the parametric system builder.
#
# Two flavours:
#   - PositionSpec      : single anchor with dx/dy ranges (arbitrary layouts)
#   - PositionGroupSpec : defined in `trifoil.jl`, grouped formations
#
# Both subtype AbstractPositionSpec so SystemBuilderSpec can store a mixed vector.
abstract type AbstractPositionSpec end

include("positionspec.jl")
include("groupspec.jl")

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
	positions::Vector{AbstractPositionSpec}
	length::Any         # (valuespec, pctspec) or scalar
	temperature::Any    # (valuespec, pctspec) or scalar
	earth::EarthSpec
	frequencies::Vector{Float64}
end

function SystemBuilderSpec(id::AbstractString, cbs::CableBuilderSpec,
	positions::Vector{<:AbstractPositionSpec};
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
	positions::AbstractVector{<:AbstractPositionSpec};
	length = 1000.0, temperature = 20.0, earth::EarthSpec, f::AbstractVector{<:Real}) = SystemBuilderSpec(id, cbs, positions; length, temperature, earth, f)

SystemBuilder(id::AbstractString, cbs::CableBuilderSpec,
	positions::AbstractPositionSpec;
	length = 1000.0, temperature = 20.0, earth::EarthSpec, f::AbstractVector{<:Real}) = SystemBuilderSpec(id, cbs, [positions]; length, temperature, earth, f)

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

_expand_earth(e::EarthSpec) = (
	(ρ, ε, μ, t)
	for ρ in _expand_pair(e.rho),
	ε in _expand_pair(e.eps_r),
	μ in _expand_pair(e.mu_r),
	t in _expand_pair(e.t)
)

# Choice count for single positions: size of the dx × dy grid
function _position_choice_count(p::PositionSpec)
	nx = length(collect(_axis(p.x0, p.dx)))
	ny = length(collect(_axis(p.y0, p.dy)))
	return nx * ny
end

# Choice count for grouped positions: number of spacing samples
function _position_choice_count(g::PositionGroupSpec)
	spec, pct = g.d
	return length(collect(_make_range(spec; pct = pct)))
end

# ─────────────────────────────────────────────────────────────────────────────
# Main iterator: yields fully-formed LineParametersProblem objects
# Overlaps are *not* emitted (skipped with warning  by catching the geometry error).
# Designs are identical per system realization (no cross-mixing).
# ─────────────────────────────────────────────────────────────────────────────
function iterate(spec::SystemBuilderSpec)
	return Channel{LineParametersProblem}(32) do ch
		produced = 0
		try
			for des in spec.builder
				for L in _expand_pair(spec.length)
					# NB: _expand_position keeps grouped spacings atomic and
					# materializes after `des` (and its outer radius) are known.
					for choice in _expand_position(spec.positions, des)
						try
							x1, y1, c1 = choice[1]
							sys = DataModel.LineCableSystem(
								spec.system_id,
								L,
								DataModel.CablePosition(des, x1, y1, c1),
							)

							for k in Iterators.drop(eachindex(choice), 1)
								xk, yk, ck = choice[k]
								sys = add!(sys, des, xk, yk, ck)
							end

							for T in _expand_pair(spec.temperature)
								for (ρ, ε, μ, t) in _expand_earth(spec.earth)
									em = EarthModel(spec.frequencies, ρ, ε, μ; t = t)
									prob = LineParametersProblem(
										sys;
										temperature = T,
										earth_props = em,
										frequencies = spec.frequencies,
									)
									put!(ch, prob)
									produced += 1
								end
							end
						catch e
							if occursin("overlap", sprint(showerror, e)) || occursin(
								"conductor resistivity must be positive",
								sprint(showerror, e),
							)
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
			@error "iterate SystemBuilderSpec failed" exception = (e, catch_backtrace())
		finally
			@debug "iterate SystemBuilderSpec finished" produced = produced upper_bound =
				cardinality(spec)
		end
	end
end

function build(spec::SystemBuilderSpec)
	problems = LineParametersProblem[]
	for prob in spec              # uses Base.iterate(spec::SystemBuilderSpec)
		push!(problems, prob)
	end
	return problems
end
