module UQ

# Export public API
export sample, mc, hist

# Module-specific dependencies
using ..Commons: BASE_FLOAT
import ..Commons: domain
using ..ParametricBuilder:
	MaterialSpec, PartSpec, CableBuilderSpec, SystemBuilderSpec, AbstractPositionSpec,
	PositionSpec, PositionGroupSpec, build, iterate, _spec, determinize
using ..Engine:
	EMTFormulation, compute!, LineParameters, PhaseDomain, ModalDomain,
	LineParamsDomain
using ..DataModel: get_outer_radius
using Measurements: Measurement, measurement, value, uncertainty
using Random, Statistics, DataFrames
using Distributions:
	Distributions, ContinuousUnivariateDistribution, Normal, Uniform, cdf, sampler
using StatsBase: fit, Histogram, normalize, quantile, ecdf
using LinearAlgebra


# Draw once from a "range-like" spec
#   spec :: Number                             → return as-is
#   spec :: AbstractVector                     → random element (uniform over indices)
#   spec :: (lo::Number, hi::Number, n::Int)   → given [lo, hi], interpret as ±1σ around μ = (lo+hi)/2, σ = (hi-lo)/2.
#   anything iterable                          → pick a random element
@inline function _rand_in(spec, distribution::Symbol)
	if spec isa Number
		return spec

	elseif spec isa AbstractVector
		@inbounds return spec[rand(1:length(spec))]

	elseif spec isa Tuple && length(spec) == 3 &&
		   spec[1] isa Number && spec[2] isa Number && spec[3] isa Integer
		lo, hi = spec[1], spec[2]
		# - Given [lo, hi], interpret as ±1σ around μ = (lo+hi)/2, σ = (hi-lo)/2.
		lo_f = float(lo)
		hi_f = float(hi)
		# TODO: handle edge case lo == hi when the nominal value is 0
		@assert hi_f > lo_f "hi must be greater than lo"
		μ = (lo_f + hi_f) / 2
		σ = (hi_f - lo_f) / 2

		if distribution === :normal
			# - :normal => Normal(μ, σ).
			return rand(Distributions.Normal(μ, σ))
		elseif distribution === :uniform
			# - :uniform  => Uniform(μ ± √3 σ) so std matches σ.
			d = √3 * σ
			return rand(Distributions.Uniform(μ - d, μ + d))
		else
			throw(
				ArgumentError(
					"unsupported distribution: $(distribution). Use :uniform or :normal",
				),
			)
		end

	else
		if Base.iterable(spec)
			vals = collect(spec)
			@inbounds return vals[rand(1:length(vals))]
		end
		return spec
	end
end

# Collapse a (spec, pct) pair → (value::Number, pct::Union{Nothing,Number})
"""
	_collapse_pair(sp::Tuple, distribution::Symbol; domain=nothing, max_tries::Int=10_000)

Collapse a (spec, pct) pair into `(value, pct_value)` by drawing once from the
"range-like" `spec` and `pct` using `_rand_in`.

If `domain !== nothing`, it must be a tuple `(lo, hi)` where each bound can be
`Real` or `nothing`. The value draw `v` is accepted only if:

	(lo === nothing || v ≥ lo) && (hi === nothing || v ≤ hi)

Otherwise a new draw is attempted, up to `max_tries`. If no feasible value is
found, an `error` is thrown.

This gives you generic rejection-sampling with minimal code, suitable for
enforcing physical domains like `(0, Inf)` for resistivity, spacing, thickness, etc.
"""
@inline function _collapse_pair(
	sp::Tuple,
	distribution::Symbol;
	domain::Union{Nothing, Tuple} = nothing,
	max_tries::Int = 10_000,
)
	spec, pct = sp

	# 1) VALUE SIDE (v): domain only for randomizable specs
	v =
		if domain === nothing || spec isa Number
			# No domain guardrail for scalars: if the user hard-codes nonsense,
			# let geometry/physics code blow up later.
			_rand_in(spec, distribution)
		else
			lo, hi = domain
			tries = 0
			accepted = nothing
			while true
				tries += 1
				tries > max_tries && error(
					"Unable to draw value in domain $domain from spec=$spec " *
					"after $max_tries attempts. Check your range and distribution.",
				)
				val = _rand_in(spec, distribution)
				if (lo === nothing || val >= lo) && (hi === nothing || val <= hi)
					accepted = val
					break
				end
			end
			accepted
		end

	# 2) PCT SIDE (u): domain only for randomizable pct-specs
	u =
		pct === nothing ? nothing :
		begin
			if pct isa Number
				# Scalar pct: pass through. If it's garbage, some other validator
				# or the physics will scream, not the domain sampler.
				_rand_in(pct, distribution)
			else
				# Range-like pct: enforce 0–100 on the *random draws*.
				lo, hi = 0.0, 100.0
				tries = 0
				accepted_pct = nothing
				while true
					tries += 1
					tries > max_tries && error(
						"Unable to draw pct ∈ [0,100] from pct-spec=$pct after $max_tries attempts.",
					)
					val = _rand_in(pct, distribution)
					if lo <= val <= hi
						accepted_pct = val
						break
					end
				end
				accepted_pct
			end
		end

	return (v, u)
end


# Collapse PartSpec.args:
#   each entry can be:
#     - scalar → keep as-is
#     - (spec, pct) → collapse to (rand_val, rand_pct)
# Treat an args entry as (spec, pct) only if first element is *not* Integer
@inline function _collapse_args(args::Tuple, distribution::Symbol)
	isempty(args) && return ()
	return tuple(
		(
			begin
				a = args[i]
				if (a isa Tuple) && (length(a) == 2) && !(a[1] isa Integer)
					_collapse_pair(a, distribution)
				elseif a isa AbstractVector
					@inbounds a[rand(1:length(a))]
				else
					a
				end
			end for i in eachindex(args)
		)...,
	)
end

# Collapse an entire MaterialSpec by collapsing each (spec, pct) field
@inline function _collapse_material(
	ms::MaterialSpec,
	distribution::Symbol,
)
	return MaterialSpec(;
		rho   = _collapse_pair(ms.rho, distribution),
		eps_r = _collapse_pair(ms.eps_r, distribution),
		mu_r  = _collapse_pair(ms.mu_r, distribution),
		T0    = _collapse_pair(ms.T0, distribution),
		alpha = _collapse_pair(ms.alpha, distribution),
	)
end

# Collapse one PartSpec → singleton PartSpec (no enumerations left)
@inline function _collapse_part(p::PartSpec, distribution::Symbol)
	new_dim  = _collapse_pair(p.dim, distribution)
	new_args = _collapse_args(p.args, distribution)
	new_mat  = _collapse_material(p.material, distribution)
	return PartSpec(p.component, p.part_type, p.n_layers;
		dim = new_dim, args = new_args, material = new_mat)
end

"""
	collapse(cbs::CableBuilderSpec; distribution::Symbol = :uniform) -> CableBuilderSpec

Return a **singleton** `CableBuilderSpec` by collapsing every range-like item
(dims, args, and material fields) into one random draw using the chosen distribution.

"""
function collapse(
	cbs::CableBuilderSpec;
	distribution::Symbol = :normal,
)
	parts = PartSpec[_collapse_part(p, distribution) for p in cbs.parts]
	return CableBuilderSpec(cbs.cable_id, parts, cbs.nominal)
end

"""
	sample(cbs::CableBuilderSpec; distribution::Symbol = :uniform) -> DataModel.CableDesign

Collapse ranges in `cbs` using `collapse` and build **one** cable design.
Useful for Monte Carlo style sampling where each call yields a new realization.
"""
function sample(
	cbs::CableBuilderSpec;
	distribution::Symbol = :normal,
)
	scbs = collapse(cbs; distribution = distribution)
	designs = build(scbs)              # with singleton choices, this yields length == 1
	@assert length(designs) == 1
	return designs[1]
end

"""
	_collapse_position(p::AbstractPositionSpec, distribution) -> PositionSpec or PositionGroupSpec

Collapse the uncertainty-bearing fields of a position specification.
No geometry is touched — grouped formations remain lazy, but their spacing
is collapsed to a concrete `(value, pct)` pair.
"""
# --- collapse for single positions -------------------------------------------------
function _collapse_position(p::PositionSpec, distribution::Symbol)
	dxc = _collapse_pair(_spec(p.dx), distribution)
	dyc = _collapse_pair(_spec(p.dy), distribution)
	return PositionSpec(
		p.x0,
		p.y0,
		dxc,
		dyc,
		p.conn,
	)
end
# --- collapse for grouped formations ----------------------------------------------
function _collapse_position(
	g::PositionGroupSpec,
	distribution::Symbol,
	d_min::Real;
	max_tries::Int = 10_000,
)
	# Physical constraint: d ≥ 2*R_out
	dspec_collapsed = _collapse_pair(
		g.d,
		distribution;
		domain = (d_min, nothing),   # (lo, hi), hi unconstrained
		max_tries = max_tries,
	)
	# value is guaranteed ≥ d_min here

	return PositionGroupSpec(
		g.arrangement,
		g.n,
		g.anchor,
		dspec_collapsed,
		g.conn,
	)
end

"""
	collapse(sbs::SystemBuilderSpec; distribution::Symbol = :uniform) -> SystemBuilderSpec

Collapse ranges in a `SystemBuilderSpec` using existing helpers.

Rules:
- Anchors `x, y` are numbers → pass through unchanged.
- `pos` are ((nom_range), (unc_range)) → `_collapse_position(..., distribution)`.
- `length`, `temperature` → `_collapse_pair(_spec(...), distribution)`.
- Earth fields (`rho`, `eps_r`, `mu_r`, `t`) → `_collapse_pair(_spec(...), distribution)`.
- Inner `builder` → `collapse(builder; distribution)`.
"""
function collapse(
	sbs::SystemBuilderSpec;
	distribution::Symbol = :normal,
)
	# 1) collapse cable builder (dims, mats, etc.)
	scbs = collapse(sbs.builder; distribution = distribution)

	# 2) build the *single* cable design and get its outer radius
	designs = build(scbs)
	@assert length(designs) == 1 "Collapsed CableBuilderSpec should yield exactly one design"
	des = designs[1]
	r_out = get_outer_radius(des)

	# 3) collapse positions: singles are collapsed generically,
	#    grouped formations are collapsed with geometry-aware rejection.
	pos = Vector{AbstractPositionSpec}(undef, length(sbs.positions))
	for (i, p) in enumerate(sbs.positions)
		if p isa PositionSpec
			pos[i] = _collapse_position(p, distribution)
		elseif p isa PositionGroupSpec
			pos[i] = _collapse_position(p, distribution, 2*r_out)
		else
			error("Unsupported position type in SystemBuilderSpec: $(typeof(p))")
		end
	end

	# 4) system-level scalars as before
	L = _collapse_pair(_spec(sbs.length), distribution)
	T = _collapse_pair(_spec(sbs.temperature), distribution)

	er = sbs.earth
	ρ = _collapse_pair(_spec(er.rho), distribution)
	ε = _collapse_pair(_spec(er.eps_r), distribution)
	μ = _collapse_pair(_spec(er.mu_r), distribution)
	t = _collapse_pair(_spec(er.t), distribution)
	earth = typeof(er)(; rho = ρ, eps_r = ε, mu_r = μ, t = t)

	return typeof(sbs)(
		sbs.system_id,
		scbs,
		pos;
		length      = L,
		temperature = T,
		earth       = earth,
		f           = sbs.frequencies,
	)
end


"""
	sample(sbs::SystemBuilderSpec; distribution::Symbol = :uniform)

Collapse ranges in `sbs` and produce one `LineParametersProblem`.
"""
function sample(
	sbs::SystemBuilderSpec;
	distribution::Symbol = :normal,
)
	ss = collapse(sbs; distribution = distribution)
	ch = iterate(ss)
	return take!(ch)
end

include("types.jl")
include("distributions.jl")
include("montecarlo.jl")
include("dataframe.jl")
include("plot.jl")

include("plotspecs/mcstatsplotspec.jl")

end # module UQ