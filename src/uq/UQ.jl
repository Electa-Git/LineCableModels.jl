module UQ

# Export public API
export sample, mc

# Module-specific dependencies
using ..ParametricBuilder:
	MaterialSpec, PartSpec, CableBuilderSpec, SystemBuilderSpec, build, iterate_problems,
	_spec, determinize
using Measurements: Measurement, value, uncertainty
using Random, Statistics, DataFrames
using Distributions: Distributions, Normal, Uniform
using StatsBase: fit, Histogram

# ─────────────────────────────────────────────────────────────────────────────
# Random collapse of ranges → singleton specs + one-shot build
# ─────────────────────────────────────────────────────────────────────────────

# - Given [lo, hi], interpret as ±1σ around μ = (lo+hi)/2, σ = (hi-lo)/2.
# - :gaussian => Normal(μ, σ).
# - :uniform  => Uniform(μ ± √3 σ) so std matches σ.
# - Any other symbol => hard error.
@inline function _bounded_draw(lo::Real, hi::Real, distribution::Symbol = :uniform)
	lo_f = float(lo)
	hi_f = float(hi)
	@assert hi_f > lo_f "hi must be greater than lo"
	μ = (lo_f + hi_f) / 2
	σ = (hi_f - lo_f) / 2

	if distribution === :gaussian || distribution === :normal
		return rand(Distributions.Normal(μ, σ))
	elseif distribution === :uniform
		d = √3 * σ
		return rand(Distributions.Uniform(μ - d, μ + d))
	else
		throw(
			ArgumentError(
				"unsupported distribution: $(distribution). Use :uniform or :gaussian",
			),
		)
	end
end

# Draw once from a "range-like" spec
#   spec :: Number                             → return as-is
#   spec :: AbstractVector                     → random element (uniform over indices)
#   spec :: (lo::Number, hi::Number, n::Int)   → scalar draw in [lo, hi] (n only defines a grid elsewhere)
#   anything iterable                          → pick a random element
@inline function _rand_in(spec, distribution::Symbol)
	if spec isa Number
		return spec
	elseif spec isa AbstractVector
		@inbounds return spec[rand(1:length(spec))]
	elseif spec isa Tuple && length(spec) == 3 &&
		   spec[1] isa Number && spec[2] isa Number && spec[3] isa Integer
		lo, hi = spec[1], spec[2]
		return _bounded_draw(lo, hi, distribution)
	else
		if Base.iterable(spec)
			vals = collect(spec)
			@inbounds return vals[rand(1:length(vals))]
		end
		return spec
	end
end

# Collapse a (spec, pct) pair → (value::Number, pct::Union{Nothing,Number})
@inline function _collapse_pair(sp::Tuple, distribution::Symbol)
	spec, pct = sp
	v = _rand_in(spec, distribution)
	u = pct === nothing ? nothing : _rand_in(pct, distribution)
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
	collapse_cbs(cbs::CableBuilderSpec; distribution::Symbol = :uniform) -> CableBuilderSpec

Return a **singleton** `CableBuilderSpec` by collapsing every range-like item
(dims, args, and material fields) into one random draw using the chosen distribution.

"""
function collapse_cbs(
	cbs::CableBuilderSpec;
	distribution::Symbol = :uniform,
)
	parts = PartSpec[_collapse_part(p, distribution) for p in cbs.parts]
	return CableBuilderSpec(cbs.cable_id, parts, cbs.nominal)
end

"""
	sample(cbs::CableBuilderSpec; distribution::Symbol = :uniform) -> DataModel.CableDesign

Collapse ranges in `cbs` using `collapse_cbs` and build **one** cable design.
Useful for Monte Carlo style sampling where each call yields a new realization.
"""
function sample(
	cbs::CableBuilderSpec;
	distribution::Symbol = :uniform,
)
	scbs = collapse_cbs(cbs; distribution = distribution)
	designs = build(scbs)              # with singleton choices, this yields length == 1
	@assert length(designs) == 1
	return designs[1]
end

"""
	collapse_sbs(sbs::SystemBuilderSpec; distribution::Symbol = :uniform) -> SystemBuilderSpec

Collapse ranges in a `SystemBuilderSpec` using existing helpers.

Rules:
- Anchors `x, y` are numbers → pass through unchanged.
- `dx, dy` are ((nom_range), (unc_range)) → `_collapse_pair(_spec(...), distribution)`.
- `length`, `temperature` → `_collapse_pair(_spec(...), distribution)`.
- Earth fields (`rho`, `eps_r`, `mu_r`, `t`) → `_collapse_pair(_spec(...), distribution)`.
- Inner `builder` → `collapse_cbs(builder; distribution)`.
"""
function collapse_sbs(
	sbs::SystemBuilderSpec;
	distribution::Symbol = :uniform,
)
	scbs = collapse_cbs(sbs.builder; distribution = distribution)

	# positions: keep anchors; collapse dx/dy pairs
	pos = [
		begin
			dxc = _collapse_pair(_spec(p.dx), distribution)
			dyc = _collapse_pair(_spec(p.dy), distribution)
			typeof(p)(p.x0, p.y0, dxc, dyc, p.conn)
		end for p in sbs.positions
	]

	# system-level scalars-as-pairs
	L = _collapse_pair(_spec(sbs.length), distribution)
	T = _collapse_pair(_spec(sbs.temperature), distribution)

	# earth
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
		length = L,
		temperature = T,
		earth = earth,
		f = sbs.frequencies,
	)
end

"""
	sample(sbs::SystemBuilderSpec; distribution::Symbol = :uniform)

Collapse ranges in `sbs` and produce one `LineParametersProblem`.
"""
function sample(
	sbs::SystemBuilderSpec;
	distribution::Symbol = :uniform,
)
	ss = collapse_sbs(sbs; distribution = distribution)
	ch = iterate_problems(ss)
	return take!(ch)
end

include("montecarlo.jl")

end # module UQ