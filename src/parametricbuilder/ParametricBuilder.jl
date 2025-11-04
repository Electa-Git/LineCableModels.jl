module ParametricBuilder

# Export public API
export make_stranded, make_screened
export Conductor, Insulator, Material, CableBuilder
export at, Earth, SystemBuilder
export determinize, sample

# Module-specific dependencies
using ..Commons
import ..Commons: add!
using ..Materials: Materials
using ..DataModel: DataModel
using ..EarthProps: EarthModel
using ..Engine: LineParametersProblem
using Measurements
using Base.Iterators: product

# normalize input to (spec, pct)
_spec(x) = (x isa Tuple && length(x)==2) ? x : (x, nothing)


_values(x::Number) = (x,)
_values(v::AbstractVector) = collect(v)
_values(t::Tuple{<:Number, <:Number, <:Integer}) = range(t[1], t[2]; length = t[3])

_pcts(::Nothing) = (0.0,)
_pcts(p::Number) = (float(p),)
_pcts(v::AbstractVector) = map(float, collect(v))
_pcts(t::Tuple{<:Number, <:Number, <:Integer}) =
	range(float(t[1]), float(t[2]); length = t[3])

function _make_range(spec; pct = nothing)
	vs, ps = collect(_values(spec)), collect(_pcts(pct))
	if all(p->p==0.0, ps)
		;
		return vs;
	end
	out = Any[]
	for v in vs, p in ps
		push!(out, measurement(v, abs(v)*(p/100)))
	end
	out
end

# expand positional args tuple → iterator of resolved tuples
function _expand_args(args::Tuple)
	spaces =
		map(a -> (a isa Tuple && length(a)==2 ? _make_range(a[1]; pct = a[2]) : (a,)), args)
	return (tuple(vals...) for vals in Iterators.product(spaces...))
end

include("materialspec.jl")
include("cablebuilderspec.jl")
include("systembuilderspec.jl")
include("determinize.jl")
include("uquant.jl")
include("base.jl")

# Submodule `WirePatterns`
include("wirepatterns/WirePatterns.jl")
using .WirePatterns

end # module ParametricBuilder
