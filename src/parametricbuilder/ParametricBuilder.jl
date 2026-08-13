module ParametricBuilder

# Export public API
export make_stranded, make_screened
export Conductor, Insulator, Material, CableBuilder
export at, trifoil, Earth, SystemBuilder
export determinize

# Module-specific dependencies
using ..Commons
import ..Commons: add!
using ..Materials: Materials
using ..DataModel: DataModel, trifoil_formation, CableDesign, get_outer_radius
using ..EarthProps: EarthModel
using ..Engine: LineParametersProblem
using ..Utils: to_nominal
using Measurements
using Base.Iterators: product



# normalize input to (spec, pct)
function _spec(x)
	if x isa Tuple && length(x) == 2
		spec, pct = x
		_validate_valuespec(spec)
		_validate_valuespec(pct)
		return (spec, pct)
	else
		_validate_valuespec(x)
		return (x, nothing)
	end
end

# THIS is the minimal validator 
@inline function _validate_valuespec(x)
	if x isa Tuple && length(x) == 3 && x[1] isa Number && x[2] isa Number &&
	   x[3] isa Integer
		lo, hi, n = x
		n >= 2 || error("Range (lo,hi,n) must have n ≥ 2, got $x")
	end
	return nothing
end


_values(x::Number) = (x,)
_values(v::AbstractVector) = collect(v)
_values(t::Tuple{<:Number, <:Number, <:Integer}) = range(t[1], t[2]; length = t[3])


_pcts(::Nothing) = (0.0,)
_pcts(p::Number) = (float(p),)
_pcts(v::AbstractVector) = map(float, collect(v))
_pcts(t::Tuple{<:Number, <:Number, <:Integer}) = range(t[1], t[2]; length = t[3])


function _make_range(spec; pct = nothing)

	if spec isa Tuple && length(spec)==3
		lo, hi, n = spec
		n < 2 && Base.error("Invalid (lo,hi,n) values range: n=$n must be ≥2")
	end

	vs, ps = collect(_values(spec)), collect(_pcts(pct))
	if all(p->p==0.0, ps)
		return vs
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
include("base.jl")

# Submodule `WirePatterns`
include("wirepatterns/WirePatterns.jl")
using .WirePatterns

end # module ParametricBuilder
