# The Concrete Target (Temporary Stub)
struct CableDesign{T <: Tuple}
	payload::T # WIP, only the payload for now
end

# Subtype it to your exact contract. Target = CableDesign.
struct CableDesignSpec{T <: Tuple} <: AbstractSpec{CableDesign}
	layers::T
end

# ---------------------------------------------------------
# Hooking into the Introspection Kernel
# ---------------------------------------------------------
# Override grid_args to unroll the tuple instead of inspecting struct fields
@inline grid_args(spec::CableDesignSpec) = spec.layers
# Specialize generator to avoid the extra layer of nesting/splat from the default generator
@inline generator(spec::AbstractSpec{CableDesign}) = (
	CableDesign(args) for args in Iterators.product(grid_args(spec)...)
)

# ---------------------------------------------------------
# Hooking into your Stochastic Kernel
# ---------------------------------------------------------
using Distributions: Distributions

# Map perfectly preserves tuple type-stability and splats into the Target
@inline function Base.rand(
	spec::CableDesignSpec;
	dist::Type{<:Distributions.ContinuousUnivariateDistribution} = Distributions.Normal,
)
	samples = map(l_spec -> rand(l_spec; dist = dist), spec.layers)
	return CableDesign(samples...)
end


# ---------------------------------------------------------
# The Cheap-Allocation Stacking Engine
# ---------------------------------------------------------

# Type-stable machinery to say that the cable center is at 0. Whether I'm being sarcastic or not depends on how much caffeine I've had.
@inline scalar_type(b) = typeof(_scalar_seed(b))

@inline _scalar_seed(b) = _scalar_seed_fields(b, Val(1), Val(fieldcount(typeof(b))))

@inline _scalar_seed_fields(::Any, ::Val{i}, ::Val{N}) where {i, N} =
	0.0  # should never happen unless no Real fields exist

# If we run out of fields, your builder has no Real-valued geometry at all.
# Which is… bold, from an engineering standpoint.
@inline _scalar_seed_fields(::Any, ::Val{0}) =
	Base.error(
		"Builder has no Real fields; cannot infer scalar type for stacking. Not to mention the engineering challenges involved in building some imaginary cable.",
	)

@inline function _scalar_seed_fields(b, ::Val{i}, ::Val{N}) where {i, N}
	x = getfield(b, i)
	x isa Real && return x
	i == N && return _scalar_seed_fields(b, Val(i+1), Val(N))  # trigger the error above
	return _scalar_seed_fields(b, Val(i+1), Val(N))
end

# The most overengineered 0 you'll see today. But it stacks cable layers, so shut up and enjoy.
@inline scalar_zero(b) = zero(scalar_type(b))

# Base case: no builders left
@inline build_layer(r, ::Tuple{}) = ()

# Recursive case: peel first + tail (Tuple recursion, because vectors are for people who enjoy allocations)
@inline function build_layer(r, builders::Tuple)
	b = first(builders)
	part = b(r)
	return (part, build_layer(r_ex(part.shape), Base.tail(builders))...)
end

@inline function CableDesign(builders::Tuple)
	z = scalar_zero(first(builders))   # "0", but in the right numeric universe
	parts = build_layer(z, builders)   # start stacking from the center like civilized cable geometry
	return CableDesign{typeof(parts)}(parts)
end

CableDesign(builders...) = CableDesign(builders)  # keep this non-inline if benchmarking