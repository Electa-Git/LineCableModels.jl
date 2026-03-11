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
# Base case: no builders left
@inline build_layer(r, ::Tuple{}) = ()

# Recursive case: peel first + tail 
@inline function build_layer(r, builders::Tuple)
	b = first(builders)
	part = b(r) # The Shape's Janitor handles promotion natively here
	return (part, build_layer(r_ex(part.shape), Base.tail(builders))...)
end

@inline function CableDesign(builders::Tuple)
	# Toss a default Float64 zero and let Julia's tuple recursion 
	# and promote_type handle the rest.
	parts = build_layer(0.0, builders)
	return CableDesign{typeof(parts)}(parts)
end

CableDesign(builders...) = CableDesign(builders)  # keep this non-inline if benchmarking