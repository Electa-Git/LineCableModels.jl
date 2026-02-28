# The Concrete Target (Temporary Stub)
struct CableDesign{T <: Tuple}
	payload::T
end

# Subtype it to your exact contract. Target = CableDesign.
struct CableDesignSpec{T <: Tuple} <: AbstractSpec{CableDesign}
	layers::T
end

# ---------------------------------------------------------
# Hooking into your Introspection Kernel
# ---------------------------------------------------------
# Override grid_args to unroll the tuple instead of inspecting struct fields
@inline grid_args(spec::CableDesignSpec) = spec.layers

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
