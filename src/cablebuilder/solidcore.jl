# ==========================================
# 1. THE VAULT
# ==========================================
# Just strictly holds the universal boundaries.
struct SolidCore{L, T <: Real} <: AbstractShape{L, T}
	r_in::T
	r_ex::T
end

function SolidCore{L}(r_in, r_ex) where {L}
	T = promote_type(typeof(r_in), typeof(r_ex))
	return SolidCore{L, T}(convert(T, r_in), convert(T, r_ex))
end

function Base.convert(::Type{<:AbstractShape{L, T}}, s::SolidCore{L}) where {L, T <: Real}
	return SolidCore{L, T}(convert(T, s.r_in), convert(T, s.r_ex))
end

# ==========================================
# 2. THE FUNCTOR SPECIALIZATION (Spatial Collapse)
# ==========================================
# We dispatch exactly on the ShapeTag. 
@inline function (b::PartBuilder{Part, SolidCore})(current_r::T) where {Part, T <: Real}
	# Unpack the payload exactly in the positional order you defined for this shape
	mat, r = b.payload

	current_r != zero(T) && error("Topological violation: Solid core must be at r=0.")

	shape = SolidCore{Concentric}(current_r, current_r + r)

	# Natively build the final role object (e.g., ConductorPart)
	return b.part(b.grp, shape, mat)
end
