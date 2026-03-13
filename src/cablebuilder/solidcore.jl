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
@inline function build_part(
	::Type{Target},
	::Type{SolidCore},
	grp::Symbol,
	current_r::T,
	payload::Tuple,
) where {Target, T}
	mat, r = payload
	current_r != zero(T) && error("Core must be at r=0.")
	shape = SolidCore{Concentric}(current_r, current_r + r)
	return Target(grp, shape, mat) # Target is a compile-time constant here!
end
