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
# 2. THE BUILDER
# ==========================================
# Holds the actual target radius grid-value.
struct SolidCoreBuilder{P, T <: Real, M <: Real}
	grp::Symbol
	r::T
	mat::Material{M}
end

@inline function SolidCoreBuilder{P}(
	grp::Symbol,
	r::T,
	mat::Material{M},
) where {P, T, M}
	return SolidCoreBuilder{P, T, M}(grp, r, mat)
end

@inline function (b::SolidCoreBuilder{P})(current_r::T) where {P, T <: Real}
	# If someone tries to stack a solid core on top of an existing layer, mock them.
	current_r != zero(T) && error(
		"Topological violation: Solid core must be at r=0. You can't put a solid core on the outside of a cable.",
	)

	shape = SolidCore{Concentric}(current_r, b.r)
	return P(b.grp, shape, b.mat)
end

# ==========================================
# 3. THE BLUEPRINT
# ==========================================
# Flat fields. Perfectly aligned with spec.jl's Introspection Kernel.
struct SolidCoreSpec{P, G, T, M <: AbstractSpec{Material}} <:
	   AbstractSpec{SolidCoreBuilder{P}}
	grp::G
	r::T
	mat::M
end

@inline function SolidCoreSpec(
	::Type{P},
	grp::G,
	r::T,
	mat::M,
) where {P, G, T, M <: AbstractSpec{Material}}
	return SolidCoreSpec{P, G, T, M}(grp, r, mat)
end
