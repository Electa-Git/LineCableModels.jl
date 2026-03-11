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
struct SolidCoreBuilder{P, Tr <: Real, Tmat <: Real}
	cmp::Symbol
	r::Tr
	mat::Material{Tmat}
end

@inline function SolidCoreBuilder{P}(
	cmp::Symbol,
	r::Tr,
	mat::Material{Tmat},
) where {P, Tr, Tmat}
	return SolidCoreBuilder{P, Tr, Tmat}(cmp, r, mat)
end

@inline function (b::SolidCoreBuilder{P})(current_r::T) where {P, T <: Real}
	# If someone tries to stack a solid core on top of an existing layer, mock them.
	current_r != zero(T) && error(
		"Topological violation: Solid core must be at r=0. You can't put a solid core on the outside of a cable.",
	)

	shape = SolidCore{Concentric}(current_r, b.r)
	return P(b.cmp, shape, b.mat)
end

# ==========================================
# 3. THE BLUEPRINT
# ==========================================
# Flat fields. Perfectly aligned with spec.jl's Introspection Kernel.
struct SolidCoreSpec{P, Tcmp, Tr, M <: AbstractSpec{Material}} <:
	   AbstractSpec{SolidCoreBuilder{P}}
	cmp::Tcmp
	r::Tr
	mat::M
end

@inline function SolidCoreSpec(
	::Type{P},
	cmp::Tcmp,
	r::Tr,
	mat::M,
) where {P, Tcmp, Tr, M <: AbstractSpec{Material}}
	return SolidCoreSpec{P, Tcmp, Tr, M}(cmp, r, mat)
end
