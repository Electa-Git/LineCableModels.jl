# ==========================================
# 1. THE VAULT
# ==========================================
struct TubularLayer{T <: Real} <: AbstractShape{T}
	r_in::T
	r_ex::T
end

function TubularLayer(r_in, r_ex)
	T = promote_type(typeof(r_in), typeof(r_ex))
	return TubularLayer{T}(convert(T, r_in), convert(T, r_ex))
end

function Base.convert(
	::Type{<:AbstractShape{T}},
	s::TubularLayer,
) where {T <: Real}
	return TubularLayer{T}(convert(T, s.r_in), convert(T, s.r_ex))
end

# ==========================================
# 2. THE BUILDER
# ==========================================
# Holds the abstract thickness value.
struct TubularLayerBuilder{P, T <: Real, M <: Real}
	cmp::Symbol
	t::T
	mat::Material{M}
end

@inline function TubularLayerBuilder{P}(
	cmp::Symbol,
	t::T,
	mat::Material{M},
) where {P, T, M}
	return TubularLayerBuilder{P, T, M}(cmp, t, mat)
end

@inline function (b::TubularLayerBuilder{P})(current_r::T) where {P, T <: Real}
	# Physics reminder: r_ex = r_in + thickness. 
	r_ex = current_r + b.t
	shape = TubularLayer(current_r, r_ex)
	return P(b.cmp, shape, b.mat)
end

# ==========================================Tr
# 3. THE BLUEPRINT
# ==========================================
# Flat fields. Measurements.jl will propagate the uncertainty of `t` cleanly.
struct TubularLayerSpec{P, G, T, M <: AbstractSpec{Material}} <:
	   AbstractSpec{TubularLayerBuilder{P}}
	cmp::G
	t::T
	mat::M
end

@inline function TubularLayerSpec(
	::Type{P},
	cmp::G,
	t::T,
	mat::M,
) where {P, G, T, M <: AbstractSpec{Material}}
	return TubularLayerSpec{P, G, T, M}(cmp, t, mat)
end
