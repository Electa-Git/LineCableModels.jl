# ==========================================
# 1. THE VAULT
# ==========================================
struct TubularLayer{L, T <: Real} <: AbstractShape{L, T}
	r_in::T
	r_ex::T
end

function TubularLayer{L}(r_in, r_ex) where {L}
	T = promote_type(typeof(r_in), typeof(r_ex))
	return TubularLayer{L, T}(convert(T, r_in), convert(T, r_ex))
end

function Base.convert(
	::Type{<:AbstractShape{L, T}},
	s::TubularLayer{L},
) where {L, T <: Real}
	return TubularLayer{L, T}(convert(T, s.r_in), convert(T, s.r_ex))
end

# ==========================================
# 2. THE BUILDER
# ==========================================
# Holds the abstract thickness value.
struct TubularLayerBuilder{P, Tt <: Real, Tmat <: Real}
	cmp::Symbol
	t::Tt
	mat::Material{Tmat}
end

@inline function TubularLayerBuilder{P}(
	cmp::Symbol,
	t::Tt,
	mat::Material{Tmat},
) where {P, Tt, Tmat}
	return TubularLayerBuilder{P, Tt, Tmat}(cmp, t, mat)
end

@inline function (b::TubularLayerBuilder{P})(current_r::T) where {P, T <: Real}
	# Physics reminder: r_ex = r_in + thickness. 
	r_ex = current_r + b.t
	shape = TubularLayer{Concentric}(current_r, r_ex)
	return P(b.cmp, shape, b.mat)
end

# ==========================================
# 3. THE BLUEPRINT
# ==========================================
# Flat fields. Measurements.jl will propagate the uncertainty of `t` cleanly.
struct TubularLayerSpec{P, Tcmp, Tt, M <: AbstractSpec{Material}} <:
	   AbstractSpec{TubularLayerBuilder{P}}
	cmp::Tcmp
	t::Tt
	mat::M
end

@inline function TubularLayerSpec(
	::Type{P},
	cmp::Tcmp,
	t::Tt,
	mat::M,
) where {P, Tcmp, Tt, M <: AbstractSpec{Material}}
	return TubularLayerSpec{P, Tcmp, Tt, M}(cmp, t, mat)
end
