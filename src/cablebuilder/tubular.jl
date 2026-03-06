struct TubularShape{L, T <: Real} <: AbstractShape{L, T}
	r_in::T
	r_ex::T
end

function TubularShape{L}(r_in, r_ex) where {L}
	T = promote_type(typeof(r_in), typeof(r_ex))
	return TubularShape{L, T}(convert(T, r_in), convert(T, r_ex))
end

function Base.convert(
	::Type{<:AbstractShape{L, T}},
	s::TubularShape{L},
) where {L, T <: Real}
	# Safely cast both boundaries to the new T, and lock them in a new Vault
	return TubularShape{L, T}(convert(T, r_in(s)), convert(T, r_ex(s)))
end

# ---------------------------------------------------------
# The Payload (Does NOT subtype AbstractSpec)
# ---------------------------------------------------------
struct TubularBuilder{P, Tgeom <: Real, Tmat <: Real}
	cmp::Symbol
	t::Tgeom
	mat::Material{Tmat}
end

@inline function TubularBuilder{P}(
	cmp::Symbol,
	t::Tgeom,
	mat::Material{Tmat},
) where {P, Tgeom, Tmat}
	return TubularBuilder{P, Tgeom, Tmat}(cmp, t, mat)
end

# # It waits peacefully until the materializer hands it current_r
@inline function (b::TubularBuilder{P})(current_r::T) where {P, T <: Real}
	r_ex = current_r + b.t
	return P(b.cmp, TubularShape{Concentric}(current_r, r_ex), b.mat)
end

# ---------------------------------------------------------
# The Blueprint (Subtypes AbstractSpec)
# ---------------------------------------------------------
struct TubularLayerSpec{P, Tcmp, Tt, M <: AbstractSpec{Material}} <:
	   AbstractSpec{TubularBuilder{P}}
	cmp::Tcmp
	t::Tt
	mat::M
end

@inline function TubularLayerSpec(
	::Type{P},
	cmp::Tcmp,
	t::Tt,
	mat::M,
) where
	{P, Tcmp, Tt, M <: AbstractSpec{Material}}
	return TubularLayerSpec{P, Tcmp, Tt, M}(cmp, t, mat)
end
