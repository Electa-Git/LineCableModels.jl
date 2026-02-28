struct TubularShape{L, T <: Real} <: AbstractShape{L, T}
	r_in::T
	r_ex::T
end

function TubularShape{L}(r_in::Real, r_ex::Real) where {L}
	# Force the two boundaries to agree on a common precision T
	p = promote(r_in, r_ex)
	T_common = typeof(first(p))

	# Hand the perfectly aligned tuple to the auto-generated strict struct
	return TubularShape{L, T_common}(p...)
end

function Base.convert(
	::Type{<:AbstractShape{L, T}},
	s::TubularShape{L},
) where {L, T <: Real}
	# Safely cast both boundaries to the new T, and lock them in a new Vault
	return TubularShape{L, T}(convert(T, s.r_in), convert(T, s.r_ex))
end

# ---------------------------------------------------------
# The Payload (Does NOT subtype AbstractSpec)
# ---------------------------------------------------------
struct TubularBuilder{P, T_geom <: Real, T_mat <: Real}
	PartType::Type{P}  # Holds ConductorPart or InsulatorPart
	tag::Symbol
	t::T_geom
	mat::Material{T_mat}
end

# It waits peacefully until the materializer hands it current_r
function (b::TubularBuilder)(current_r::Real)
	r_ex = current_r + b.t

	# Stomp out the final physical struct using the role
	return b.PartType(b.tag, TubularShape{Concentric}(current_r, r_ex), b.mat)
end

# ---------------------------------------------------------
# The Blueprint (Subtypes AbstractSpec)
# ---------------------------------------------------------
struct TubularPartSpec{P, T, Th, M <: AbstractSpec{Material}} <:
	   AbstractSpec{TubularBuilder}
	PartType::P  # Will be a DeterministicGrid holding the Type
	tag::T
	t::Th
	mat::M
end