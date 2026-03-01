# The Vault (The Struct): Blunt, violently strict, purely parametric (T). Zero logic. Zero inner constructors. It just dictates the memory layout.
struct SolidCore{L, T <: Real} <: AbstractShape{L, T}
	r_ex::T
end

# The Janitor (The Outer Constructor): Takes the messy, mixed-type garbage from the user, sweeps it through promote(), and hands a perfectly clean, uniform tuple to the Vault.
# If the closure passes a Float64, T becomes Float64.
# If the closure passes a Measurement, T becomes Measurement.
SolidCore{L}(r_ex::T) where {L, T <: Real} = SolidCore{L, T}(r_ex)

# The Diplomat (Base.convert): When a parent container (like ConductorPart or your physics solver) demands that an already-built struct upgrade its precision (e.g., Float64 to Measurement), this hook effortlessly translates the fields and returns a new Vault.
function Base.convert(::Type{<:AbstractShape{L, T}}, s::SolidCore{L}) where {L, T <: Real}
	# Safely upgrade r_ex to the new T, and call the strictly-typed Vault
	return SolidCore{L, T}(convert(T, r_ex(s)))
end

@inline r_in(s::SolidCore) = zero(typeof(s.r_ex))
@inline r_ex(s::SolidCore) = s.r_ex

# ---------------------------------------------------------
# The Generic Solid Builder
# ---------------------------------------------------------
struct SolidCoreBuilder{P, T_geom <: Real, T_mat <: Real}
	tag::Symbol
	r_ex::T_geom
	mat::Material{T_mat}
end

@inline function SolidCoreBuilder{P}(
	tag::Symbol,
	r_ex::Tgeom,
	mat::Material{Tmat},
) where {P, Tgeom, Tmat}
	return SolidCoreBuilder{P, Tgeom, Tmat}(tag, r_ex, mat)
end

@inline function (b::SolidCoreBuilder{P})(current_r::T) where {P, T <: Real}
	current_r != zero(T) && error("Topological violation: Solid core must be at r=0.")
	return P(b.tag, SolidCore{Concentric}(b.r_ex), b.mat)
end

# ---------------------------------------------------------
# The  Solid blueprint
# ---------------------------------------------------------
struct SolidCoreSpec{P, Ttag, Tr, M <: AbstractSpec{Material}} <:
	   AbstractSpec{SolidCoreBuilder{P}}
	tag::Ttag
	r_ex::Tr
	mat::M
end

@inline function SolidCoreSpec(
	::Type{P},
	tag::Ttag,
	r_ex::Tr,
	mat::M,
) where
	{P, Ttag, Tr, M <: AbstractSpec{Material}}
	return SolidCoreSpec{P, Ttag, Tr, M}(tag, r_ex, mat)
end
