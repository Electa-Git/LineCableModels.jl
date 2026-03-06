# The Vault (The Struct): Blunt, violently strict, purely parametric (T). Zero logic. Zero inner constructors. It just dictates the memory layout.
struct SolidShape{L, T <: Real} <: AbstractShape{L, T}
	r_ex::T
end

# The Janitor (The Outer Constructor): Takes the messy, mixed-type garbage from the user, sweeps it through promote(), and hands a perfectly clean, uniform tuple to the Vault.
# If the closure passes a Float64, T becomes Float64.
# If the closure passes a Measurement, T becomes Measurement.
SolidShape{L}(r_ex::T) where {L, T <: Real} = SolidShape{L, T}(r_ex)

# The Diplomat (Base.convert): When a parent container (like ConductorPart or your physics solver) demands that an already-built struct upgrade its precision (e.g., Float64 to Measurement), this hook effortlessly translates the fields and returns a new Vault.
function Base.convert(::Type{<:AbstractShape{L, T}}, s::SolidShape{L}) where {L, T <: Real}
	# Safely upgrade r_ex to the new T, and call the strictly-typed Vault
	return SolidShape{L, T}(convert(T, r_ex(s)))
end

# Override the universal accessor because SolidShape has no r_in field
@inline r_in(s::SolidShape) = zero(typeof(s.r_ex))
# r_ex(s) implicitly falls back to s.r_ex defined on AbstractShape

# ---------------------------------------------------------
# The Generic Solid Builder
# ---------------------------------------------------------
struct SolidCoreBuilder{P, T_geom <: Real, T_mat <: Real}
	cmp::Symbol
	r_ex::T_geom
	mat::Material{T_mat}
end

@inline function SolidCoreBuilder{P}(
	cmp::Symbol,
	r_ex::Tgeom,
	mat::Material{Tmat},
) where {P, Tgeom, Tmat}
	return SolidCoreBuilder{P, Tgeom, Tmat}(cmp, r_ex, mat)
end

@inline function (b::SolidCoreBuilder{P})(current_r::T) where {P, T <: Real}
	current_r != zero(T) && error("Topological violation: Solid core must be at r=0.")
	return P(b.cmp, SolidShape{Concentric}(b.r_ex), b.mat)
end

# ---------------------------------------------------------
# The  Solid blueprint
# ---------------------------------------------------------
struct SolidCoreSpec{P, Tcmp, Tr, M <: AbstractSpec{Material}} <:
	   AbstractSpec{SolidCoreBuilder{P}}
	cmp::Tcmp
	r_ex::Tr
	mat::M
end

@inline function SolidCoreSpec(
	::Type{P},
	cmp::Tcmp,
	r_ex::Tr,
	mat::M,
) where {P, Tcmp, Tr, M <: AbstractSpec{Material}}
	return SolidCoreSpec{P, Tcmp, Tr, M}(cmp, r_ex, mat)
end
