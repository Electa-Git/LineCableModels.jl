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
	return SolidCore{L, T}(convert(T, s.r_ex))
end

# ---------------------------------------------------------
# The Generic Solid Builder
# ---------------------------------------------------------
struct SolidCoreBuilder{P, T_geom <: Real, T_mat <: Real}
	PartType::Type{P}
	tag::Symbol
	r_ex::T_geom
	mat::Material{T_mat}
end

function (b::SolidCoreBuilder)(current_r::Real)
	if current_r != 0.0
		error("Topological violation: Solid core must be at r=0.")
	end
	return b.PartType(b.tag, SolidCore{Concentric}(b.r_ex), b.mat)
end

# ---------------------------------------------------------
# The  Solid blueprint
# ---------------------------------------------------------
struct SolidCoreSpec{P, T, R, M <: AbstractSpec{Material}} <: AbstractSpec{SolidCoreBuilder}
	PartType::P
	tag::T
	r_ex::R
	mat::M
end