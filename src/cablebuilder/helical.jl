# ==========================================
# 1. THE VAULT (Fully resolved analytical geometry)
# ==========================================
struct HelicalPath{T <: Real, U <: Integer}
	ratio::T
	pitch::T
	angle::T
	overlength::T
	dir::U
end

function HelicalPath(ratio, pitch, angle, overlength, dir::U) where {U <: Integer}
	T = promote_type(typeof(ratio), typeof(pitch), typeof(angle), typeof(overlength))
	return HelicalPath{T, U}(
		convert(T, ratio), convert(T, pitch), convert(T, angle), convert(T, overlength), dir,
	)
end

# ==========================================
# 2. THE BUILDERS (Materialize the math using mean_diam)
# ==========================================
struct LayRatioBuilder{T <: Real, U <: Integer}
	ratio::T
	dir::U
end

@inline function (b::LayRatioBuilder)(mean_diam::T) where {T <: Real}
	pitch = b.ratio * mean_diam
	angle = atan(pi * mean_diam / pitch)
	overlength = sqrt(1 + (pi * mean_diam / pitch)^2)
	return HelicalPath(b.ratio, pitch, angle, overlength, b.dir)
end

# ==========================================
# 3. THE BLUEPRINTS (Specs for Combinatorics)
# ==========================================
struct LayRatioSpec{T, U} <: AbstractSpec{LayRatioBuilder}
	ratio::T
	dir::U
end