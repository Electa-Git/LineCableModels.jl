# The Vault
struct CircStrandedShape{L, T <: Real, U <: Integer} <: AbstractShape{L, T}
	r_in::T
	r_ex::T
	r_w::T
	n_w::U
	lay_r::T
	lay_d::U
end

# The Janitor
function CircStrandedShape{L}(
	r_in, r_ex, r_w, n_w::U, lay_r, lay_d::U,
) where {L, U <: Integer}
	T = promote_type(typeof(r_in), typeof(r_ex), typeof(r_w), typeof(lay_r))
	return CircStrandedShape{L, T, U}(
		convert(T, r_in),
		convert(T, r_ex),
		convert(T, r_w),
		n_w,
		convert(T, lay_r),
		lay_d,
	)
end

# The Diplomat
function Base.convert(
	::Type{<:AbstractShape{L, T}},
	s::CircStrandedShape{L},
) where {L, T <: Real}
	return CircStrandedShape{L, T, typeof(s.n_w)}(
		convert(T, s.r_in),
		convert(T, s.r_ex),
		convert(T, s.r_w),
		s.n_w,
		convert(T, s.lay_r),
		s.lay_d,
	)
end

struct CircStrandedBuilder{P, Tgeom <: Real, Tmat <: Real, U <: Integer}
	cmp::Symbol
	r_w::Tgeom
	n_w::U
	lay_r::Tgeom
	lay_d::U
	mat::Material{Tmat}
end

@inline function CircStrandedBuilder{P}(
	cmp::Symbol, r_w::Tgeom, n_w::U, lay_r::Tgeom, lay_d::U, mat::Material{Tmat},
) where {P, Tgeom, Tmat, U <: Integer}
	return CircStrandedBuilder{P, Tgeom, Tmat, U}(cmp, r_w, n_w, lay_r, lay_d, mat)
end

@inline function (b::CircStrandedBuilder{P})(current_r::T) where {P, T <: Real}
	# Strict geometric continuation. r_ex is fully determined by current_r + diameter.
	r_ex = current_r + 2 * b.r_w
	shape = CircStrandedShape{Concentric}(current_r, r_ex, b.r_w, b.n_w, b.lay_r, b.lay_d)

	# P is the ConductorPart. This call will trigger the T = promote_type(...) rule.
	return P(b.cmp, shape, b.mat)
end

struct CircStrandedSpec{P, Tcmp, Trw, Tnw, Tlayr, Tlayd, M <: AbstractSpec{Material}} <:
	   AbstractSpec{CircStrandedBuilder{P}}
	cmp::Tcmp
	r_w::Trw
	n_w::Tnw
	lay_r::Tlayr
	lay_d::Tlayd
	mat::M
end

@inline function CircStrandedSpec(
	::Type{P}, cmp::Tcmp, r_w::Trw, n_w::Tnw, lay_r::Tlayr, lay_d::Tlayd, mat::M,
) where {P, Tcmp, Trw, Tnw, Tlayr, Tlayd, M <: AbstractSpec{Material}}
	return CircStrandedSpec{P, Tcmp, Trw, Tnw, Tlayr, Tlayd, M}(
		cmp, r_w, n_w, lay_r, lay_d, mat,
	)
end
