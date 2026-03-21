# ==========================================
# 1. THE UNIVERSAL VAULT
# ==========================================
struct StrandedLayer{
	L,
	T <: Real,
	U <: Integer,
	P <: AbstractWire{T},
	H <: HelicalPath,
} <: AbstractShape{T}
	r_in::T
	r_ex::T
	n_w::U
	wire::P
	pitch::H
end

function StrandedLayer(
	r_in,
	r_ex,
	n_w::Integer,
	wire::AbstractWire,
	pitch::HelicalPath,
)
	T = promote_type(typeof(r_in), typeof(r_ex))
	p = convert(AbstractWire{T}, wire)
	return StrandedLayer{T, typeof(n_w), typeof(p), typeof(pitch)}(
		convert(T, r_in), convert(T, r_ex), n_w, p, pitch,
	)
end

function Base.convert(
	::Type{<:AbstractShape{T}},
	s::StrandedLayer,
) where {T <: Real}
	p_converted = convert(AbstractWire{T}, s.wire)
	return StrandedLayer{T, typeof(s.n_w), typeof(p_converted), typeof(s.pitch)}(
		convert(T, s.r_in), convert(T, s.r_ex), s.n_w, p_converted, s.pitch,
	)
end

# ==========================================
# 2. THE UNIVERSAL BUILDER
# ==========================================
struct StrandedBuilder{P, U <: Integer, W, H, T <: Real}
	cmp::Symbol
	n_w::U
	wire_builder::W
	pitch_builder::H
	mat::Material{T}
end

@inline function (b::StrandedBuilder{P})(current_r::T) where {P, T <: Real}
	# If someone tries to put a stranded armor at the exact center of the universe, mock them.
	current_r <= zero(T) && error(
		"Topological violation: Stranded layers cannot exist at r=0. Use a SolidCore.",
	)

	# 1. Materialize the physical entity
	wire = b.wire_builder()

	# 2. Extract its radial footprint via dispatch
	thick = char_len(wire)

	r_ex = current_r + thick
	mean_diam = current_r + (thick / 2)

	# 3. Pass context to the nested helical builder
	pitch_profile = b.pitch_builder(mean_diam)

	# 4. Lock it into the unified layer
	shape = StrandedLayer(current_r, r_ex, b.n_w, wire, pitch_profile)

	return P(b.cmp, shape, b.mat)
end

# ==========================================
# 3. THE UNIVERSAL BLUEPRINT
# ==========================================
struct StrandedSpec{
	P,
	G,
	U,
	W <: AbstractSpec,
	H <: AbstractSpec,
	M <: AbstractSpec{Material},
} <: AbstractSpec{StrandedBuilder{P}}
	cmp::G
	n_w::U
	wire_spec::W
	pitch_spec::H
	mat::M
end