"""
$(TYPEDEF)

Holds the pure geometric layout for a concentric layer of rectangular/flat strands.
"""
struct RectStrandsShape{T <: REALSCALAR, U <: Int} <: AbstractShapeGeometry
	"Thickness of the individual rectangular strip \\[m\\]."
	thickness::T
	"Width of the individual rectangular strip \\[m\\]."
	width::T
	"Number of wires in the layer \\[dimensionless\\]."
	num_wires::U
	"Ratio defining the lay length of the wires \\[dimensionless\\]."
	lay_ratio::T
	"Twisting direction of the strands (1 = unilay, -1 = contralay) \\[dimensionless\\]."
	lay_direction::U
	"Mean diameter of the wire layer \\[m\\]."
	mean_diameter::T
	"Pitch length of the wire layer \\[m\\]."
	pitch_length::T
	"Total cross-sectional area of the conductive metal in this layer \\[m²\\]."
	cross_section::T
end

struct RectStrands{T <: REALSCALAR, S <: RectStrandsShape} <: AbstractStrandsLayer{T}
	"Internal radial boundary \\[m\\]."
	radius_in::T
	"External radial boundary \\[m\\]."
	radius_ext::T
	material_props::Material{T}
	temperature::T

	# Fundamental electrical properties for ConductorGroup
	resistance::T
	gmr::T

	# Shape payload, defines how resistance and gmr are computed
	shape::S
end

function RectStrands(
	radius_in::T,
	thickness::T,
	width::T,
	num_wires::U,
	lay_ratio::T,
	material_props::Material{T},
	temperature::T,
	lay_direction::U,
) where {T <: REALSCALAR, U <: Int}

	# Geometry calculations
	radius_ext =
		num_wires == 1 ? Base.error("num_wires must be > 1") : radius_in + thickness
	mean_diameter, pitch_length, overlength =
		calc_helical_params(radius_in, radius_ext, lay_ratio)

	cross_section = num_wires * (thickness * width)

	shape_payload = RectStrandsShape(
		thickness, width, num_wires, lay_ratio, lay_direction,
		mean_diameter, pitch_length, cross_section,
	)

	# Electrical properties
	rho = material_props.rho
	T0 = material_props.T0
	alpha = material_props.alpha

	R_wire =
		calc_strip_resistance(thickness, width, rho, alpha, T0, temperature) * overlength
	R_layer = R_wire / num_wires

	gmr_layer = calc_tubular_gmr(radius_ext, radius_in, material_props.mu_r)

	# 3. Instantiate the concrete struct
	return RectStrands(
		radius_in,
		radius_ext,
		material_props,
		temperature,
		R_layer,
		gmr_layer,
		shape_payload,
	)
end

const _REQ_RECTSTRANDS =
	(:radius_in, :thickness, :width, :num_wires, :lay_ratio, :material_props)
const _OPT_RECTSTRANDS = (:temperature, :lay_direction)
const _DEFS_RECTSTRANDS = (T₀, 1)

Validation.has_radii(::Type{RectStrands}) = false
Validation.has_temperature(::Type{RectStrands}) = true
Validation.required_fields(::Type{RectStrands}) = _REQ_RECTSTRANDS
Validation.keyword_fields(::Type{RectStrands}) = _OPT_RECTSTRANDS
Validation.keyword_defaults(::Type{RectStrands}) = _DEFS_RECTSTRANDS

Validation.coercive_fields(::Type{RectStrands}) =
	(:radius_in, :thickness, :width, :lay_ratio, :material_props, :temperature)

Validation.is_radius_input(
	::Type{RectStrands},
	::Val{:radius_in},
	x::AbstractCablePart,
) = true

Validation.is_radius_input(::Type{RectStrands}, ::Val{:radius_in}, x::Thickness) = true

# Specific for rectangular strands
Validation.maxfill(::Type{RectStrands}, rin::Real, w::Real) =
	floor(Int, 2 * π * rin / w)

Validation.extra_rules(::Type{RectStrands}) = (
	Normalized(:radius_in), Finite(:radius_in), Nonneg(:radius_in),
	Normalized(:thickness), Finite(:thickness), Positive(:thickness),
	IntegerField(:num_wires), Positive(:num_wires),
	Finite(:lay_ratio), Nonneg(:lay_ratio),
	IsA{Material}(:material_props),
	OneOf(:lay_direction, (-1, 1)), Finite(:width),
	Positive(:width),
	PhysicalFillLimit(:num_wires, (:radius_in, :width)), # THE BOUNCER
)

Validation.parse(::Type{RectStrands}, nt) = begin
	rin, rw = _normalize_radii(RectStrands, nt.radius_in, nt.thickness)

	# Resolves to Int using the generic interface
	n_wires = _resolve_strands(nt.num_wires, RectStrands, rin, nt.width)

	return (; nt..., radius_in = rin, thickness = rw, num_wires = n_wires)
end

@construct RectStrands _REQ_RECTSTRANDS _OPT_RECTSTRANDS _DEFS_RECTSTRANDS