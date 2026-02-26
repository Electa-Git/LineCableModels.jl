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

"""
$(TYPEDEF)

Represents a concentric layer of rectangular strands with defined geometric, material, and electrical properties:

$(TYPEDFIELDS)
"""
struct RectStrands{T <: REALSCALAR, S <: RectStrandsShape} <: AbstractStrandsLayer{T}
	"Internal radial boundary \\[m\\]."
	radius_in::T
	"External radial boundary \\[m\\]."
	radius_ext::T
	"Material properties of the conductive strands."
	material_props::Material{T}
	"Operating temperature of the layer \\[°C\\]."
	temperature::T
	"Equivalent electrical resistance of the layer \\[Ω/m\\]."
	resistance::T
	"Geometric mean radius (GMR) of the layer \\[m\\]."
	gmr::T
	"Shape payload defining the internal geometric layout."
	shape::S
end

# struct SectorCore{T <: REALSCALAR, S <: SectorShape} <: AbstractStrandsLayer{T}
# 	"Internal radial boundary \\[m\\]."
# 	radius_in::T
# 	"External radial boundary \\[m\\]."
# 	radius_ext::T
# 	"Material properties of the conductive strands."
# 	material_props::Material{T}
# 	"Operating temperature of the layer \\[°C\\]."
# 	temperature::T
# 	"Equivalent electrical resistance of the layer \\[Ω/m\\]."
# 	resistance::T
# 	"Geometric mean radius (GMR) of the layer \\[m\\]."
# 	gmr::T
# 	"Shape payload defining the internal geometric layout."
# 	shape::S
# end

# struct CircCore{T <: REALSCALAR, S <: Concentric} <: AbstractStrandsLayer{T}
# 	"Internal radial boundary \\[m\\]."
# 	radius_in::T
# 	"External radial boundary \\[m\\]."
# 	radius_ext::T
# 	"Material properties of the conductive strands."
# 	material_props::Material{T}
# 	"Operating temperature of the layer \\[°C\\]."
# 	temperature::T
# 	"Equivalent electrical resistance of the layer \\[Ω/m\\]."
# 	resistance::T
# 	"Geometric mean radius (GMR) of the layer \\[m\\]."
# 	gmr::T
# 	"Shape payload defining the internal geometric layout."
# 	shape::S
# end

"""
$(TYPEDSIGNATURES)

Constructs a [`RectStrands`](@ref) instance.

# Arguments

- `radius_in`: Internal radius of the layer \\[m\\].
- `thickness`: Radial thickness of the strands \\[m\\].
- `width`: Width of the individual rectangular strip \\[m\\].
- `num_wires`: Number of strands in the layer \\[dimensionless\\].
- `lay_ratio`: Ratio defining the lay length of the strands \\[dimensionless\\].
- `material_props`: A [`Material`](@ref) object containing physical properties.
- `temperature`: Operating temperature \\[°C\\].
- `lay_direction`: Twisting direction (1 = unilay, -1 = contralay) \\[dimensionless\\].

# Returns

- A [`RectStrands`](@ref) object with calculated geometric and electrical properties.

# Examples

```julia
material = Material(1.724e-8, 1.0, 1.0, 20.0, 0.00393)
layer = $(FUNCTIONNAME)(0.01, 0.002, 0.005, 10, 12.0, material, 25.0, 1)
```
"""
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


	# Target Area (the 'invariant' property)
	A0 = width * thickness

	# Area-Preserving Radial Expansion
	# This solves A_total = pi * (r_ext^2 - r_in^2)
	radius_ext =
		num_wires == 1 ? Base.error("num_wires must be > 1") :
		sqrt(radius_in^2 + (num_wires * A0) / T(π))
	thickness_effective=radius_ext-radius_in
	@info "Calculating outer radius to preserve total cross-sectional area of strands." radius_ext radius_ext_without_correction=radius_in+thickness thickness_effective thickness num_wires
	mean_diameter, pitch_length, overlength =
		calc_helical_params(radius_in, radius_ext, lay_ratio)

	cross_section = num_wires * A0

	shape_payload = RectStrandsShape(
		thickness_effective, width, num_wires, lay_ratio, lay_direction,
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