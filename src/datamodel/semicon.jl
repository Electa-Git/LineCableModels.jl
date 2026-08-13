"""
$(TYPEDEF)

Represents a semiconducting layer with defined geometric, material, and electrical properties given by the attributes:

$(TYPEDFIELDS)
"""
struct Semicon{T <: REALSCALAR} <: AbstractInsulatorPart{T}
	"Internal radius of the semiconducting layer \\[m\\]."
	r_in::T
	"External radius of the semiconducting layer \\[m\\]."
	r_ex::T
	"Material properties of the semiconductor."
	material_props::Material{T}
	"Operating temperature of the semiconductor \\[°C\\]."
	temperature::T
	"Cross-sectional area of the semiconducting layer \\[m²\\]."
	cross_section::T
	"Electrical resistance of the semiconducting layer \\[Ω/m\\]."
	resistance::T
	"Geometric mean radius of the semiconducting layer \\[m\\]."
	gmr::T
	"Shunt capacitance per unit length of the semiconducting layer \\[F/m\\]."
	shunt_capacitance::T
	"Shunt conductance per unit length of the semiconducting layer \\[S·m\\]."
	shunt_conductance::T
end

"""
$(TYPEDSIGNATURES)

Constructs a [`Semicon`](@ref) instance with calculated electrical and geometric properties.

# Arguments

- `r_in`: Internal radius of the semiconducting layer \\[m\\].
- `r_ex`: External radius or thickness of the layer \\[m\\].
- `material_props`: Material properties of the semiconducting material.
- `temperature`: Operating temperature of the layer \\[°C\\] (default: T₀).

# Returns

- A [`Semicon`](@ref) object with initialized properties.

# Examples

```julia
material_props = Material(1e6, 2.3, 1.0, 20.0, 0.00393)
semicon_layer = $(FUNCTIONNAME)(0.01, Thickness(0.002), material_props, temperature=25)
println(semicon_layer.cross_section)      # Expected output: ~6.28e-5 [m²]
println(semicon_layer.resistance)         # Expected output: Resistance in [Ω/m]
println(semicon_layer.gmr)                # Expected output: GMR in [m]
println(semicon_layer.shunt_capacitance)  # Expected output: Capacitance in [F/m]
println(semicon_layer.shunt_conductance)  # Expected output: Conductance in [S·m]
```
"""
function Semicon(
	r_in::T,
	r_ex::T,
	material_props::Material{T},
	temperature::T,
) where {T <: REALSCALAR}

	rho = material_props.rho
	T0 = material_props.T0
	alpha = material_props.alpha
	epsr_r = material_props.eps_r

	cross_section = π * (r_ex^2 - r_in^2)

	resistance =
		calc_tubular_resistance(r_in, r_ex, rho, alpha, T0, temperature)
	gmr = calc_tubular_gmr(r_ex, r_in, material_props.mu_r)
	shunt_capacitance = calc_shunt_capacitance(r_in, r_ex, epsr_r)
	shunt_conductance = calc_shunt_conductance(r_in, r_ex, rho)

	# Initialize object
	return Semicon(
		r_in,
		r_ex,
		material_props,
		temperature,
		cross_section,
		resistance,
		gmr,
		shunt_capacitance,
		shunt_conductance,
	)
end

const _REQ_SEMICON = (:r_in, :r_ex, :material_props)
const _OPT_SEMICON = (:temperature,)
const _DEFS_SEMICON = (T₀,)

Validation.has_radii(::Type{Semicon}) = true
Validation.has_temperature(::Type{Semicon}) = true
Validation.required_fields(::Type{Semicon}) = _REQ_SEMICON
Validation.keyword_fields(::Type{Semicon}) = _OPT_SEMICON
Validation.keyword_defaults(::Type{Semicon}) = _DEFS_SEMICON

# accept proxies for radii
Validation.is_radius_input(::Type{Semicon}, ::Val{:r_in}, x::AbstractCablePart) = true
Validation.is_radius_input(::Type{Semicon}, ::Val{:r_in}, x::Thickness) = true
Validation.is_radius_input(::Type{Semicon}, ::Val{:r_ex}, x::Thickness) = true
Validation.is_radius_input(::Type{Semicon}, ::Val{:r_ex}, x::Diameter) = true

Validation.extra_rules(::Type{Semicon}) = (IsA{Material}(:material_props),)

# normalize proxies -> numbers
Validation.parse(::Type{Semicon}, nt) = begin
	rin, rex = _normalize_radii(Semicon, nt.r_in, nt.r_ex)
	(; nt..., r_in = rin, r_ex = rex)
end

# This macro expands to a weakly-typed constructor for Semicon
@construct Semicon _REQ_SEMICON _OPT_SEMICON _DEFS_SEMICON
