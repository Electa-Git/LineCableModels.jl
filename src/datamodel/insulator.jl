"""
$(TYPEDEF)

Represents an insulating layer with defined geometric, material, and electrical properties given by the attributes:

$(TYPEDFIELDS)
"""
struct Insulator{T <: REALSCALAR} <: AbstractInsulatorPart{T}
	"Internal radius of the insulating layer \\[m\\]."
	r_in::T
	"External radius of the insulating layer \\[m\\]."
	r_ex::T
	"Material properties of the insulator."
	material_props::Material{T}
	"Operating temperature of the insulator \\[°C\\]."
	temperature::T
	"Cross-sectional area of the insulating layer \\[m²\\]."
	cross_section::T
	"Electrical resistance of the insulating layer \\[Ω/m\\]."
	resistance::T
	"Geometric mean radius of the insulator \\[m\\]."
	gmr::T
	"Shunt capacitance per unit length of the insulating layer \\[F/m\\]."
	shunt_capacitance::T
	"Shunt conductance per unit length of the insulating layer \\[S·m\\]."
	shunt_conductance::T
end

"""
$(TYPEDSIGNATURES)

Constructs an [`Insulator`](@ref) object with specified geometric and material parameters.

# Arguments

- `r_in`: Internal radius of the insulating layer \\[m\\].
- `r_ex`: External radius or thickness of the layer \\[m\\].
- `material_props`: Material properties of the insulating material.
- `temperature`: Operating temperature of the insulator \\[°C\\].

# Returns

- An [`Insulator`](@ref) object with calculated electrical properties.

# Examples

```julia
material_props = Material(1e10, 3.0, 1.0, 20.0, 0.0)
insulator_layer = $(FUNCTIONNAME)(0.01, 0.015, material_props, temperature=25)
```
"""
function Insulator(
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
	return Insulator(
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

const _REQ_INSULATOR = (:r_in, :r_ex, :material_props)
const _OPT_INSULATOR = (:temperature,)
const _DEFS_INSULATOR = (T₀,)

Validation.has_radii(::Type{Insulator}) = true
Validation.has_temperature(::Type{Insulator}) = true
Validation.required_fields(::Type{Insulator}) = _REQ_INSULATOR
Validation.keyword_fields(::Type{Insulator}) = _OPT_INSULATOR
Validation.keyword_defaults(::Type{Insulator}) = _DEFS_INSULATOR

# accept proxies for radii
Validation.is_radius_input(::Type{Insulator}, ::Val{:r_in}, x::AbstractCablePart) =
	true
Validation.is_radius_input(::Type{Insulator}, ::Val{:r_in}, x::Thickness) = true
Validation.is_radius_input(::Type{Insulator}, ::Val{:r_ex}, x::Thickness) = true
Validation.is_radius_input(::Type{Insulator}, ::Val{:r_ex}, x::Diameter) = true

Validation.extra_rules(::Type{Insulator}) = (IsA{Material}(:material_props),)

# normalize proxies -> numbers
Validation.parse(::Type{Insulator}, nt) = begin
	rin, rex = _normalize_radii(Insulator, nt.r_in, nt.r_ex)
	(; nt..., r_in = rin, r_ex = rex)
end

# This macro expands to a weakly-typed constructor for Insulator
@construct Insulator _REQ_INSULATOR _OPT_INSULATOR _DEFS_INSULATOR
