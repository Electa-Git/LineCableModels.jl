"""
$(TYPEDEF)

Represents a tubular or solid (`r_in=0`) conductor with geometric and material properties defined as:

$(TYPEDFIELDS)
"""
struct Tubular{T <: REALSCALAR} <: AbstractConductorPart{T}
	"Internal radius of the tubular conductor \\[m\\]."
	r_in::T
	"External radius of the tubular conductor \\[m\\]."
	r_ex::T
	"A [`Material`](@ref) object representing the physical properties of the conductor material."
	material_props::Material{T}
	"Temperature at which the properties are evaluated \\[°C\\]."
	temperature::T
	"Cross-sectional area of the tubular conductor \\[m²\\]."
	cross_section::T
	"Electrical resistance (DC) of the tubular conductor \\[Ω/m\\]."
	resistance::T
	"Geometric mean radius of the tubular conductor \\[m\\]."
	gmr::T
end

"""
$(TYPEDSIGNATURES)

Initializes a [`Tubular`](@ref) object with specified geometric and material parameters.

# Arguments

- `r_in`: Internal radius of the tubular conductor \\[m\\].
- `r_ex`: External radius of the tubular conductor \\[m\\].
- `material_props`: A [`Material`](@ref) object representing the physical properties of the conductor material.
- `temperature`: Temperature at which the properties are evaluated \\[°C\\]. Defaults to [`T₀`](@ref).

# Returns

- An instance of [`Tubular`](@ref) initialized with calculated geometric and electrical properties.

# Examples

```julia
material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
tubular = $(FUNCTIONNAME)(0.01, 0.02, material_props, temperature=25)
println(tubular.cross_section) # Output: 0.000942 [m²]
println(tubular.resistance)    # Output: Resistance value [Ω/m]
```

# See also

- [`Material`](@ref)
- [`calc_tubular_resistance`](@ref)
- [`calc_tubular_gmr`](@ref)
"""
function Tubular(
	r_in::T,
	r_ex::T,
	material_props::Material{T},
	temperature::T,
) where {T <: REALSCALAR}

	rho = material_props.rho
	T0 = material_props.T0
	alpha = material_props.alpha

	cross_section = π * (r_ex^2 - r_in^2)

	R0 = calc_tubular_resistance(r_in, r_ex, rho, alpha, T0, temperature)

	gmr = calc_tubular_gmr(r_ex, r_in, material_props.mu_r)

	# Initialize object
	return Tubular(
		r_in,
		r_ex,
		material_props,
		temperature,
		cross_section,
		R0,
		gmr,
	)
end

const _REQ_TUBULAR = (:r_in, :r_ex, :material_props)
const _OPT_TUBULAR = (:temperature,)
const _DEFS_TUBULAR = (T₀,)

Validation.has_radii(::Type{Tubular}) = true
Validation.has_temperature(::Type{Tubular}) = true
Validation.required_fields(::Type{Tubular}) = _REQ_TUBULAR
Validation.keyword_fields(::Type{Tubular}) = _OPT_TUBULAR
Validation.keyword_defaults(::Type{Tubular}) = _DEFS_TUBULAR

# accept proxies for radii
Validation.is_radius_input(::Type{Tubular}, ::Val{:r_in}, x::AbstractCablePart) = true
Validation.is_radius_input(::Type{Tubular}, ::Val{:r_in}, x::Thickness) = true
Validation.is_radius_input(::Type{Tubular}, ::Val{:r_ex}, x::Thickness) = true
Validation.is_radius_input(::Type{Tubular}, ::Val{:r_ex}, x::Diameter) = true

Validation.extra_rules(::Type{Tubular}) = (IsA{Material}(:material_props),)

# normalize proxies -> numbers
Validation.parse(::Type{Tubular}, nt) = begin
	rin, rex = _normalize_radii(Tubular, nt.r_in, nt.r_ex)
	(; nt..., r_in = rin, r_ex = rex)
end

# This macro expands to a weakly-typed constructor for Tubular
@construct Tubular _REQ_TUBULAR _OPT_TUBULAR _DEFS_TUBULAR


