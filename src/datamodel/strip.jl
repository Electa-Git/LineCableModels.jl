"""
$(TYPEDEF)

Represents a flat conductive strip with defined geometric and material properties given by the attributes:

$(TYPEDFIELDS)
"""
struct Strip{T<:REALSCALAR} <: AbstractConductorPart{T}
    "Internal radius of the strip \\[m\\]."
    radius_in::T
    "External radius of the strip \\[m\\]."
    radius_ext::T
    "Thickness of the strip \\[m\\]."
    thickness::T
    "Width of the strip \\[m\\]."
    width::T
    "Ratio defining the lay length of the strip (twisting factor) \\[dimensionless\\]."
    lay_ratio::T
    "Mean diameter of the strip's helical path \\[m\\]."
    mean_diameter::T
    "Pitch length of the strip's helical path \\[m\\]."
    pitch_length::T
    "Twisting direction of the strip (1 = unilay, -1 = contralay) \\[dimensionless\\]."
    lay_direction::Int
    "Material properties of the strip."
    material_props::Material{T}
    "Temperature at which the properties are evaluated \\[°C\\]."
    temperature::T
    "Cross-sectional area of the strip \\[m²\\]."
    cross_section::T
    "Electrical resistance of the strip \\[Ω/m\\]."
    resistance::T
    "Geometric mean radius of the strip \\[m\\]."
    gmr::T
end

"""
$(TYPEDSIGNATURES)

Constructs a [`Strip`](@ref) object with specified geometric and material parameters.

# Arguments

- `radius_in`: Internal radius of the strip \\[m\\].
- `radius_ext`: External radius or thickness of the strip \\[m\\].
- `width`: Width of the strip \\[m\\].
- `lay_ratio`: Ratio defining the lay length of the strip \\[dimensionless\\].
- `material_props`: Material properties of the strip.
- `temperature`: Temperature at which the properties are evaluated \\[°C\\]. Defaults to [`T₀`](@ref).
- `lay_direction`: Twisting direction of the strip (1 = unilay, -1 = contralay) \\[dimensionless\\]. Defaults to 1.

# Returns

- A [`Strip`](@ref) object with calculated geometric and electrical properties.

# Examples

```julia
material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
strip = $(FUNCTIONNAME)(0.01, Thickness(0.002), 0.05, 10, material_props, temperature=25)
println(strip.cross_section) # Output: 0.0001 [m²]
println(strip.resistance)    # Output: Resistance value [Ω/m]
```

# See also

- [`Material`](@ref)
- [`ConductorGroup`](@ref)
- [`calc_strip_resistance`](@ref)
- [`calc_tubular_gmr`](@ref)
- [`calc_helical_params`](@ref)
"""
function Strip(
    radius_in::T,
    radius_ext::T,
    width::T,
    lay_ratio::T,
    material_props::Material{T},
    temperature::T,
    lay_direction::Int,
) where {T<:REALSCALAR}

    thickness = radius_ext - radius_in
    rho = material_props.rho
    T0 = material_props.T0
    alpha = material_props.alpha

    mean_diameter, pitch_length, overlength = calc_helical_params(
        radius_in,
        radius_ext,
        lay_ratio,
    )

    cross_section = thickness * width

    R_strip =
        calc_strip_resistance(thickness, width, rho, alpha, T0, temperature) *
        overlength

    gmr = calc_tubular_gmr(radius_ext, radius_in, material_props.mu_r)

    # Initialize object
    return Strip(
        radius_in,
        radius_ext,
        thickness,
        width,
        lay_ratio,
        mean_diameter,
        pitch_length,
        lay_direction,
        material_props,
        temperature,
        cross_section,
        R_strip,
        gmr,
    )
end

const _REQ_STRIP = (:radius_in, :radius_ext, :width, :lay_ratio, :material_props,)
const _OPT_STRIP = (:temperature, :lay_direction,)
const _DEFS_STRIP = (T₀, 1,)

Validation.has_radii(::Type{Strip}) = true
Validation.has_temperature(::Type{Strip}) = true
Validation.required_fields(::Type{Strip}) = _REQ_STRIP
Validation.keyword_fields(::Type{Strip}) = _OPT_STRIP
Validation.keyword_defaults(::Type{Strip}) = _DEFS_STRIP

Validation.coercive_fields(::Type{Strip}) = (:radius_in, :radius_ext, :width, :lay_ratio, :material_props, :temperature)  # not :lay_direction
# accept proxies for radii

Validation.is_radius_input(::Type{Strip}, ::Val{:radius_in}, x::AbstractCablePart) = true
Validation.is_radius_input(::Type{Strip}, ::Val{:radius_in}, x::Thickness) = true
Validation.is_radius_input(::Type{Strip}, ::Val{:radius_ext}, x::Thickness) = true
Validation.is_radius_input(::Type{Strip}, ::Val{:radius_ext}, x::Diameter) = true

Validation.extra_rules(::Type{Strip}) = (IsA{Material}(:material_props), OneOf(:lay_direction, (-1, 1)), Finite(:lay_ratio), Nonneg(:lay_ratio), Finite(:width), Positive(:width),)

# normalize proxies -> numbers
Validation.parse(::Type{Strip}, nt) = begin
    rin, rex = _normalize_radii(Strip, nt.radius_in, nt.radius_ext)
    (; nt..., radius_in=rin, radius_ext=rex)
end

# This macro expands to a weakly-typed constructor for Strip
@_ctor Strip _REQ_STRIP _OPT_STRIP _DEFS_STRIP
