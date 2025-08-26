"""
$(TYPEDEF)

Represents an array of wires equally spaced around a circumference of arbitrary radius, with attributes:

$(TYPEDFIELDS)
"""
struct WireArray{T<:REALSCALAR,U<:Int} <: AbstractWireArray
    "Internal radius of the wire array \\[m\\]."
    radius_in::T
    "External radius of the wire array \\[m\\]."
    radius_ext::T
    "Radius of each individual wire \\[m\\]."
    radius_wire::T
    "Number of wires in the array \\[dimensionless\\]."
    num_wires::U
    "Ratio defining the lay length of the wires (twisting factor) \\[dimensionless\\]."
    lay_ratio::T
    "Mean diameter of the wire array \\[m\\]."
    mean_diameter::T
    "Pitch length of the wire array \\[m\\]."
    pitch_length::T
    "Twisting direction of the strands (1 = unilay, -1 = contralay) \\[dimensionless\\]."
    lay_direction::U
    "Material object representing the physical properties of the wire material."
    material_props::Material{T}
    "Temperature at which the properties are evaluated \\[°C\\]."
    temperature::T
    "Cross-sectional area of all wires in the array \\[m²\\]."
    cross_section::T
    "Electrical resistance per wire in the array \\[Ω/m\\]."
    resistance::T
    "Geometric mean radius of the wire array \\[m\\]."
    gmr::T
end

"""
$(TYPEDSIGNATURES)

Constructs a [`WireArray`](@ref) instance based on specified geometric and material parameters.

# Arguments

- `radius_in`: Internal radius of the wire array \\[m\\].
- `radius_wire`: Radius of each individual wire \\[m\\].
- `num_wires`: Number of wires in the array \\[dimensionless\\].
- `lay_ratio`: Ratio defining the lay length of the wires (twisting factor) \\[dimensionless\\].
- `material_props`: A [`Material`](@ref) object representing the material properties.
- `temperature`: Temperature at which the properties are evaluated \\[°C\\].
- `lay_direction`: Twisting direction of the strands (1 = unilay, -1 = contralay) \\[dimensionless\\].

# Returns

- A [`WireArray`](@ref) object with calculated geometric and electrical properties.

# Examples

```julia
material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
wire_array = $(FUNCTIONNAME)(0.01, Diameter(0.002), 7, 10, material_props, temperature=25)
println(wire_array.mean_diameter)  # Outputs mean diameter in m
println(wire_array.resistance)     # Outputs resistance in Ω/m
```

# See also

- [`Material`](@ref)
- [`ConductorGroup`](@ref)
- [`calc_tubular_resistance`](@ref)
- [`calc_wirearray_gmr`](@ref)
- [`calc_helical_params`](@ref)
"""
function WireArray(
    radius_in::T,
    radius_wire::T,
    num_wires::U,
    lay_ratio::T,
    material_props::Material{T},
    temperature::T,
    lay_direction::U,
) where {T<:REALSCALAR,U<:Int}

    rho = material_props.rho
    T0 = material_props.T0
    alpha = material_props.alpha
    radius_ext = num_wires == 1 ? radius_wire : radius_in + 2 * radius_wire

    mean_diameter, pitch_length, overlength = calc_helical_params(
        radius_in,
        radius_ext,
        lay_ratio,
    )

    cross_section = num_wires * (π * radius_wire^2)

    R_wire =
        calc_tubular_resistance(0.0, radius_wire, rho, alpha, T0, temperature) *
        overlength
    R_all_wires = R_wire / num_wires

    gmr = calc_wirearray_gmr(
        radius_in + radius_wire,
        num_wires,
        radius_wire,
        material_props.mu_r,
    )

    # Initialize object
    return WireArray(
        radius_in,
        radius_ext,
        radius_wire,
        num_wires,
        lay_ratio,
        mean_diameter,
        pitch_length,
        lay_direction,
        material_props,
        temperature,
        cross_section,
        R_all_wires,
        gmr,
    )
end

const _REQ_WIREARRAY = (:radius_in, :radius_wire, :num_wires, :lay_ratio, :material_props,)
const _OPT_WIREARRAY = (:temperature, :lay_direction,)
const _DEFS_WIREARRAY = (T₀, 1,)

Validation.has_radii(::Type{WireArray}) = false
Validation.has_temperature(::Type{WireArray}) = true
Validation.required_fields(::Type{WireArray}) = _REQ_WIREARRAY
Validation.keyword_fields(::Type{WireArray}) = _OPT_WIREARRAY
Validation.coercive_fields(::Type{WireArray}) = (:radius_in, :radius_wire, :lay_ratio, :material_props, :temperature)  # not :num_wires, :lay_direction
# accept proxies for radii
Validation.is_radius_input(::Type{WireArray}, ::Val{:radius_in},
    x::AbstractCablePart) = true
Validation.is_radius_input(::Type{WireArray}, ::Val{:radius_ext},
    x::Diameter) = true

Validation.extra_rules(::Type{WireArray}) = (
    # radii (post-parse they must be numeric)
    Normalized(:radius_in), Finite(:radius_in), Nonneg(:radius_in),
    Normalized(:radius_wire), Finite(:radius_wire), Positive(:radius_wire),

    # counts and geometry params
    IntegerField(:num_wires), Positive(:num_wires),
    Finite(:lay_ratio), Nonneg(:lay_ratio),

    # material type
    IsA{Material}(:material_props),

    # lay direction constraint (pin to -1 or +1)
    OneOf(:lay_direction, (-1, 1)),
)

# normalize proxies -> numbers
Validation.parse(::Type{WireArray}, nt) = begin
    rin, rw = _normalize_radii(WireArray, nt.radius_in, nt.radius_wire)
    (; nt..., radius_in=rin, radius_wire=rw)
end

# This macro expands to a weakly-typed constructor for WireArray
@_ctor WireArray _REQ_WIREARRAY _OPT_WIREARRAY _DEFS_WIREARRAY

