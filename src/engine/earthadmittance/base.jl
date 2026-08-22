"""
$(TYPEDEF)

Reference each cable's outermost insulation directly to an ideal ground plane.

The environmental potential coefficient is

```math
P_{e,ij}(s) = 0 \\quad \\mathrm{m/F}.
```

The cable's outermost insulation remains part of the assembled potential-coefficient
matrix. Inter-cable shunt coupling through soil is omitted. This formulation matches
line-constant models that calculate cable insulation admittance against an ideal ground
reference while applying an independent earth-return formulation to series impedance.
"""
struct IdealGround <: EarthAdmittanceFormulation end

description(::IdealGround) = "Ideal ground reference"

@inline function (::IdealGround)(
        ::Val{:self},
        heights,
        separation,
        earth_resistivity,
        earth_permittivity,
        earth_permeability,
        frequency_point::Complex{T}
) where {T <: Real}
    return zero(frequency_point)
end

@inline function (::IdealGround)(
        ::Val{:mutual},
        heights,
        separation,
        earth_resistivity,
        earth_permittivity,
        earth_permeability,
        frequency_point::Complex{T}
) where {T <: Real}
    return zero(frequency_point)
end

@inline function Base.getproperty(f::Homogeneous, name::Symbol)
    if name === :s || name === :t || name === :Γx || name === :γ1 || name === :γ2 ||
       name === :μ2
        return getproperty(from_kernel(f), name)
    end
    return getfield(f, name)  # subtype-specific fields (if any)
end
