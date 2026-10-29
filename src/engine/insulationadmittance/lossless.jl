defaults(::Type{Formula{:lossless}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the admittivity of an insulator layer of `material` at the frequency `f`, by the
formula `Formula{:lossless}`, whose equation is in the section "Formulas" of
[`InsulationAdmittance`](@ref).

# Arguments

- `Formula{:lossless}`: dielectric law without conduction or polarization loss.
- `material`: material of the layer.
- `f`: frequency \\[Hz\\].

# Returns

- Admittivity \\[S/m\\].
"""
function admittivity(::Formula{:lossless}, material::Material{T}, f::T) where {T <: Real}
    return im * 2π * f * vacuum_permittivity * material.eps_r
end
