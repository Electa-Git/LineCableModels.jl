defaults(::Type{Formula{:lossless}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the admittivity of an insulator layer of `material` at the frequency `f`, without
conduction or polarization loss:

```math
\\kappa = \\mathrm{j} 2 \\pi f \\varepsilon_0 \\varepsilon_r,
```

Here ``f`` is frequency, ``\\kappa`` is admittivity, and ``\\varepsilon_r`` is the relative
permittivity `eps_r` of the material \\[dimensionless\\]. The constant
``\\varepsilon_0`` is [`vacuum_permittivity`](@ref).

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
