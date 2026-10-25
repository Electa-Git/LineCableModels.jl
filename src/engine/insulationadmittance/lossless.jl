defaults(::Type{Formula{:lossless}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the admittivity of an insulating layer of `material` at the frequency `f`, without
conduction or polarization loss:

```math
\\kappa = \\mathrm{j} 2 \\pi f \\varepsilon_0 \\varepsilon_r,
```

where ``\\varepsilon_0`` is the permittivity of vacuum and ``\\varepsilon_r`` the relative
permittivity `eps_r` of the material.

# Arguments

- `material`: material of the layer.
- `f`: frequency \\[Hz\\].

# Returns

- Admittivity \\[S/m\\].
"""
function admittivity(::Formula{:lossless}, material::Material{T}, f::T) where {T <: Real}
    return im * 2π * f * vacuum_permittivity * material.eps_r
end
