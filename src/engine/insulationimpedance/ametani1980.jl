defaults(::Type{Formula{:ametani1980}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the magnetic term of the series impedance across an annular dielectric layer:

```math
z = \\frac{s \\mu_0 \\mu_r}{2 \\pi} \\ln \\frac{r_\\mathrm{ex}}{r_\\mathrm{in}},
```

Here ``s`` is complex frequency, ``r_\\mathrm{in}`` and ``r_\\mathrm{ex}`` are the layer
radii, and ``\\mu_r`` is relative permeability. The constant ``\\mu_0`` is
[`vacuum_permeability`](@ref). The term is zero when `r_in` is zero or equals `r_ex`.

# Arguments

- `Formula{:ametani1980}`: magnetic annular-dielectric formula.
- `r_in`: inner radius \\[m\\].
- `r_ex`: outer radius, in \\[m\\].
- `mu_r`: relative permeability \\[dimensionless\\].
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].

# Returns

- Impedance \\[Ω/m\\].
"""
function impedance(::Formula{:ametani1980}, r_in::T, r_ex::T, mu_r::T,
        s::Complex{T}) where {T <: Real}
    factor = s * vacuum_permeability * mu_r / 2π
    if iszero(r_in) || r_in == r_ex
        return zero(factor)
    end
    return factor * log(r_ex / r_in)
end
