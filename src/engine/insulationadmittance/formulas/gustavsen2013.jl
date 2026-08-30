assumptions(::Val{:Gustavsen2013}) = (;)

function description(::Formula{:Gustavsen2013})
    "Gustavsen conventional lossless-insulation model (2013)"
end

"""
$(TYPEDSIGNATURES)

Apply the conventional lossless-insulation approximation by suppressing
material conductivity while retaining permittivity:

```math
\\kappa=j\\omega\\varepsilon,
\\qquad
Y=j\\omega C.
```

# Arguments

- `material`: Static insulation properties.
- `frequency`: Evaluation frequency \\[Hz\\].
- `temperature`: Operating temperature \\[°C\\].
- `values`: Formula assumptions.

# Returns

- An insulation material with infinite electrical resistivity and all other
  properties unchanged.

# References

B. Gustavsen, H. K. Høidalen, and T. M. Ohnstad, “Field Measurement and
Simulation of 132 kV Oil-Filled Submarine Cables,” IPST, 2013.
"""
@inline function (::Functor{:Gustavsen2013})(
        material::Material{T},
        frequency::Real,
        temperature::Real,
        values::NamedTuple
) where {T <: Real}
    return Material{T}(
        material.kind,
        convert(T, Inf),
        material.eps_r,
        material.mu_r,
        material.T0,
        material.alpha,
        material.rho_thermal,
        material.theta_max,
        material.tan_delta,
        material.sigma_solar
    )
end

:Gustavsen2013
