assumptions(::Val{:Gustavsen2013}) = (;)

"""
$(TYPEDSIGNATURES)

**Identification.** Conventional lossless-insulation approximation. The
material permittivity is retained and its finite conductivity is suppressed.

**Expression.**

```math
\\kappa=j\\omega\\varepsilon_0\\varepsilon_r,\\qquad
G=0,\\qquad
C=\\frac{2\\pi\\varepsilon_0\\varepsilon_r}{\\ln(b/a)}.
```

**Reference.** B. Gustavsen, H. K. Høidalen, and T. M. Ohnstad, “Field
Measurement and Simulation of 132 kV Oil-Filled Submarine Cables,”
*International Conference on Power Systems Transients*, 2013.
"""
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
