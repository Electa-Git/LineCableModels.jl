assumptions(::Val{:Marti2001}) = (;)

description(::Formula{:Marti2001}) = "Yu–Martí parallel-conductance insulation model (2001)"

"""
$(TYPEDSIGNATURES)

Retain the material conductivity and permittivity used by the Yu–Marti
parallel-conductance insulation model:

```math
\\kappa=\\sigma+j\\omega\\varepsilon,
\\qquad
Y=G+j\\omega C.
```

The common Coaxial Engine operator applies annular geometry to ``\\kappa`` and
combines physical dielectric layers radially in series.

# Arguments

- `material`: Static insulation properties.
- `frequency`: Evaluation frequency \\[Hz\\].
- `temperature`: Operating temperature \\[°C\\].
- `values`: Formula assumptions.

# Returns

- The unchanged material, retaining both finite resistivity and permittivity.

# References

T. C. Yu and J. R. Martí, “zCable Model for Frequency Dependent Modelling of
Cable Transmission Systems,” IPST, 2001.
"""
@inline function (::Functor{:Marti2001})(
        material::Material,
        frequency::Real,
        temperature::Real,
        values::NamedTuple
)
    return material
end

:Marti2001
