"Return the 100 Hz boundary assumption of the Visacro–Alipio empirical relation."
assumptions(::Val{:VisacroAlipio2012}) = (frequency_boundary = 100.0,)

"""
$(TYPEDSIGNATURES)

**Identification.** Empirical causal fit with a 100 Hz lower boundary.

**Expression.** With ``f_e=\\max(f,100)`` and ``\\sigma_0=1/\\rho_0``,

```math
\\varepsilon_r(f)=1.3+7.6\\times10^3 f_e^{-0.4},
\\qquad
\\sigma(f)=\\sigma_0+1.2\\times10^{-6}\\sigma_0^{0.27}(f_e-100)^{0.65}.
```

**Reference.** S. Visacro and R. Alipio, “Frequency Dependence of Soil
Parameters: Experimental Results, Predicting Formula and Influence on the
Lightning Response of Grounding Electrodes,” *IEEE Transactions on Power
Delivery*, 27(2), 2012. DOI: 10.1109/TPWRD.2011.2179070.
"""
description(::Formula{:VisacroAlipio2012}) =
    "Visacro–Alipio empirical soil dispersion (2012)"

#=
Evaluate the Visacro–Alipio empirical soil relation.

The fit is real-valued from 100 Hz upward. Below that boundary the 100 Hz
properties are retained instead of relying on MATLAB complex promotion followed
by `abs`. Input and returned resistivity use Ω·m.

# Reference

S. Visacro and R. Alipio, *Frequency Dependence of Soil Parameters:
Experimental Results, Predicting Formula and Influence on the Lightning
Response of Grounding Electrodes*, IEEE Transactions on Power Delivery,
27(2), 2012. DOI: 10.1109/TPWRD.2011.2179070.
=#
function earth_material(
        ::Val{:VisacroAlipio2012},
        material::EarthMaterial{T},
        frequency::T,
        values::NamedTuple
) where {T <: Real}
    frequency_boundary = convert(T, values.frequency_boundary)
    evaluated_frequency = frequency < frequency_boundary ?
                          frequency_boundary : frequency
    conductivity_reference = inv(material.rho)
    relative_permittivity = convert(T, 1.3) +
                            convert(T, 7.6e3) *
                            evaluated_frequency^convert(T, -0.4)
    conductivity = conductivity_reference +
                   convert(T, 1.2e-6) *
                   conductivity_reference^convert(T, 0.27) *
                   (evaluated_frequency - frequency_boundary)^convert(T, 0.65)
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

# Return the stable discovery identifier.
:VisacroAlipio2012
