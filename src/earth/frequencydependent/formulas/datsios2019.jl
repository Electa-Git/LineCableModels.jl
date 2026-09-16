"Return the dry-soil parameter of the Datsios–Mikropoulos relation."
assumptions(::Val{:datsios2019}) = (dry_permittivity = 3.5,)

"""
$(TYPEDSIGNATURES)

Datsios–Mikropoulos two-limit fit for sandy soil, with relative permittivity
held at its 3 kHz value below the fitted boundary.

**Expression.** With ``\\widehat\\sigma_{42}=10^4/\\rho_0`` in μS/cm and
dry permittivity ``\\varepsilon_d=3.5``,

```math
p=0.537\\widehat\\sigma_{42}^{0.16},\\quad
\\varepsilon_\\infty=1.24\\widehat\\sigma_{42}^{0.415}\\varepsilon_d,
\\quad
\\varepsilon_{3k}=4\\widehat\\sigma_{42}^{0.463}(2.9\\varepsilon_d-3.8).
```

The conductivity interpolation is evaluated between the measured 42 Hz and
high-frequency limits.

**Reference.** Z. G. Datsios and P. N. Mikropoulos, *IEEE Transactions on
Dielectrics and Electrical Insulation*, 26(3), 2019.
"""
function description(::Type{<:Formula{:datsios2019}}; compact::Bool=false)
    compact ? "Datsios" : "Datsios–Mikropoulos two-limit soil fit (2019)"
end

function earth_material(
        ::Formula{:datsios2019}, material::EarthMaterial{T}, frequency::T,
        values::NamedTuple, options::FormulationOptions, workspace
) where {T <: Real}
    conductivity_low = convert(T, 1e4) / material.rho
    frequency_low = convert(T, 42)
    frequency_boundary = convert(T, 3000)
    dry_permittivity = convert(T, values.dry_permittivity)

    permittivity_exponent = convert(T, 0.537) * conductivity_low^convert(T, 0.16)
    dry_permittivity_3khz = convert(T, 2.9) * dry_permittivity - convert(T, 3.8)
    permittivity_high = convert(T, 1.24) *
                        conductivity_low^convert(T, 0.415) * dry_permittivity
    permittivity_3khz = convert(T, 4) *
                        conductivity_low^convert(T, 0.463) * dry_permittivity_3khz
    permittivity_frequency = frequency < frequency_boundary ?
                             frequency_boundary : frequency
    relative_permittivity = permittivity_high +
                            (frequency_boundary / permittivity_frequency)^permittivity_exponent *
                            (permittivity_3khz - permittivity_high)

    conductivity_high = conductivity_low * (
        one(frequency) + convert(T, 0.65) / conductivity_low^convert(T, 0.57)
    )
    micro = convert(T, 1e-6)
    conductivity_micro_siemens_per_centimetre =
        conductivity_high * frequency * micro +
        (frequency_low - frequency_low * (frequency - frequency_low) * micro) *
        (conductivity_low / frequency_low - conductivity_high * micro)
    conductivity = convert(T, 1e-4) * conductivity_micro_siemens_per_centimetre
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

formulation_options(::FormulaMethod{<:Formula{:datsios2019}, typeof(earth_material)}) = FormulationOptions()

:datsios2019
