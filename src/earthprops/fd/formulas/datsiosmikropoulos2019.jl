"Return the dry-soil assumption of the Datsios–Mikropoulos relation."
assumptions(::Val{:DatsiosMikropoulos2019}) = (dry_permittivity = 3.5,)

description(::Formula{:DatsiosMikropoulos2019}) =
    "Datsios–Mikropoulos 2019 (two-limit soil fit)"

#=
Evaluate the Datsios–Mikropoulos soil relation.

The published conductivity fit uses μS/cm and is converted to S/m before the
result is returned as resistivity in Ω·m. Relative permittivity is held at its
3 kHz value below the fitted 3 kHz boundary.

# Reference

Z. G. Datsios and P. N. Mikropoulos, *Characterization of the Frequency
Dependence of the Electrical Properties of Sandy Soil With Variable Grain Size
and Water Content*, IEEE Transactions on Dielectrics and Electrical Insulation,
26(3), 2019. DOI: 10.1109/TDEI.2018.007864.
=#
function (::Functor{:DatsiosMikropoulos2019})(
        material::EarthMaterial{T},
        frequency::T,
        values::NamedTuple
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
        one(frequency) + convert(T, 0.65) /
        conductivity_low^convert(T, 0.57)
    )
    micro = convert(T, 1e-6)
    conductivity_micro_siemens_per_centimetre =
        conductivity_high * frequency * micro +
        (frequency_low - frequency_low * (frequency - frequency_low) * micro) *
        (conductivity_low / frequency_low - conductivity_high * micro)
    conductivity = convert(T, 1e-4) *
                   conductivity_micro_siemens_per_centimetre
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

# Return the stable discovery identifier.
:DatsiosMikropoulos2019
