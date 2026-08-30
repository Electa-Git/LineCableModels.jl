assumptions(::Val{:Scott1967}) = (;)

"""
$(TYPEDSIGNATURES)

**Identification.** Empirical moist-rock relation over the measured range
100 Hz to 1 MHz.

**Expression.** Let ``s=\\log_{10}(1000/\\rho_0)`` for the 100 Hz conductivity
in mS/m and ``x=\\log_{10}f``. Then

```math
\\log_{10}\\varepsilon_r=5.491+0.946s-1.097x+0.069s^2-0.114sx+0.067x^2,
```

```math
\\log_{10}\\widehat\\sigma=0.028+1.098s-0.068x+0.036s^2-0.046sx+0.018x^2,
\\qquad \\sigma=10^{-3}\\widehat\\sigma.
```

**Reference.** J. H. Scott, R. D. Carroll, and D. R. Cunningham, “Dielectric
Constant and Electrical Conductivity Measurements of Moist Rock: A New
Laboratory Method,” *Journal of Geophysical Research*, 72(20), 1967.
DOI: 10.1029/JZ072i020p05101.
"""
description(::Formula{:Scott1967}) =
    "Scott–Carroll–Cunningham 1967 (empirical moist-soil fit)"

#=
Evaluate the Scott–Carroll–Cunningham empirical soil relation.

The polynomial fit uses conductivity at 100 Hz in mS/m and frequency in Hz.
The returned resistivity is converted to Ω·m and relative permittivity is
dimensionless. The source data cover 100 Hz to 1 MHz.

# Reference

J. H. Scott, R. D. Carroll, and D. R. Cunningham, *Dielectric Constant and
Electrical Conductivity Measurements of Moist Rock: A New Laboratory Method*,
Journal of Geophysical Research, 72(20), 1967. DOI: 10.1029/JZ072i020p05101.
=#
function (::Functor{:Scott1967})(
        material::EarthMaterial{T},
        frequency::T,
        ::NamedTuple
) where {T <: Real}
    thousand = convert(T, 1000)
    conductivity_100hz = thousand / material.rho
    conductivity_log = log10(conductivity_100hz)
    frequency_log = log10(frequency)
    permittivity_log = convert(T, 5.491) +
                       convert(T, 0.946) * conductivity_log -
                       convert(T, 1.097) * frequency_log +
                       convert(T, 0.069) * conductivity_log^2 -
                       convert(T, 0.114) * conductivity_log * frequency_log +
                       convert(T, 0.067) * frequency_log^2
    output_conductivity_log = convert(T, 0.028) +
                              convert(T, 1.098) * conductivity_log -
                              convert(T, 0.068) * frequency_log +
                              convert(T, 0.036) * conductivity_log^2 -
                              convert(T, 0.046) * conductivity_log * frequency_log +
                              convert(T, 0.018) * frequency_log^2
    ten = convert(T, 10)
    relative_permittivity = ten^permittivity_log
    conductivity = ten^output_conductivity_log
    return EarthMaterial{T}(
        thousand / conductivity,
        relative_permittivity,
        material.mu_r
    )
end

# Return the stable discovery identifier.
:Scott1967
