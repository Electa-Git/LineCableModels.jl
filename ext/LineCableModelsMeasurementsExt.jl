"""
    LineCableModelsMeasurementsExt

Compute linear propagation of uncertainty. A `Measurement` input gives series impedances
`Z` \\[Ω/m\\] and shunt admittances `Y` \\[S/m\\] with uncertain real and imaginary parts.
Define the `nominal` method and four methods of SpecialFunctions Bessel functions on an
uncertain complex argument. The engine calls these functions for Schelkunoff surface
impedances and the kernels of `:martins2026`.

For ``z = x + iy``, write the scaled functions as
``I^s_\\nu(z) = \\mathrm{e}^{-|x|} I_\\nu(z)`` (`besselix`),
``K^s_\\nu(z) = \\mathrm{e}^{z} K_\\nu(z)`` (`besselkx`), and
``J^s_\\nu(z) = \\mathrm{e}^{-|y|} J_\\nu(z)`` (`besseljx`), with ``i^2 = -1``.
The unscaled ``K_\\nu(z)`` is `besselk`.
The Cartesian partial derivatives and their order recurrences are

```math
\\begin{aligned}
\\partial_x I^s_\\nu &= (I^s_{\\nu-1}+I^s_{\\nu+1})/2 - \\mathrm{sgn}(x)I^s_\\nu, \\\\
\\partial_y I^s_\\nu &= i(I^s_{\\nu-1}+I^s_{\\nu+1})/2, \\\\
\\partial_x K^s_\\nu &= K^s_\\nu - (K^s_{\\nu-1}+K^s_{\\nu+1})/2, &
\\partial_y K^s_\\nu &= i\\partial_x K^s_\\nu, \\\\
\\partial_x J^s_\\nu &= (J^s_{\\nu-1}-J^s_{\\nu+1})/2, &
\\partial_y J^s_\\nu &= i\\partial_x J^s_\\nu - \\mathrm{sgn}(y)J^s_\\nu, \\\\
\\partial_x K_\\nu &= -(K_{\\nu-1}+K_{\\nu+1})/2, &
\\partial_y K_\\nu &= i\\partial_x K_\\nu, \\\\
2I_\\nu' &= I_{\\nu-1}+I_{\\nu+1}, &
2K_\\nu' &= -(K_{\\nu-1}+K_{\\nu+1}), & 2J_\\nu' &= J_{\\nu-1}-J_{\\nu+1}.
\\end{aligned}
```

The modified Bessel recurrences are equation 10.29.1 of the NIST Digital Library of
Mathematical Functions. The ordinary Bessel recurrence is equation 10.6.1.
The derivative at the cusp of the absolute value takes `sign(0) = 0`.
"""
module LineCableModelsMeasurementsExt

import LineCableModels.Engine: nominal
using Measurements: Measurement, result, value
using SpecialFunctions: SpecialFunctions, besselix, besselkx, besseljx, besselk

nominal(quantity::Measurement) = value(quantity)

function SpecialFunctions.besselix(ν::Real, z::Complex{<:Measurement})
    x = nominal(z)
    w = besselix(ν, x)
    dw = (besselix(ν - 1, x) + besselix(ν + 1, x)) / 2
    dx = dw - sign(real(x)) * w
    dy = im * dw
    return result(w, (real(dx), real(dy), imag(dx), imag(dy)), z)
end

function SpecialFunctions.besselkx(ν::Real, z::Complex{<:Measurement})
    x = nominal(z)
    w = besselkx(ν, x)
    dx = w - (besselkx(ν - 1, x) + besselkx(ν + 1, x)) / 2
    dy = im * dx
    return result(w, (real(dx), real(dy), imag(dx), imag(dy)), z)
end

function SpecialFunctions.besseljx(ν::Real, z::Complex{<:Measurement})
    x = nominal(z)
    w = besseljx(ν, x)
    dx = (besseljx(ν - 1, x) - besseljx(ν + 1, x)) / 2
    dy = im * dx - sign(imag(x)) * w
    return result(w, (real(dx), real(dy), imag(dx), imag(dy)), z)
end

function SpecialFunctions.besselk(ν::Real, z::Complex{<:Measurement})
    x = nominal(z)
    w = besselk(ν, x)
    dx = -(besselk(ν - 1, x) + besselk(ν + 1, x)) / 2
    dy = im * dx
    return result(w, (real(dx), real(dy), imag(dx), imag(dy)), z)
end

end
