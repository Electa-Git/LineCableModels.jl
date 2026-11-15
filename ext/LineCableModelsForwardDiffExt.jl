"""
    LineCableModelsForwardDiffExt

Compute derivatives with the dual numbers of ForwardDiff. A `Dual` input gives series
impedances `Z` \\[Ω/m\\] and shunt admittances `Y` \\[S/m\\] with dual real and imaginary
parts. Their partials are the derivatives with respect to the inputs. Define the `nominal`
method and four methods of SpecialFunctions Bessel functions on a complex dual argument.

For ``z = x + iy``, write the scaled functions as
``I^s_\\nu(z) = \\mathrm{e}^{-|x|} I_\\nu(z)`` (`besselix`),
``K^s_\\nu(z) = \\mathrm{e}^{z} K_\\nu(z)`` (`besselkx`), and
``J^s_\\nu(z) = \\mathrm{e}^{-|y|} J_\\nu(z)`` (`besseljx`), with ``i^2 = -1``.
The unscaled ``K_\\nu(z)`` is `besselk`.
At the argument's value ``z_0 = x_0 + i y_0``, the order recurrences give the expansions

```math
\\begin{aligned}
I^s_\\nu(z) &\\approx I^s_\\nu(z_0)
    + \\tfrac12\\left(I^s_{\\nu-1}+I^s_{\\nu+1}\\right)(z_0)\\,(z-z_0)
    - I^s_\\nu(z_0)\\,(|x|-|x_0|), \\\\
K^s_\\nu(z) &\\approx K^s_\\nu(z_0)
    + \\left(K^s_\\nu-\\tfrac12\\left(K^s_{\\nu-1}+K^s_{\\nu+1}\\right)\\right)(z_0)\\,(z-z_0), \\\\
J^s_\\nu(z) &\\approx J^s_\\nu(z_0)
    + \\tfrac12\\left(J^s_{\\nu-1}-J^s_{\\nu+1}\\right)(z_0)\\,(z-z_0)
    - J^s_\\nu(z_0)\\,(|y|-|y_0|), \\\\
K_\\nu(z) &\\approx K_\\nu(z_0)
    - \\tfrac12\\left(K_{\\nu-1}+K_{\\nu+1}\\right)(z_0)\\,(z-z_0).
\\end{aligned}
```

The modified Bessel recurrences are equation 10.29.1 of the NIST Digital Library of
Mathematical Functions. The ordinary Bessel recurrence is equation 10.6.1.
The absolute value has a cusp at ``x = 0`` for ``I^s_\\nu`` and at ``y = 0`` for ``J^s_\\nu``.
`besselix` and `besseljx` evaluate the absolute value through the number type's own `abs`.
A scaled function times the factor that removes its scaling keeps the derivative of the
unscaled function.

In dual expansions, ``z - z_0`` and differences of absolute values have the value zero
and the partial derivatives of the argument and the number type's own `abs`, respectively.
A nested dual passes its inner dual to the same method, so derivatives of
higher order follow.

The internal impedance `:schelkunoff1934`, the direct and image means of the earth-return
formulas and the cable states of `:martins2026` pass a complex argument. So do the closed
forms of `:pollaczek1926`, `:papadopoulos2010`, `:saad1996` and `:xue2018`. The kernels of
`:martins2026` call `besseljx` with the real argument ``r \\lambda``, which takes the rules
that ForwardDiff defines for these functions on a real dual.
"""
module LineCableModelsForwardDiffExt

import LineCableModels.Commons: nominal
using ForwardDiff: Dual, value
using SpecialFunctions: SpecialFunctions, besselix, besselkx, besseljx, besselk

nominal(quantity::Dual) = nominal(value(quantity))

function SpecialFunctions.besselix(ν::Real, z::Complex{<:Dual})
    z₀ = complex(value(real(z)), value(imag(z)))
    w = besselix(ν, z₀)
    dw = (besselix(ν - 1, z₀) + besselix(ν + 1, z₀)) / 2
    return w + dw * (z - z₀) - w * (abs(real(z)) - abs(real(z₀)))
end

function SpecialFunctions.besselkx(ν::Real, z::Complex{<:Dual})
    z₀ = complex(value(real(z)), value(imag(z)))
    w = besselkx(ν, z₀)
    dw = w - (besselkx(ν - 1, z₀) + besselkx(ν + 1, z₀)) / 2
    return w + dw * (z - z₀)
end

function SpecialFunctions.besseljx(ν::Real, z::Complex{<:Dual})
    z₀ = complex(value(real(z)), value(imag(z)))
    w = besseljx(ν, z₀)
    dw = (besseljx(ν - 1, z₀) - besseljx(ν + 1, z₀)) / 2
    return w + dw * (z - z₀) - w * (abs(imag(z)) - abs(imag(z₀)))
end

function SpecialFunctions.besselk(ν::Real, z::Complex{<:Dual})
    z₀ = complex(value(real(z)), value(imag(z)))
    w = besselk(ν, z₀)
    dw = -(besselk(ν - 1, z₀) + besselk(ν + 1, z₀)) / 2
    return w + dw * (z - z₀)
end

end
