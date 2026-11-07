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

Each method evaluates the function and its two partial derivatives at the value
``z_0 = x_0 + i y_0`` of the argument, then returns the first-order expansion

```math
f(z) = f(z_0) + \\partial_x f(z_0)\\,(x - x_0) + \\partial_y f(z_0)\\,(y - y_0)
```

in dual arithmetic. The differences ``x - x_0`` and ``y - y_0`` have the value zero and
the partial derivatives of the argument. A nested dual passes its inner dual to the same
method, so derivatives of higher order follow.

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
    x = complex(value(real(z)), value(imag(z)))
    w = besselix(ν, x)
    dw = (besselix(ν - 1, x) + besselix(ν + 1, x)) / 2
    dx = dw - sign(real(x)) * w
    dy = im * dw
    return w + dx * (real(z) - real(x)) + dy * (imag(z) - imag(x))
end

function SpecialFunctions.besselkx(ν::Real, z::Complex{<:Dual})
    x = complex(value(real(z)), value(imag(z)))
    w = besselkx(ν, x)
    dx = w - (besselkx(ν - 1, x) + besselkx(ν + 1, x)) / 2
    dy = im * dx
    return w + dx * (real(z) - real(x)) + dy * (imag(z) - imag(x))
end

function SpecialFunctions.besseljx(ν::Real, z::Complex{<:Dual})
    x = complex(value(real(z)), value(imag(z)))
    w = besseljx(ν, x)
    dx = (besseljx(ν - 1, x) - besseljx(ν + 1, x)) / 2
    dy = im * dx - sign(imag(x)) * w
    return w + dx * (real(z) - real(x)) + dy * (imag(z) - imag(x))
end

function SpecialFunctions.besselk(ν::Real, z::Complex{<:Dual})
    x = complex(value(real(z)), value(imag(z)))
    w = besselk(ν, x)
    dx = -(besselk(ν - 1, x) + besselk(ν + 1, x)) / 2
    dy = im * dx
    return w + dx * (real(z) - real(x)) + dy * (imag(z) - imag(x))
end

end
