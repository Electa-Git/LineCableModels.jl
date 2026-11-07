"""
    LineCableModelsIntervalArithmeticExt

Enclose the Bessel functions of SpecialFunctions on an interval argument with the ball
arithmetic of Arb. Define the `nominal` method, the midpoint of an interval, and five
methods of SpecialFunctions Bessel functions. Four take a complex interval argument, and
`besseljx` also takes a real one.

Each method forms the Arb ball of its argument at 128 bits, which holds the box of the
intervals exactly. It evaluates the function in Arb and returns the result as intervals
rounded outward, in the number type of the argument's intervals. The conversions to and from
Arb are those of IntervalArithmetic, which marks each enclosure as guaranteed.

For ``z = x + iy``, the functions are those of SpecialFunctions:
``\\mathrm{e}^{-|x|} I_\\nu(z)`` (`besselix`), ``\\mathrm{e}^{z} K_\\nu(z)`` (`besselkx`),
``\\mathrm{e}^{-|y|} J_\\nu(z)`` (`besseljx`) and ``K_\\nu(z)`` (`besselk`). Arb gives the
second and the fourth. The methods multiply the unscaled ``I_\\nu`` and ``J_\\nu`` of Arb by
their scaling, in ball arithmetic.

The engine passes a complex argument to these functions, except in the kernels of
`:martins2026`, whose `besseljx` takes the real argument ``r \\lambda``. Kernel branches and
integration breakpoints decide on the midpoint of an interval, `nominal`. A computation with
interval inputs takes the branches and the breakpoints of its midpoints.
"""
module LineCableModelsIntervalArithmeticExt

import LineCableModels.Commons: nominal
using Arblib: Acb, Arb
using IntervalArithmetic: Interval, interval, mid
using SpecialFunctions: SpecialFunctions, besseli, besselj, besselk, besselkx

nominal(quantity::Interval) = mid(quantity)

function SpecialFunctions.besselix(ν::Real, z::Complex{Interval{T}}) where {T}
    w = Acb(z; prec = 128)
    return interval(T, besseli(Acb(ν), w) * exp(-abs(real(w))))
end

function SpecialFunctions.besselkx(ν::Real, z::Complex{Interval{T}}) where {T}
    w = Acb(z; prec = 128)
    return interval(T, besselkx(Acb(ν), w))
end

function SpecialFunctions.besseljx(ν::Real, z::Complex{Interval{T}}) where {T}
    w = Acb(z; prec = 128)
    return interval(T, besselj(Acb(ν), w) * exp(-abs(imag(w))))
end

function SpecialFunctions.besselk(ν::Real, z::Complex{Interval{T}}) where {T}
    w = Acb(z; prec = 128)
    return interval(T, besselk(Acb(ν), w))
end

function SpecialFunctions.besseljx(ν::Real, x::Interval{T}) where {T}
    w = Arb(x; prec = 128)
    return interval(T, besselj(Arb(ν), w))
end

end
