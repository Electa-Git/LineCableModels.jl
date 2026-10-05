# Homogeneous solid-cylinder references for exp(j*omega*t - Gamma*z), SI units.
# The exterior reference is infinity for buried receivers and the point h
# vertically below the conductor for an overhead receiver.
using SpecialFunctions: besselk

function transverse_root(value)
    q=sqrt(complex(value))
    real(q)<0 || (iszero(real(q)) && imag(q)<0) ? -q : q
end

function cylinder_admittance(frequency,radius,sigma,epsilon,mu,gamma;reference_distance=nothing)
    omega=2pi*frequency
    admittivity=complex(sigma,omega*epsilon)
    q=transverse_root(im*omega*mu*admittivity-gamma^2)
    iszero(q) && throw(ArgumentError("exact air cutoff is singular"))
    x=q*radius
    denominator=besselk(0,x)
    reference_distance===nothing || (denominator-=besselk(0,q*reference_distance))
    2pi*admittivity*x*besselk(1,x)/denominator
end
