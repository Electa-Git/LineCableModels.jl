# Homogeneous solid-cylinder references for exp(j*omega*t - Gamma*z), SI units.
# The exterior reference is infinity for buried receivers and the point h
# vertically below the conductor for an overhead receiver.
using SpecialFunctions: besselk, besselix

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

# The voltage-reference correction follows Z = K + Gamma^2 P, with native P = 1/Y.
# Scaled interior Bessel functions avoid overflow for highly conducting cylinders.
function cylinder_impedance(frequency,radius,sigma,epsilon,mu,gamma,
        conductor_sigma,conductor_epsilon,conductor_mu;reference_distance=nothing)
    omega=2pi*frequency
    q=transverse_root(im*omega*mu*complex(sigma,omega*epsilon)-gamma^2)
    iszero(q) && throw(ArgumentError("exact transverse cutoff is singular"))
    x=q*radius
    exterior=im*omega*mu*besselk(0,x)/(2pi*x*besselk(1,x))
    if reference_distance!==nothing
        exterior-=gamma^2*besselk(0,q*reference_distance)/
            (2pi*complex(sigma,omega*epsilon)*x*besselk(1,x))
    end
    admittivity=complex(conductor_sigma,omega*conductor_epsilon)
    kc=transverse_root(im*omega*conductor_mu*admittivity-gamma^2)
    interior=kc/(2pi*radius*admittivity)*besselix(0,kc*radius)/besselix(1,kc*radius)
    exterior+interior
end
