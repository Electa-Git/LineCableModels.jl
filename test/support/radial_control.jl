# Independent radial diffusion control. This solves E''+E'/r-k²E=0 via
# Frobenius recurrences, not the production Bessel functions or surface ratios.
module RadialControl
using LinearAlgebra
export fundamental, surfaces, annular_admittance
function fundamental(z)
    iszero(z) && throw(ArgumentError("use regularity at the axis separately"))
    t=one(z); u=t; du=zero(z); harmonic=zero(real(z))
    weighted=zero(z); dweighted=zero(z); magnitude=abs(t)
    remainder=oftype(abs(z),Inf)
    for n in 1:4096
        t*=z^2/(4n^2)
        harmonic+=inv(oftype(real(z),n))
        u+=t; du+=2n*t/z
        weighted+=harmonic*t; dweighted+=2n*harmonic*t/z
        magnitude+=abs(t)*(1+abs(harmonic))*(1+2n/abs(z))
        ratio=abs(z)^2/(4(n+1)^2)*(n+2)/(n+1)
        if ratio<1/2
            # The ratios of subsequent terms, harmonic weights and derivative
            # weights decrease; this bounds all four omitted series tails.
            next=abs(t)*abs(z)^2/(4(n+1)^2)
            remainder=next*(1+harmonic+inv(oftype(harmonic,n+1)))*(1+2(n+1)/abs(z))/(1-ratio)
            remainder < eps(real(z))*magnitude && break
        end
        n==4096 && error("inconclusive radial reference: series term limit")
    end
    logarithm=log(z/2)+oftype(real(z),Base.MathConstants.eulergamma)
    v=weighted-logarithm*u
    dv=dweighted-u/z-logarithm*du
    # Roundoff allowance accounts for summed absolute terms, log and derivatives.
    bound=(remainder+4096eps(real(z))*magnitude)*(2+abs(logarithm)+inv(abs(z)))
    return (;u,du,v,dv,bound)
end
function surfaces(ri,ro,rho,mu,s)
    pi=oftype(rho,π)
    k=sqrt(s*mu/rho)
    outer=fundamental(k*ro)
    if iszero(ri)
        value=outer.u/outer.du*s*mu/(2pi*ro*k)
        # Denominator perturbation bound for the regular radial solution.
        e=outer.bound
        e<abs(outer.du)/2 || error("inconclusive radial reference denominator")
        bound=abs(s*mu/(2pi*ro*k))*2e*(1+abs(outer.u/outer.du))/abs(outer.du)
        return (inner=zero(value),outer=value,transfer=zero(value),bound)
    end
    inner=fundamental(k*ri)
    derivatives=k.*[inner.du inner.dv; outer.du outer.dv]
    currents=diagm([-s*mu/(2pi*ri),s*mu/(2pi*ro)])
    coefficients=derivatives\currents
    values=[inner.u inner.v;outer.u outer.v]*coefficients
    e=max(inner.bound,outer.bound)
    perturbation=2abs(k)*e*opnorm(inv(derivatives),Inf)
    perturbation<1/2 || error("inconclusive radial reference conditioning")
    bound=4e*opnorm(coefficients,Inf)+2perturbation*opnorm(values,Inf)+4096eps(real(k))*opnorm(values,Inf)
    return (inner=values[1,1],outer=values[2,2],transfer=values[2,1],bound)
end
function annular_admittance(a,b,kappa,mu,s)
    pi=oftype(a,π)
    k=sqrt(s*mu*kappa)
    inner,outer=fundamental(k*a),fundamental(k*b)
    matrix=[inner.u inner.v;outer.u outer.v]
    coefficients=matrix\[one(k),zero(k)]
    derivative=k*(inner.du*coefficients[1]+inner.dv*coefficients[2])
    value=-2pi*a*kappa*derivative
    e=max(inner.bound,outer.bound)
    perturbation=2e*opnorm(inv(matrix),Inf)
    perturbation<1/2 || error("inconclusive annular reference conditioning")
    bound=abs(2pi*a*kappa*k)*(2e*norm(coefficients)+
        2perturbation*norm([inner.du,inner.dv])*norm(coefficients))+4096eps(real(k))*abs(value)
    return (;value,bound)
end
end
