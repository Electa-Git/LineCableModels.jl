# Alternative evaluations of one Pollaczek integral, not new earth models.
# Higher precision protects the cancellation in the source decomposition.
# Arithmetic uses the caller's current BigFloat precision; no global precision
# setting is changed during a threaded line-parameter calculation.
function _half_besseli(n,z)
    ν=BigFloat(n)+BigFloat(1)/2
    term=one(z); total=term
    for k in 1:10000
        term*=z^2/(4k*(ν+k)); total+=term
        abs(term)<=eps(BigFloat)*max(one(BigFloat),abs(total)) &&
            return (z/2)^ν/gamma(ν+1)*total
    end
    throw(ErrorException("half-order Bessel I series did not converge"))
end

function _half_besselk(n,z)
    term=one(z); total=term
    for k in 1:n
        term*=(n+k)*(n-k+1)/(2k*z); total+=term
    end
    return sqrt(BigFloat(π)/(2z))*exp(-z)*total
end

function _hyp1f1_one(b,z)
    term=one(z); total=term
    for n in 1:10000
        term*=z/(b+n-1); total+=term
        abs(term)<=eps(BigFloat)*max(one(BigFloat),abs(total)) && return total
    end
    throw(ErrorException("confluent hypergeometric series did not converge"))
end

function _pollaczek_auxiliary(::Val{:finite_integral},k,H,x)
    R=hypot(H,x)
    base=(H^2-x^2)/R^4*exp(-k*H)*(1+k*H)
    iszero(x) && return base
    # t=cos(theta) removes the square-root endpoint in source equation (22).
    θmax=atan(abs(x),H)
    integral=quadgk(θ->cos(2θ)*exp(-k*R*cos(θ)),
        zero(H),θmax;rtol=max(sqrt(eps(BigFloat)),BigFloat("1e-35")))[1]
    return base-k^2*abs(x)*H/R^2*integral
end

function _pollaczek_auxiliary(::Val{:bessel_product},k,H,x)
    iszero(x) && return exp(-k*H)*(1+k*H)/H^2
    R=hypot(H,x)
    abs(x)>4H && return _pollaczek_auxiliary(Val(:finite_integral),k,H,x)
    z=k*x^2/(2*(R+H)); Z=k*(R+H)/2
    total=zero(k); tolerance=max(sqrt(eps(BigFloat)),BigFloat("1e-35"))
    for n in 0:1000
        term=(-1)^n*(2n+1)^2*(z*_half_besseli(n+1,z)*_half_besselk(n,Z)+
            Z*_half_besseli(n,z)*_half_besselk(n+1,Z))
        total+=term
        n>2 && abs(term)<=tolerance*abs(total) && return 2total/(abs(x)*R)
    end
    return _pollaczek_auxiliary(Val(:finite_integral),k,H,x)
end

function _pollaczek_auxiliary(::Val{:single_bessel},k,H,x)
    R=hypot(H,x); z=k*R
    abs(x)>H && return _pollaczek_auxiliary(Val(:finite_integral),k,H,x)
    # The n=0 term is required by the source parent (9a), including x=0.
    factor=-one(k)/z^(BigFloat(3)/2); total=zero(k)
    tolerance=max(sqrt(eps(BigFloat)),BigFloat("1e-35"))
    for n in 0:1000
        term=factor*_half_besselk(n+1,z); total+=term
        n>2 && abs(term)<=tolerance*abs(total) &&
            return -k^2*sqrt(2/BigFloat(π))*k*H*total
        factor*=(k*x)^2/(2(n+1)*z)*(2n-1)/(2n+1)
    end
    return _pollaczek_auxiliary(Val(:finite_integral),k,H,x)
end

function _pollaczek_auxiliary(::Val{:hypergeometric},k,H,x)
    R=hypot(H,x); δ=x^2/(R*(R+H))
    base=(H^2-x^2)/R^4*exp(-k*H)*(1+k*H)
    iszero(x) && return base
    coefficient=-one(BigFloat); total=zero(k)
    tolerance=max(sqrt(eps(BigFloat)),BigFloat("1e-35"))
    for n in 0:2000
        term=coefficient*δ^(n+BigFloat(1)/2)*
            (_hyp1f1_one(BigFloat(n)+BigFloat(3)/2,-k*(R-H))*
            (1-k*R/2*(2n-1)/(2n+1))-1)
        total+=term
        if n>2 && abs(term)<=tolerance*max(abs(total),one(BigFloat))
            # The common exponential follows from the finite integral (22).
            return base+exp(-k*H)*2sqrt(BigFloat(2))*k*abs(x)*H/R^3*total
        end
        coefficient*=(2n-1)/(4BigFloat(n+1))
    end
    return _pollaczek_auxiliary(Val(:finite_integral),k,H,x)
end


function _pollaczek_auxiliary(::Val{:recursive},k,H,x)
    R=hypot(H,x); a=H/R; b=abs(x)/R; θ=atan(abs(x),H)
    base=(H^2-x^2)/R^4*exp(-k*H)*(1+k*H)
    (iszero(x) || abs(k*R)>30) && return base
    A0=(θ-a*b)/2; A1=b^3/3; B0=θ; B1=b
    power=-k*R
    total=2A0-B0+power*(2A1-B1)
    apower=a; small_terms=0
    tolerance=max(sqrt(eps(BigFloat)),BigFloat("1e-35"))
    for n in 2:2000
        # Integration-by-parts recurrences of the source-defined moments (5).
        An=(apower*b^3+(n-1)*A0)/(n+2)
        Bn=(apower*b+(n-1)*B0)/n
        power*=-k*R/n
        term=power*(2An-Bn); total+=term
        small_terms=abs(term)<=tolerance*max(one(BigFloat),abs(total)) ? small_terms+1 : 0
        small_terms>=3 && return base+k^2*abs(x)*H/R^2*total
        A0,A1=A1,An; B0,B1=B1,Bn; apower*=a
    end
    return _pollaczek_auxiliary(Val(:finite_integral),k,H,x)
end

function _pollaczek_series_coefficient(state,pair,method)
    geometry=_geometry(pair)
    k=state.gamma[2]
    # Keep the spectral evaluation for argument ranges where the decomposition
    # would subtract exponentially different terms, and for non-floating types.
    T=typeof(real(state.jω))
    if !(T<:AbstractFloat) || precision(BigFloat)<128
        return nothing
    end
    if method===Val(:recursive) && abs(k*geometry.D_ij)>30
        # Source-prescribed zero residual above |D/p|=30.
        bracket=_sunde_image_terms(k,geometry.d_ij,geometry.H,geometry.y_ij)
        return state.jω*state.mu[1]/(2*(one(geometry.H)*π))*bracket
    end
    threshold=method===Val(:recursive) ? 30 : 20
    abs(k*geometry.D_ij)>threshold && return nothing
    kb=Complex{BigFloat}(k); H=BigFloat(geometry.H); x=BigFloat(geometry.y_ij)
    R=hypot(H,x); d=BigFloat(geometry.d_ij)
    auxiliary=_pollaczek_auxiliary(method,kb,H,x)
    J=(H/R)^2*special_besselk(0,kb*R)+
        (2*(H/R)^2-1)/(kb*R)*special_besselk(1,kb*R)-auxiliary/kb^2
    bracket=special_besselk(0,kb*d)-special_besselk(0,kb*R)+2J
    value=Complex{BigFloat}(state.jω)*BigFloat(state.mu[1])/(2BigFloat(π))*bracket
    return Complex{T}(value)
end
