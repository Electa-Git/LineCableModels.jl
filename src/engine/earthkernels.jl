"Return electrical conductivity, including the open-circuit limit."
@inline conductivity(rho) = isinf(rho) ? zero(rho) : inv(rho)

# Arithmetic with a zero longitudinal term can change a signed zero on the
# lossless propagating cut. Select its outgoing boundary value from frequency.
@inline function spectral_root(value::Complex{T},jω::Complex{T}) where {T <: Real}
    if iszero(imag(value)) && real(value)<zero(T)
        sense=signbit(imag(jω)) ? -one(T) : one(T)
        return complex(zero(T),sense*sqrt(-real(value)))
    end
    return sqrt(value)
end

@inline special_besselix(order::Integer, value) = SpecialFunctions.besselix(order, value)
@inline special_besselkx(order::Integer, value) = SpecialFunctions.besselkx(order, value)
@inline special_besselk(order::Integer, value) = SpecialFunctions.besselk(order, value)

# SpecialFunctions omits complex BigFloat Bessel functions. The local methods
# retain the caller's working precision for unsupported argument types.
function special_besselix(order::Integer, value::Complex{BigFloat})
    order >= 0 || throw(DomainError(order, "Bessel order must be nonnegative"))
    bits=precision(BigFloat)
    # Complex power-series terms grow with |z|, while the retained function
    # grows with |Re(z)|. Preserve the digits lost in their cancellation.
    extra=max(0,ceil(Int,(abs(value)-abs(real(value)))/log(BigFloat(2))))+24
    result=setprecision(BigFloat,bits+extra) do
        _big_besselix_series(order,Complex{BigFloat}(value))
    end
    return complex(BigFloat(real(result);precision=bits),
        BigFloat(imag(result);precision=bits))
end

function _big_besselix_series(order,value)
    half = value / BigFloat(2)
    term = half^order / BigFloat(factorial(big(order)))
    result = term
    for index in 1:100_000
        term *= half^2 / (BigFloat(index) * BigFloat(index + order))
        next = result + term
        if next == result || abs(term) <= eps(BigFloat) * abs(next)
            return exp(-abs(real(value))) * next
        end
        result = next
    end
    throw(ErrorException("complex BigFloat besseli series did not converge"))
end

# DLMF 10.31.2 and K0'=-K1 (10.29.3), without a small-argument truncation.
# The convergent series avoids resolving a very long hyperbolic-integral tail.
function _small_besselk(order,z::Complex{BigFloat})
    half=z/2; square=half^2
    term=one(z); harmonic=zero(BigFloat)
    sum0=zero(z); sum1=zero(z)
    for k in 1:10000
        term*=square/(BigFloat(k)^2)
        harmonic+=inv(BigFloat(k))
        addition=harmonic*term
        sum0+=addition; sum1+=k*addition
        if abs(addition)<=eps(BigFloat)*max(abs(sum0),abs(sum1))
            i0=special_besselix(0,z)*exp(abs(real(z)))
            logarithm=log(half)+BigFloat(Base.MathConstants.eulergamma)
            if order==0
                return sum0-logarithm*i0
            end
            i1=special_besselix(1,z)*exp(abs(real(z)))
            return i0/z+logarithm*i1-2sum1/z
        end
    end
    throw(ErrorException("complex BigFloat small-argument K series did not converge"))
end

function special_besselk(order::Integer, value::Complex{BigFloat})
    order >= 0 || throw(DomainError(order, "Bessel order must be nonnegative"))
    real(value)<0 && return exp(-value)*_continued_besselkx(order,value)
    real(value) >= 0 && !iszero(value) || throw(DomainError(
        value, "complex BigFloat besselk requires a nonzero right-half-plane argument"
    ))
    order<=1 && abs(value)<=1 && return _small_besselk(order,value)
    real(value)<abs(imag(value))/10 && return _rotated_besselk(order,value)
    tolerance = max(sqrt(eps(BigFloat)), BigFloat("1e-30"))
    integrand(t) = exp(-value * cosh(t)) * cosh(BigFloat(order) * t)
    result, _ = quadgk(
        integrand, zero(BigFloat), BigFloat(Inf); rtol = tolerance
    )
    return result
end

# DLMF 10.32.8, t=1+u/z followed by contour rotation and u=x².
# The deformation avoids the negative-real cut of (1+u/(2z)) and
# extends the decaying integral to the imaginary-axis boundary.
# https://dlmf.nist.gov/10.32.E8
function _rotated_besselk(order::Integer,value::Complex{T}) where {T <: AbstractFloat}
    order>=0 && real(value)>=0 && !iszero(value) ||
        throw(DomainError((order,value),"rotated Bessel integral requires nonnegative order and a nonzero right-half-plane argument"))
    half=one(T)/2
    integrand(x)=2exp(-x^2)*x^(2order)*(1+x^2/(2value))^(order-half)
    transition=min(sqrt(2abs(value)),one(T))
    integral=quadgk(integrand,zero(T),transition,T(Inf);
        rtol=max(sqrt(eps(T)),T(1e-30)))[1]
    return sqrt(T(π))*exp(-value)/sqrt(2value)/
        SpecialFunctions.gamma(T(order)+half)*integral
end

function special_besselkx(order::Integer, value::Complex{BigFloat})
    real(value)<0 && return _continued_besselkx(order,value)
    return exp(value) * special_besselk(order, value)
end

# DLMF 10.34.5, m=+/-1, in scaled form. Required when a fitted
# complex-image distance puts its Bessel argument in the left half-plane.
# https://dlmf.nist.gov/10.34.E5
function _continued_besselkx(order::Integer,value::Complex{BigFloat})
    order>=0 && !iszero(imag(value)) || throw(DomainError((order,value),
        "Bessel continuation needs a nonnegative order and an off-cut argument"))
    parity=isodd(order) ? -1 : 1
    sense=sign(imag(value))
    return parity*exp(2value)*special_besselkx(order,-value)-
        complex(zero(BigFloat),sense*BigFloat(π))*
        exp(complex(zero(BigFloat),imag(value)))*special_besselix(order,-value)
end

# Q(w)=exp(-w) E1(-w), with the principal E1 branch off its cut.
# Small arguments use its convergent series, retaining complex BigFloat.
function scaled_expint_negative(w::Complex{T}) where {T <: AbstractFloat}
    iszero(w) && throw(DomainError(w,"scaled exponential integral is singular at zero"))
    if abs(w)<=one(T)
        term=w; series=w
        for k in 2:10000
            term*=w/k
            addition=term/k
            series+=addition
            if abs(addition)<=eps(T)*max(one(T),abs(series))
                return exp(-w)*(-T(Base.MathConstants.eulergamma)-log(-w)-series)
            end
        end
        throw(ErrorException("scaled exponential-integral series did not converge"))
    end
    if T <: Union{Float32,Float64} && abs(real(w))<500
        wide=ComplexF64(w)
        return Complex{T}(exp(-wide)*SpecialFunctions.expint(-wide))
    end
    iszero(imag(w)) && real(w)>zero(T) &&
        throw(DomainError(w,"the exponential-integral cut requires a specified boundary value"))
    integrand(t)=exp(-t)/(t-w)
    points=real(w)>zero(T) ? (zero(T),real(w),T(Inf)) : (zero(T),T(Inf))
    result,_=quadgk(integrand,points...;rtol=max(sqrt(eps(T)),T(1e-35)))
    return Complex{T}(result)
end

"Evaluate the shared earth-return modified-Bessel difference."
@inline function bessel_difference(gamma_s, inner, outer)
    maximum_argument = max(abs(gamma_s) * inner, abs(gamma_s) * outer)
    isapprox(nominal(maximum_argument), 0; atol = 1.0e-6) &&
        return log(outer / inner)
    return special_besselk(0, gamma_s * inner) -
           special_besselk(0, gamma_s * outer)
end
