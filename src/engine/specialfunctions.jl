@inline special_besselix(order::Integer, value) = SpecialFunctions.besselix(order, value)
@inline special_besselkx(order::Integer, value) = SpecialFunctions.besselkx(order, value)
@inline special_besselk(order::Integer, value) = SpecialFunctions.besselk(order, value)
@inline special_besseljx(order::Integer, value) = SpecialFunctions.besseljx(order, value)
@inline special_besselix(order::Integer,
    value::Complex{Float32}) = ComplexF32(SpecialFunctions.besselix(Float32(order), value))
@inline special_besselkx(order::Integer,
    value::Complex{Float32}) = ComplexF32(SpecialFunctions.besselkx(Float32(order), value))
@inline special_besselk(order::Integer,
    value::Complex{Float32}) = ComplexF32(SpecialFunctions.besselk(Float32(order), value))
@inline special_besseljx(order::Integer,
    value::Complex{Float32}) = ComplexF32(SpecialFunctions.besseljx(Float32(order), value))
function special_besseljx(order::Integer, value::Complex{BigFloat})
    order>=0 || throw(DomainError(order, "Bessel order must be nonnegative"))
    abs(value)>max(128, precision(BigFloat), order^2) ||
        return (-im)^order*special_besselix(order, im*value)
    real(value)<0 && return (-1)^order*special_besseljx(order, -value)
    # DLMF 10.17.1, 10.17.5–6: combine the two scaled Hankel expansions.
    # This avoids an angular grid proportional to a large complex argument.
    plus=one(value)
    minus=one(value)
    tp=one(value)
    tm=one(value)
    for n in 1:100_000
        factor=(4BigFloat(order)^2-BigFloat(2n-1)^2)/(8n*value)
        tp*=im*factor
        tm*=-im*factor
        plus+=tp
        minus+=tm
        max(abs(tp), abs(tm))<=eps(BigFloat)*max(abs(plus), abs(minus)) && break
        n==100_000 &&
            throw(ErrorException("scaled BigFloat Bessel asymptotic did not converge"))
    end
    phase=cis(BigFloat(π)*(BigFloat(order)/2+BigFloat(1)/4))
    growth=abs(imag(value))
    return (exp(im*value-growth)*conj(phase)*plus +
            exp(-im*value-growth)*phase*minus)/sqrt(2BigFloat(π)*value)
end

# SpecialFunctions omits complex BigFloat Bessel functions. The local methods
# retain the caller's working precision for unsupported argument types.
function special_besselix(order::Integer, value::Complex{BigFloat})
    order >= 0 || throw(DomainError(order, "Bessel order must be nonnegative"))
    abs(value)>max(128, precision(BigFloat), order^2) &&
        return (-im)^order*special_besseljx(order, im*value)
    if abs(value)>8
        # DLMF 10.32.3, with the scaling inside the exponential. Seed the
        # angular phase scale rather than relying on two aliased rules.
        count=max(8, ceil(Int, abs(imag(value))+order))
        points=collect(range(zero(BigFloat), BigFloat(π); length = count+1))
        f=θ->exp(value*cos(θ)-abs(real(value)))*cos(order*θ)
        result, error=quadgk(f, points; rtol = sqrt(eps(BigFloat)))
        return result/BigFloat(π)
    end
    half = value / BigFloat(2)
    term = half^order / BigFloat(factorial(big(order)))
    result = term
    for index in 1:100_000
        term *= half^2 / (BigFloat(index) * BigFloat(index + order))
        next = result + term
        if next == result || abs(term) <= eps(BigFloat) * max(abs(next), one(BigFloat))
            return exp(-abs(real(value))) * next
        end
        result = next
    end
    throw(ErrorException("complex BigFloat besseli series did not converge"))
end

function special_besselk(order::Integer, value::Complex{BigFloat})
    return exp(-value)*special_besselkx(order, value)
end

function special_besselkx(order::Integer, value::Complex{BigFloat})
    order>=0 || throw(DomainError(order, "Bessel order must be nonnegative"))
    !iszero(value)&&real(value)>=0 || throw(DomainError(value,
        "complex BigFloat K requires a nonzero argument on the outgoing right half-plane"))
    # DLMF 10.32.8, w=z(t−1), rotate the w contour to the positive real
    # axis, then w=u². This remains exponentially decaying in the lossless
    # imaginary-argument limit; exp(-z*cosh(t)) does not.
    exponent=BigFloat(order)-BigFloat(1)/2
    f=u->exp(-u*u)*u^(2order)*(1+u*u/(2value))^exponent
    feature=min(sqrt(abs(value)), BigFloat(1)/2)
    result,
    error=quadgk(f, zero(BigFloat), feature, one(BigFloat), BigFloat(Inf);
        rtol = sqrt(eps(BigFloat)))
    return 2sqrt(BigFloat(π)/(2value))/SpecialFunctions.gamma(BigFloat(order)+BigFloat(1)/2)*result
end
