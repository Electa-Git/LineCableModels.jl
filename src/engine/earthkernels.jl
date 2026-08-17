"Return electrical conductivity, including the open-circuit limit."
@inline conductivity(rho) = isinf(rho) ? zero(rho) : inv(rho)

@inline special_besselix(order::Integer, value) = SpecialFunctions.besselix(order, value)
@inline special_besselkx(order::Integer, value) = SpecialFunctions.besselkx(order, value)
@inline special_besselk(order::Integer, value) = SpecialFunctions.besselk(order, value)

# SpecialFunctions intentionally does not implement complex BigFloat Bessel
# functions. These methods isolate that external numerical boundary without
# reducing the caller's working precision.
function special_besselix(order::Integer, value::Complex{BigFloat})
    order >= 0 || throw(DomainError(order, "Bessel order must be nonnegative"))
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
    order >= 0 || throw(DomainError(order, "Bessel order must be nonnegative"))
    real(value) > 0 || throw(DomainError(
        value, "complex BigFloat besselk requires a positive real part"
    ))
    tolerance = max(sqrt(eps(BigFloat)), BigFloat("1e-30"))
    integrand(t) = exp(-value * cosh(t)) * cosh(BigFloat(order) * t)
    result, _ = quadgk(
        integrand, zero(BigFloat), BigFloat(Inf); rtol = tolerance
    )
    return result
end

function special_besselkx(order::Integer, value::Complex{BigFloat})
    exp(value) * special_besselk(order, value)
end

"Evaluate the shared earth-return modified-Bessel difference."
@inline function bessel_difference(gamma_s, inner, outer)
    maximum_argument = max(abs(gamma_s) * inner, abs(gamma_s) * outer)
    isapprox(nominal(maximum_argument), 0; atol = 1.0e-6) &&
        return log(outer / inner)
    return special_besselk(0, gamma_s * inner) -
           special_besselk(0, gamma_s * outer)
end
