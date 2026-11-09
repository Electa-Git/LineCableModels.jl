defaults(::Type{Formula{:saad1996}}) = (;)
domain(::Formula{:saad1996}) = (mu_r_earth = 1.0 .. 1.0, Γ = 0.0 .. 0.0)

validate(formula::Formula{:saad1996}, ::typeof(impedance!)) = formula
function validate(formula::Formula{:saad1996}, t::Int, s::Int)
    t > 1 && s > 1 || invoke(validate, Tuple{Formula, Int, Int}, formula, t, s)
    return formula
end

function (f::Formula{:saad1996})(::Val{:Z}, ::Val{:self}, i, buffers, input, media, s, k)
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    m = sqrt(s * μ₀ / media[2].rho)
    return s * μ₀ / (2π) *
           (besselk(0, m * r_i) + 2exp(-2m * h_i) / (4 + m^2 * r_i^2))
end

function (f::Formula{:saad1996})(
        ::Val{:Z}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    h_i = abs(input.y[i])
    h_j = abs(input.y[j])
    y = abs(input.x[i] - input.x[j])
    H = h_i + h_j
    D = hypot(y, h_i - h_j)
    m = sqrt(s * μ₀ / media[2].rho)
    return s * μ₀ / (2π) *
           (besselk(0, m * D) + 2exp(-m * H) / (4 + m^2 * y^2))
end
