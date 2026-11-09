defaults(::Type{Formula{:wedepohl1973}}) = (;)
domain(::Formula{:wedepohl1973}) = (Γ = 0.0 .. 0.0, mD = 0.0 .. 0.25)

validate(formula::Formula{:wedepohl1973}, ::typeof(impedance!)) = formula
function validate(formula::Formula{:wedepohl1973}, t::Int, s::Int)
    t > 1 && s > 1 || invoke(validate, Tuple{Formula, Int, Int}, formula, t, s)
    return formula
end

function (f::Formula{:wedepohl1973})(
        ::Val{:Z}, ::Val{:self}, i, buffers, input, media, s, k)
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    μ₂ = μ₀ * media[2].mu_r
    m = sqrt(s * μ₂ / media[2].rho)
    return s * μ₂ / (2π) *
           (-log(exp(MathConstants.γ) * m * r_i / 2) + 1 / 2 - 4m * h_i / 3)
end

function (f::Formula{:wedepohl1973})(
        ::Val{:Z}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    h_i = abs(input.y[i])
    h_j = abs(input.y[j])
    y = abs(input.x[i] - input.x[j])
    H = h_i + h_j
    D = hypot(y, h_i - h_j)
    μ₂ = μ₀ * media[2].mu_r
    m = sqrt(s * μ₂ / media[2].rho)
    return s * μ₂ / (2π) *
           (-log(exp(MathConstants.γ) * m * D / 2) + 1 / 2 - 2m * H / 3)
end
