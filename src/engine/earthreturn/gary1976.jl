defaults(::Type{Formula{:gary1976}}) = (;)
domain(::Formula{:gary1976}) = (mu_r_earth = 1.0 .. 1.0, Γ = 0.0 .. 0.0)

validate(formula::Formula{:gary1976}, ::typeof(impedance!)) = formula
function validate(formula::Formula{:gary1976}, t::Int, s::Int)
    t == 1 && s == 1 || invoke(validate, Tuple{Formula, Int, Int}, formula, t, s)
    return formula
end

function (f::Formula{:gary1976})(::Val{:Z}, ::Val{:self}, i, buffers, input, media, s, k)
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    σ₂ = inv(media[2].rho)
    p = inv(sqrt(s * μ₀ * σ₂))
    return s * μ₀ / (2π) *
           (Λ(f, Val(:self), zero(s), r_i, h_i) + log((h_i + p) / h_i))
end

function (f::Formula{:gary1976})(
        ::Val{:Z}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    h_j = abs(input.y[j])
    y = abs(input.x[i] - input.x[j])
    H = h_i + h_j
    D = hypot(y, h_i - h_j)
    D′ = hypot(y, H)
    σ₂ = inv(media[2].rho)
    p = inv(sqrt(s * μ₀ * σ₂))
    return s * μ₀ / (2π) *
           (Λ(f, Val(:mutual), zero(s), r_i, D, D′) +
            log(sqrt(y^2 + (H + 2p)^2) / D′))
end
