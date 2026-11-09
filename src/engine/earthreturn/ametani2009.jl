defaults(::Type{Formula{:ametani2009}}) = (;)
domain(::Formula{:ametani2009}) = (mu_r_earth = 1.0 .. 1.0, Γ = 0.0 .. 0.0)

validate(formula::Formula{:ametani2009}, ::typeof(impedance!)) = formula
function validate(formula::Formula{:ametani2009}, t::Int, s::Int)
    (t == 1) != (s == 1) || invoke(validate, Tuple{Formula, Int, Int}, formula, t, s)
    return formula
end

function (f::Formula{:ametani2009})(
        ::Val{:Z}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    h_a = max(input.y[i], input.y[j])
    h_g = -min(input.y[i], input.y[j])
    y = abs(input.x[i] - input.x[j])
    m = sqrt(s * μ₀ / media[2].rho)
    p = inv(m)
    S = sqrt((h_a + h_g + 2p)^2 + y^2)
    D = hypot(h_a + h_g, y)
    return s * μ₀ / (2π) * exp(-h_g / p) * log(S / D)
end
