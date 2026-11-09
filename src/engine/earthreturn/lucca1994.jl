defaults(::Type{Formula{:lucca1994}}) = (;)
domain(::Formula{:lucca1994}) = (mu_r_earth = 1.0 .. 1.0, Γ = 0.0 .. 0.0)

validate(formula::Formula{:lucca1994}, ::typeof(impedance!)) = formula
function validate(formula::Formula{:lucca1994}, t::Int, s::Int)
    (t == 1) != (s == 1) || invoke(validate, Tuple{Formula, Int, Int}, formula, t, s)
    return formula
end

function (f::Formula{:lucca1994})(
        ::Val{:Z}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    h_a = max(input.y[i], input.y[j])
    h_g = -min(input.y[i], input.y[j])
    y = abs(input.x[i] - input.x[j])
    V = h_a + h_g
    σ₂ = inv(media[2].rho)
    p = inv(sqrt(s * μ₀ * σ₂))
    Hc = V + 2p
    S² = Hc^2 + y^2
    D = hypot(V, y)
    return s * μ₀ / (2π) *
           (log(sqrt(S²) / D) - (2 / 3) * (p / S²)^3 * Hc * (Hc^2 - 3y^2))
end
