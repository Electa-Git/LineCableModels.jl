defaults(::Type{Formula{:carson1926}}) = (integral = integrate => :quad,)
domain(::Formula{:carson1926}) = (mu_r_earth = 1.0 .. 1.0, Γ = 0.0 .. 0.0)

validate(formula::Formula{:carson1926}, ::typeof(impedance!)) = formula
function validate(formula::Formula{:carson1926}, t::Int, s::Int)
    t == 1 && s == 1 || invoke(validate, Tuple{Formula, Int, Int}, formula, t, s)
    return formula
end

function initialize_buffers(f::Formula{:carson1926}, input::CoaxialInput{T},
        buffers) where {T}
    return initialize_buffers(SpectralIntegral, f.options.integral.method, T, buffers)
end

function (::Formula{:carson1926})(::Val{:J}, λ, γ², H, y)
    return (exp(-(H + im * y) * λ) + exp(-(H - im * y) * λ)) / (2 * (λ + sqrt(λ^2 + γ²)))
end

function (f::Formula{:carson1926})(::Val{:Z}, ::Val{:self}, i, buffers, input, media, s, k)
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    H = 2h_i
    y = zero(H)
    γ² = s * μ₀ * inv(media[2].rho)
    section = f.options.integral
    J = first(integrate(
        SpectralIntegral(
            λ -> f(Val(:J), λ, γ², H, y), buffers, (zero(s), sqrt(γ²)), (), (), H, y),
        section.method, section, buffers))
    return s * μ₀ / (2π) * (Λ(f, Val(:self), zero(s), r_i, h_i) + 2J)
end

function (f::Formula{:carson1926})(
        ::Val{:Z}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    h_j = abs(input.y[j])
    y = abs(input.x[i] - input.x[j])
    H = h_i + h_j
    D = hypot(y, h_i - h_j)
    D′ = hypot(y, H)
    γ² = s * μ₀ * inv(media[2].rho)
    section = f.options.integral
    J = first(integrate(
        SpectralIntegral(
            λ -> f(Val(:J), λ, γ², H, y), buffers, (zero(s), sqrt(γ²)), (), (), H, y),
        section.method, section, buffers))
    return s * μ₀ / (2π) * (Λ(f, Val(:mutual), zero(s), r_i, D, D′) + 2J)
end
