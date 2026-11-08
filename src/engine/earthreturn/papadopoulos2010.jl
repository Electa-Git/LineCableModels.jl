defaults(::Type{Formula{:papadopoulos2010}}) = (integral = integrate => :quad,)

validate(formula::Formula{:papadopoulos2010}, ::typeof(impedance!)) = formula
validate(formula::Formula{:papadopoulos2010}, ::typeof(potential!)) = formula
function validate(formula::Formula{:papadopoulos2010}, t::Int, s::Int)
    t > 1 && s > 1 || invoke(validate, Tuple{Formula, Int, Int}, formula, t, s)
    return formula
end

function initialize_buffers(f::Formula{:papadopoulos2010}, input::CoaxialInput{T},
        buffers) where {T}
    return initialize_buffers(SpectralIntegral, f.options.integral.method, T, buffers)
end

function Λ(::Formula{:papadopoulos2010}, ::Val{:self}, κ, r, h)
    besselk(0, κ * r) - besselk(0, κ * hypot(r, 2h))
end
function Λ(::Formula{:papadopoulos2010}, ::Val{:mutual}, κ, r, D, D′)
    besselk(0, κ * D) - besselk(0, κ * D′)
end

function (::Formula{:papadopoulos2010})(::Val{:J}, λ, κ₁², κ₂², μ₂, H, y)
    a₁ = sqrt(λ^2 + κ₁²)
    a₂ = sqrt(λ^2 + κ₂²)
    return μ₀ * (exp(-(H + im * y) * λ) + exp(-(H - im * y) * λ)) *
           exp(-H * κ₂² / (a₂ + λ)) / (2 * (a₂ * μ₀ + a₁ * μ₂))
end

function (::Formula{:papadopoulos2010})(::Val{:P̃}, λ, κ₁², κ₂², γ₁², γ₂², μ₂, H, y)
    a₁ = sqrt(λ^2 + κ₁²)
    a₂ = sqrt(λ^2 + κ₂²)
    return μ₀ * (exp(-(H + im * y) * λ) + exp(-(H - im * y) * λ)) *
           exp(-H * κ₂² / (a₂ + λ)) / (2 * (a₂ * μ₀ + a₁ * μ₂)) *
           (1 + μ₂ * a₂ * (γ₂² - γ₁²) / (a₂ * γ₁² * μ₂ + a₁ * γ₂² * μ₀))
end

function (f::Formula{:papadopoulos2010})(
        ::Val{:Z}, ::Val{:self}, i, buffers, input, media, s, k)
    μ₂ = buffers.μ[2]
    σ̃₂ = buffers.σ̃[2]
    γ₁² = s^2 * μ₀ * (ε₀ * media[1].eps_r)
    γ₂² = s * μ₂ * σ̃₂
    κ₁² = γ₁² - s^2 * μ₂ * (ε₀ * media[2].eps_r)
    κ₂² = γ₂² - s^2 * μ₂ * (ε₀ * media[2].eps_r)
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    H = 2h_i
    y = r_i
    section = f.options.integral
    J = first(integrate(
        SpectralIntegral(λ -> f(Val(:J), λ, κ₁², κ₂², μ₂, H, y), buffers,
            (sqrt(κ₁²), sqrt(κ₂²)), (), (), H, y),
        section.method, section, buffers))
    return s * μ₂ / (2π) * (Λ(f, Val(:self), sqrt(κ₂²), r_i, h_i) + 2J)
end

function (f::Formula{:papadopoulos2010})(
        ::Val{:Z}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    μ₂ = buffers.μ[2]
    σ̃₂ = buffers.σ̃[2]
    γ₁² = s^2 * μ₀ * (ε₀ * media[1].eps_r)
    γ₂² = s * μ₂ * σ̃₂
    κ₁² = γ₁² - s^2 * μ₂ * (ε₀ * media[2].eps_r)
    κ₂² = γ₂² - s^2 * μ₂ * (ε₀ * media[2].eps_r)
    h_i = abs(input.y[i])
    h_j = abs(input.y[j])
    y = abs(input.x[i] - input.x[j])
    H = h_i + h_j
    D = hypot(y, h_i - h_j)
    D′ = hypot(y, H)
    section = f.options.integral
    J = first(integrate(
        SpectralIntegral(λ -> f(Val(:J), λ, κ₁², κ₂², μ₂, H, y), buffers,
            (sqrt(κ₁²), sqrt(κ₂²)), (), (), H, y),
        section.method, section, buffers))
    return s * μ₂ / (2π) * (Λ(f, Val(:mutual), sqrt(κ₂²), input.radius[i], D, D′) + 2J)
end

function (f::Formula{:papadopoulos2010})(
        ::Val{:P}, ::Val{:self}, i, buffers, input, media, s, k)
    μ₂ = buffers.μ[2]
    σ̃₂ = buffers.σ̃[2]
    γ₁² = s^2 * μ₀ * (ε₀ * media[1].eps_r)
    γ₂² = s * μ₂ * σ̃₂
    κ₁² = γ₁² - s^2 * μ₂ * (ε₀ * media[2].eps_r)
    κ₂² = γ₂² - s^2 * μ₂ * (ε₀ * media[2].eps_r)
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    H = 2h_i
    y = r_i
    section = f.options.integral
    P̃ = first(integrate(
        SpectralIntegral(λ -> f(Val(:P̃), λ, κ₁², κ₂², γ₁², γ₂², μ₂, H, y), buffers,
            (sqrt(κ₁²), sqrt(κ₂²)), (), (), H, y),
        section.method, section, buffers))
    return s / (2π * σ̃₂) * (Λ(f, Val(:self), sqrt(κ₂²), r_i, h_i) + 2P̃)
end

function (f::Formula{:papadopoulos2010})(
        ::Val{:P}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    μ₂ = buffers.μ[2]
    σ̃₂ = buffers.σ̃[2]
    γ₁² = s^2 * μ₀ * (ε₀ * media[1].eps_r)
    γ₂² = s * μ₂ * σ̃₂
    κ₁² = γ₁² - s^2 * μ₂ * (ε₀ * media[2].eps_r)
    κ₂² = γ₂² - s^2 * μ₂ * (ε₀ * media[2].eps_r)
    h_i = abs(input.y[i])
    h_j = abs(input.y[j])
    y = abs(input.x[i] - input.x[j])
    H = h_i + h_j
    D = hypot(y, h_i - h_j)
    D′ = hypot(y, H)
    section = f.options.integral
    P̃ = first(integrate(
        SpectralIntegral(λ -> f(Val(:P̃), λ, κ₁², κ₂², γ₁², γ₂², μ₂, H, y), buffers,
            (sqrt(κ₁²), sqrt(κ₂²)), (), (), H, y),
        section.method, section, buffers))
    return s / (2π * σ̃₂) * (Λ(f, Val(:mutual), sqrt(κ₂²), input.radius[i], D, D′) + 2P̃)
end
