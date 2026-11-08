defaults(::Type{Formula{:xue2018}}) = (integral = integrate => :quad,)

validate(formula::Formula{:xue2018}, ::typeof(impedance!)) = formula
validate(formula::Formula{:xue2018}, ::typeof(potential!)) = formula
function validate(formula::Formula{:xue2018}, t::Int, s::Int)
    t > 1 && s > 1 || invoke(validate, Tuple{Formula, Int, Int}, formula, t, s)
    return formula
end

validate(Γ::Union{Number, AbstractVector}, formula::Formula{:xue2018}) = Γ

"""
$(TYPEDSIGNATURES)

Refuse a magnetic earth layer: the formula assumes the permeability of vacuum
on both sides of the interface.

# Arguments

- `problem`: earth layers of the computation.
- `formula`: earth-return formula `Formula{:xue2018}` that consumes the problem.

# Returns

- `problem`, unchanged.

# Errors

- Throws `DomainError` for an earth layer whose relative permeability differs from 1.
"""
function validate(problem::LineParametersProblem, formula::Formula{:xue2018})
    for (k, layer) in enumerate(Base.tail(problem.earth.layers))
        isone(layer.mu_r) || throw(DomainError(layer.mu_r,
            "the earth-return formula :xue2018 holds for a nonmagnetic earth; the earth layer $k has mu_r = $(layer.mu_r)"))
    end
    return problem
end

function initialize_buffers(f::Formula{:xue2018}, input::CoaxialInput{T},
        buffers) where {T}
    return initialize_buffers(SpectralIntegral, f.options.integral.method, T, buffers)
end

function Λ(::Formula{:xue2018}, ::Val{:self}, κ, r, h)
    besselk(0, κ * r) - besselk(0, κ * hypot(r, 2h))
end
Λ(::Formula{:xue2018}, ::Val{:mutual}, κ, r, D, D′) = besselk(0, κ * D) - besselk(0, κ * D′)

function (::Formula{:xue2018})(::Val{:J}, λ, κ₁², κ₂², H, y)
    a₁ = sqrt(λ^2 + κ₁²)
    a₂ = sqrt(λ^2 + κ₂²)
    return (exp(-(H + im * y) * λ) + exp(-(H - im * y) * λ)) *
           exp(-H * κ₂² / (a₂ + λ)) / (2 * (a₁ + a₂))
end

function (::Formula{:xue2018})(::Val{:M}, λ, κ₁², κ₂², σ̃₁, σ̃₂, H, y)
    a₁ = sqrt(λ^2 + κ₁²)
    a₂ = sqrt(λ^2 + κ₂²)
    return (exp(-(H + im * y) * λ) + exp(-(H - im * y) * λ)) *
           exp(-H * κ₂² / (a₂ + λ)) * (σ̃₂ - σ̃₁) /
           (2 * a₂ * (a₁ + a₂) * (σ̃₂ * a₁ + σ̃₁ * a₂))
end

function (::Formula{:xue2018})(::Val{:P̃}, λ, κ₁², κ₂², σ̃₁, σ̃₂, H, y)
    a₁ = sqrt(λ^2 + κ₁²)
    a₂ = sqrt(λ^2 + κ₂²)
    return (exp(-(H + im * y) * λ) + exp(-(H - im * y) * λ)) *
           exp(-H * κ₂² / (a₂ + λ)) * σ̃₂ * a₁ / (2 * a₂ * (σ̃₂ * a₁ + σ̃₁ * a₂))
end

function (f::Formula{:xue2018})(
        ::Val{:Z}, ::Val{:self}, i, buffers, input, media, s, k)
    Γ = buffers.gamma[k]
    σ̃₁ = buffers.σ̃[1]
    σ̃₂ = buffers.σ̃[2]
    κ₁² = s * μ₀ * σ̃₁ - Γ^2
    κ₂² = s * μ₀ * σ̃₂ - Γ^2
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    H = 2h_i
    y = r_i
    section = f.options.integral
    J = first(integrate(
        SpectralIntegral(λ -> f(Val(:J), λ, κ₁², κ₂², H, y), buffers,
            (sqrt(κ₁²), sqrt(κ₂²)), (), (abs(sqrt(κ₂²) * σ̃₁ / σ̃₂),), H, y),
        section.method, section, buffers))
    M = iszero(Γ) ? zero(J) :
        first(integrate(
        SpectralIntegral(λ -> f(Val(:M), λ, κ₁², κ₂², σ̃₁, σ̃₂, H, y), buffers,
            (sqrt(κ₁²), sqrt(κ₂²)), (), (abs(sqrt(κ₂²) * σ̃₁ / σ̃₂),), H, y),
        section.method, section, buffers))
    return s * μ₀ / (2π) * (Λ(f, Val(:self), sqrt(κ₂²), r_i, h_i) + 2J - 2Γ^2 * M)
end

function (f::Formula{:xue2018})(
        ::Val{:Z}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    Γ = buffers.gamma[k]
    σ̃₁ = buffers.σ̃[1]
    σ̃₂ = buffers.σ̃[2]
    κ₁² = s * μ₀ * σ̃₁ - Γ^2
    κ₂² = s * μ₀ * σ̃₂ - Γ^2
    h_i = abs(input.y[i])
    h_j = abs(input.y[j])
    y = abs(input.x[i] - input.x[j])
    H = h_i + h_j
    D = hypot(y, h_i - h_j)
    D′ = hypot(y, H)
    section = f.options.integral
    J = first(integrate(
        SpectralIntegral(λ -> f(Val(:J), λ, κ₁², κ₂², H, y), buffers,
            (sqrt(κ₁²), sqrt(κ₂²)), (), (abs(sqrt(κ₂²) * σ̃₁ / σ̃₂),), H, y),
        section.method, section, buffers))
    M = iszero(Γ) ? zero(J) :
        first(integrate(
        SpectralIntegral(λ -> f(Val(:M), λ, κ₁², κ₂², σ̃₁, σ̃₂, H, y), buffers,
            (sqrt(κ₁²), sqrt(κ₂²)), (), (abs(sqrt(κ₂²) * σ̃₁ / σ̃₂),), H, y),
        section.method, section, buffers))
    return s * μ₀ / (2π) *
           (Λ(f, Val(:mutual), sqrt(κ₂²), input.radius[i], D, D′) + 2J - 2Γ^2 * M)
end

function (f::Formula{:xue2018})(
        ::Val{:P}, ::Val{:self}, i, buffers, input, media, s, k)
    Γ = buffers.gamma[k]
    σ̃₁ = buffers.σ̃[1]
    σ̃₂ = buffers.σ̃[2]
    κ₁² = s * μ₀ * σ̃₁ - Γ^2
    κ₂² = s * μ₀ * σ̃₂ - Γ^2
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    H = 2h_i
    y = r_i
    section = f.options.integral
    P̃ = first(integrate(
        SpectralIntegral(λ -> f(Val(:P̃), λ, κ₁², κ₂², σ̃₁, σ̃₂, H, y), buffers,
            (sqrt(κ₁²), sqrt(κ₂²)), (), (abs(sqrt(κ₂²) * σ̃₁ / σ̃₂),), H, y),
        section.method, section, buffers))
    return s / (2π * σ̃₂) * (Λ(f, Val(:self), sqrt(κ₂²), r_i, h_i) + 2P̃)
end

function (f::Formula{:xue2018})(
        ::Val{:P}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    Γ = buffers.gamma[k]
    σ̃₁ = buffers.σ̃[1]
    σ̃₂ = buffers.σ̃[2]
    κ₁² = s * μ₀ * σ̃₁ - Γ^2
    κ₂² = s * μ₀ * σ̃₂ - Γ^2
    h_i = abs(input.y[i])
    h_j = abs(input.y[j])
    y = abs(input.x[i] - input.x[j])
    H = h_i + h_j
    D = hypot(y, h_i - h_j)
    D′ = hypot(y, H)
    section = f.options.integral
    P̃ = first(integrate(
        SpectralIntegral(λ -> f(Val(:P̃), λ, κ₁², κ₂², σ̃₁, σ̃₂, H, y), buffers,
            (sqrt(κ₁²), sqrt(κ₂²)), (), (abs(sqrt(κ₂²) * σ̃₁ / σ̃₂),), H, y),
        section.method, section, buffers))
    return s / (2π * σ̃₂) * (Λ(f, Val(:mutual), sqrt(κ₂²), input.radius[i], D, D′) + 2P̃)
end
