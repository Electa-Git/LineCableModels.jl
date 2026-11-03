defaults(::Type{Formula{:pollaczek1926}}) = (integral = integrate => :quad,)

function validate(options::FormulationOptions{Formula{:pollaczek1926}})
    validate(options.integral, integrate, options.integral.method)
    return options
end

validate(formula::Formula{:pollaczek1926}, ::typeof(impedance!)) = formula
function validate(formula::Formula{:pollaczek1926}, t::Int, s::Int)
    t > 1 || s > 1 || invoke(validate, Tuple{Formula, Int, Int}, formula, t, s)
    return formula
end

"""
$(TYPEDSIGNATURES)

Refuse a magnetic earth layer: the formula assumes the permeability of vacuum
on both sides of the interface.

# Arguments

- `problem`: earth layers of the computation.
- `formula`: earth-return formula `Formula{:pollaczek1926}` that consumes the problem.

# Returns

- `problem`, unchanged.

# Errors

- Throws `DomainError` for an earth layer whose relative permeability differs from 1.
"""
function validate(problem::LineParametersProblem, formula::Formula{:pollaczek1926})
    for (k, layer) in enumerate(Base.tail(problem.earth.layers))
        isone(layer.mu_r) || throw(DomainError(layer.mu_r,
            "the earth-return formula :pollaczek1926 holds for a nonmagnetic earth; the earth layer $k has mu_r = $(layer.mu_r)"))
    end
    return problem
end

function initialize_buffers(f::Formula{:pollaczek1926}, input::CoaxialInput{T},
        buffers) where {T}
    return initialize_buffers(SpectralIntegral, f.options.integral.method, T, buffers)
end

Λ(::Formula{:pollaczek1926}, ::Val{:self}, κ, r, h) = besselk(0, κ * r) - besselk(0, 2κ * h)
function Λ(::Formula{:pollaczek1926}, ::Val{:mutual}, κ, r, D, D′)
    besselk(0, κ * D) - besselk(0, κ * D′)
end

function (::Formula{:pollaczek1926})(::Val{:J}, ::Val{:earth}, λ, γ², H, y)
    a = sqrt(λ^2 + γ²)
    return (exp(-(H + im * y) * λ) + exp(-(H - im * y) * λ)) * exp(-H * γ² / (a + λ)) /
           (2 * (λ + a))
end

function (::Formula{:pollaczek1926})(::Val{:J}, ::Val{:mixed}, λ, γ², h_g, H, y)
    a = sqrt(λ^2 + γ²)
    return (exp(-(H + im * y) * λ) + exp(-(H - im * y) * λ)) * exp(-h_g * γ² / (a + λ)) /
           (2 * (λ + a))
end

function (f::Formula{:pollaczek1926})(
        ::Val{:Z}, ::Val{:self}, i, buffers, input, media, s, k)
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    H = 2h_i
    y = zero(H)
    γ² = s * μ₀ * inv(media[2].rho)
    section = f.options.integral
    J = first(integrate(
        SpectralIntegral(
            λ -> f(Val(:J), Val(:earth), λ, γ², H, y), buffers,
            (zero(s), sqrt(γ²)), (), (), H, y),
        section.method, section, buffers))
    return s * μ₀ / (2π) * (Λ(f, Val(:self), sqrt(γ²), r_i, h_i) + 2J)
end

function (f::Formula{:pollaczek1926})(
        ::Val{:Z}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    (input.earth_layer[i] == 1) == (input.earth_layer[j] == 1) &&
        return f(Val(:Z), Val(:earth), i, j, buffers, input, media, s)
    return f(Val(:Z), Val(:mixed), i, j, buffers, input, media, s)
end

function (f::Formula{:pollaczek1926})(
        ::Val{:Z}, ::Val{:earth}, i, j, buffers, input, media, s)
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
            λ -> f(Val(:J), Val(:earth), λ, γ², H, y), buffers,
            (zero(s), sqrt(γ²)), (), (), H, y),
        section.method, section, buffers))
    return s * μ₀ / (2π) * (Λ(f, Val(:mutual), sqrt(γ²), input.radius[i], D, D′) + 2J)
end

function (f::Formula{:pollaczek1926})(
        ::Val{:Z}, ::Val{:mixed}, i, j, buffers, input, media, s)
    h_a = max(input.y[i], input.y[j])
    h_g = -min(input.y[i], input.y[j])
    y = abs(input.x[i] - input.x[j])
    H = h_a + h_g
    γ² = s * μ₀ * inv(media[2].rho)
    section = f.options.integral
    J = first(integrate(
        SpectralIntegral(
            λ -> f(Val(:J), Val(:mixed), λ, γ², h_g, H, y), buffers,
            (zero(s), sqrt(γ²)), (), (), H, y),
        section.method, section, buffers))
    return s * μ₀ / (2π) * 2J
end
