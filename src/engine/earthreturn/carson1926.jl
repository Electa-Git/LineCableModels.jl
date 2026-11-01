defaults(::Type{Formula{:carson1926}}) = (integral = integrate => :quad,)

function validate(options::FormulationOptions{Formula{:carson1926}})
    validate(options.integral, integrate, options.integral.method)
    return options
end

validate(formula::Formula{:carson1926}, ::typeof(impedance!)) = formula
validate(formula::Formula{:carson1926}, ::Val{1}, ::Val{1}) = formula

"""
$(TYPEDSIGNATURES)

Refuse a magnetic earth layer: the formula assumes the permeability of vacuum
on both sides of the interface.

# Arguments

- `problem`: earth layers of the computation.
- `formula`: earth-return formula `Formula{:carson1926}` that consumes the problem.

# Returns

- `problem`, unchanged.

# Errors

- Throws `DomainError` for an earth layer whose relative permeability differs from 1.
"""
function validate(problem::LineParametersProblem, formula::Formula{:carson1926})
    for (k, layer) in enumerate(Base.tail(problem.earth.layers))
        isone(layer.mu_r) || throw(DomainError(layer.mu_r,
            "the earth-return formula :carson1926 holds for a nonmagnetic earth; the earth layer $k has mu_r = $(layer.mu_r)"))
    end
    return problem
end

function initialize_buffers(f::Formula{:carson1926}, ::LineParametersProblem,
        input::CoaxialInput{T}, buffers) where {T}
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
