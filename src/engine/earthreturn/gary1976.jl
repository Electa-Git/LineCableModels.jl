defaults(::Type{Formula{:gary1976}}) = (;)

validate(formula::Formula{:gary1976}, ::typeof(impedance!)) = formula
validate(formula::Formula{:gary1976}, ::Val{1}, ::Val{1}) = formula

"""
$(TYPEDSIGNATURES)

Refuse a magnetic earth layer: the formula assumes the permeability of vacuum
on both sides of the interface.

# Arguments

- `problem`: earth layers of the computation.
- `formula`: earth-return formula `Formula{:gary1976}` that consumes the problem.

# Returns

- `problem`, unchanged.

# Errors

- Throws `DomainError` for an earth layer whose relative permeability differs from 1.
"""
function validate(problem::LineParametersProblem, formula::Formula{:gary1976})
    for (k, layer) in enumerate(Base.tail(problem.earth.layers))
        isone(layer.mu_r) || throw(DomainError(layer.mu_r,
            "the earth-return formula :gary1976 holds for a nonmagnetic earth; the earth layer $k has mu_r = $(layer.mu_r)"))
    end
    return problem
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
