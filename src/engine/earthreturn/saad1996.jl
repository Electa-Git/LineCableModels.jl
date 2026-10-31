defaults(::Type{Formula{:saad1996}}) = (;)

validate(formula::Formula{:saad1996}, ::typeof(impedance!)) = formula
validate(formula::Formula{:saad1996}, ::Val{2}, ::Val{2}) = formula

"""
$(TYPEDSIGNATURES)

Refuse a magnetic earth layer: the formula assumes the permeability of vacuum
on both sides of the interface.

# Arguments

- `problem`: earth layers of the computation.
- `formula`: earth-return formula `Formula{:saad1996}` that consumes the problem.

# Returns

- `problem`, unchanged.

# Errors

- Throws `DomainError` for an earth layer whose relative permeability differs from 1.
"""
function validate(problem::LineParametersProblem, formula::Formula{:saad1996})
    for (k, layer) in enumerate(Base.tail(problem.earth.layers))
        isone(layer.mu_r) || throw(DomainError(layer.mu_r,
            "the earth-return formula :saad1996 holds for a nonmagnetic earth; the earth layer $k has mu_r = $(layer.mu_r)"))
    end
    return problem
end

function (f::Formula{:saad1996})(::Val{:Z}, ::Val{:self}, i, buffers, input, media, s, k)
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    m = sqrt(s * μ₀ / media[2].rho)
    return s * μ₀ / (2π) *
           (besselk(0, m * r_i) + 2exp(-2m * h_i) / (4 + m^2 * r_i^2))
end

function (f::Formula{:saad1996})(
        ::Val{:Z}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    h_i = abs(input.y[i])
    h_j = abs(input.y[j])
    y = abs(input.x[i] - input.x[j])
    H = h_i + h_j
    D = hypot(y, h_i - h_j)
    m = sqrt(s * μ₀ / media[2].rho)
    return s * μ₀ / (2π) *
           (besselk(0, m * D) + 2exp(-m * H) / (4 + m^2 * y^2))
end
