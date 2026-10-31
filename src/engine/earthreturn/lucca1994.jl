defaults(::Type{Formula{:lucca1994}}) = (;)

validate(formula::Formula{:lucca1994}, ::typeof(impedance!)) = formula
validate(formula::Formula{:lucca1994}, ::Val{1}, ::Val{2}) = formula
validate(formula::Formula{:lucca1994}, ::Val{2}, ::Val{1}) = formula

"""
$(TYPEDSIGNATURES)

Refuse a magnetic earth layer: the formula assumes the permeability of vacuum
on both sides of the interface.

# Arguments

- `problem`: earth layers of the computation.
- `formula`: earth-return formula `Formula{:lucca1994}` that consumes the problem.

# Returns

- `problem`, unchanged.

# Errors

- Throws `DomainError` for an earth layer whose relative permeability differs from 1.
"""
function validate(problem::LineParametersProblem, formula::Formula{:lucca1994})
    for (k, layer) in enumerate(Base.tail(problem.earth.layers))
        isone(layer.mu_r) || throw(DomainError(layer.mu_r,
            "the earth-return formula :lucca1994 holds for a nonmagnetic earth; the earth layer $k has mu_r = $(layer.mu_r)"))
    end
    return problem
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
