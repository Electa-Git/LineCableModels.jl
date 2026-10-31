defaults(::Type{Formula{:ametani2009}}) = (;)

validate(formula::Formula{:ametani2009}, ::typeof(impedance!)) = formula
validate(formula::Formula{:ametani2009}, ::Val{1}, ::Val{2}) = formula
validate(formula::Formula{:ametani2009}, ::Val{2}, ::Val{1}) = formula

"""
$(TYPEDSIGNATURES)

Refuse a magnetic earth layer: the formula assumes the permeability of vacuum
on both sides of the interface.

# Arguments

- `problem`: earth layers of the computation.
- `formula`: earth-return formula `Formula{:ametani2009}` that consumes the problem.

# Returns

- `problem`, unchanged.

# Errors

- Throws `DomainError` for an earth layer whose relative permeability differs from 1.
"""
function validate(problem::LineParametersProblem, formula::Formula{:ametani2009})
    for (k, layer) in enumerate(Base.tail(problem.earth.layers))
        isone(layer.mu_r) || throw(DomainError(layer.mu_r,
            "the earth-return formula :ametani2009 holds for a nonmagnetic earth; the earth layer $k has mu_r = $(layer.mu_r)"))
    end
    return problem
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
