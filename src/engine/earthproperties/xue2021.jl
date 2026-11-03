defaults(::Type{Formula{:xue2021}}) = (;)

validate(formula::Formula{:xue2021}, ::Type{EarthModel}) = formula

function (::Formula{:xue2021})(layers::Tuple{Vararg{EarthLayer{T}}}, f::T) where {T <: Real}
    μ = convert(T, μ₀)
    ε = convert(T, ε₀)
    ω = convert(T, 2) * convert(T, π) * f
    s = complex(zero(T), ω)
    γ = map(layer -> sqrt(s * μ * (inv(layer.rho) + s * ε * layer.eps_r) - s * s * μ * ε), layers)
    γ_eq = last(γ)
    for k in (length(layers) - 1):-1:1
        γ_k = γ[k]
        h_k = layers[k].thickness
        decay = exp(-convert(T, 2) * h_k * γ_k)
        γ_eq = γ_k *
               (γ_k + γ_eq - (γ_k - γ_eq) * decay) /
               (γ_k + γ_eq + (γ_k - γ_eq) * decay)
    end
    σ_eq = imag(γ_eq * γ_eq) / (ω * μ)
    ε_eq = one(T) - real(γ_eq * γ_eq) / (ω * ω * μ * ε)
    return EarthLayer{T}(inv(σ_eq), ε_eq, last(layers).mu_r, convert(T, Inf))
end

"""
$(TYPEDSIGNATURES)

Refuse a magnetic earth layer: the recursion of the formula assumes a nonmagnetic earth.

# Arguments

- `problem`: earth layers of the computation.
- `formula`: model `Formula{:xue2021}` that consumes the problem.

# Returns

- `problem`, unchanged.

# Errors

- Throws `DomainError` for an earth layer whose relative permeability differs from 1.
"""
function validate(problem::LineParametersProblem, formula::Formula{:xue2021})
    for (k, layer) in enumerate(Base.tail(problem.earth.layers))
        isone(layer.mu_r) || throw(DomainError(layer.mu_r,
            "the equivalent earth :xue2021 holds for a nonmagnetic earth; the earth layer $k has mu_r = $(layer.mu_r)"))
    end
    return problem
end
