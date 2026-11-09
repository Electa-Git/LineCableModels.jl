defaults(::Type{Formula{:xue2021}}) = (;)
domain(::Formula{:xue2021}) = (mu_r_earth = 1.0 .. 1.0,)

validate(formula::Formula{:xue2021}, ::Type{EarthModel}) = formula

function (::Formula{:xue2021})(layers::Tuple{Vararg{EarthLayer{T}}}, f::T) where {T <: Real}
    μ = convert(T, μ₀)
    ε = convert(T, ε₀)
    ω = convert(T, 2) * convert(T, π) * f
    s = complex(zero(T), ω)
    γ = map(layer -> sqrt(s * μ * admittivity(layer, s) - s * s * μ * ε), layers)
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
