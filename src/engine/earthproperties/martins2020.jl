defaults(::Type{Formula{:martins2020}}) = (;)

validate(formula::Formula{:martins2020}, ::Type{EarthModel}) = formula

function (::Formula{:martins2020})(layers::Tuple{Vararg{EarthLayer{T}}}, f::T) where {T <:
                                                                                      Real}
    σ_eq = inv(last(layers).rho)
    for k in (length(layers) - 1):-1:1
        σ_k = inv(layers[k].rho)
        μ_k = convert(T, μ₀) * layers[k].mu_r
        h_k = layers[k].thickness
        decay = exp(-convert(T, 2) * h_k * sqrt(convert(T, π) * f * μ_k * σ_k))
        σ_eq = σ_k *
               ((sqrt(σ_k) + sqrt(σ_eq) - (sqrt(σ_k) - sqrt(σ_eq)) * decay) /
                (sqrt(σ_k) + sqrt(σ_eq) + (sqrt(σ_k) - sqrt(σ_eq)) * decay))^convert(T, 2)
    end
    return EquivalentEarthLayer{T}(inv(σ_eq), last(layers).eps_r, last(layers).mu_r)
end
