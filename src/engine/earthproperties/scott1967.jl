defaults(::Type{Formula{:scott1967}}) = (;)
domain(::Formula{:scott1967}) = (frequencies = 100.0 .. 1.0e6,)

validate(formula::Formula{:scott1967}, ::Type{EarthLayer}) = formula

function (formula::Formula{:scott1967})(layer::EarthLayer{T}, f::T) where {T <: Real}
    u = log10(convert(T, 1000) / layer.rho)
    x = log10(f)
    ε_r = convert(T,
        10)^(convert(T, 5.491) + convert(T, 0.946) * u -
             convert(T, 1.097) * x + convert(T, 0.069) * u^convert(T, 2) -
             convert(T, 0.114) * u * x + convert(T, 0.067) * x^convert(T, 2))
    σ = convert(T,
        10)^(convert(T, 0.028) + convert(T, 1.098) * u -
             convert(T, 0.068) * x + convert(T, 0.036) * u^convert(T, 2) -
             convert(T, 0.046) * u * x + convert(T, 0.018) * x^convert(T, 2))
    ρ = convert(T, 1000) / σ
    return EarthLayer{T}(ρ, ε_r, layer.mu_r, layer.thickness)
end
