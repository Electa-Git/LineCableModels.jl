defaults(::Type{Formula{:messier1985}}) = (epsilon_high = 8.0,)
domain(::Formula{:messier1985}) = (frequencies = 100.0 .. 1.0e6,)

validate(formula::Formula{:messier1985}, ::Type{EarthLayer}) = formula

function (formula::Formula{:messier1985})(layer::EarthLayer{T}, f::T) where {T <: Real}
    ε_high = convert(T, formula.options.epsilon_high)
    σ_low = inv(layer.rho)
    ε_r = ε_high +
          sqrt(σ_low * ε_high / (convert(T, π) * f * convert(T, ε₀)))
    σ = σ_low + sqrt(convert(T, 4) * convert(T, π) * f * σ_low *
                     convert(T, ε₀) * ε_high)
    ρ = inv(σ)
    return EarthLayer{T}(ρ, ε_r, layer.mu_r, layer.thickness)
end

"""
$(TYPEDSIGNATURES)

Check the options of the law.

# Arguments

- `options`: options of `Formula{:messier1985}`.

# Returns

- `options`, unchanged.

# Errors

- Throws `ArgumentError` for a negative `epsilon_high`.
"""
function validate(options::FormulationOptions{Formula{:messier1985}})
    options.epsilon_high >= 0 || throw(ArgumentError(
        "Messier epsilon_high must be nonnegative for the real square roots"))
    return options
end
