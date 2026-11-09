defaults(::Type{Formula{:portela1999}}) = (beta = 0.1, exponent = 0.72)
domain(::Formula{:portela1999}) = (frequencies = 100.0 .. 2.0e6,)

validate(formula::Formula{:portela1999}, ::Type{EarthLayer}) = formula

function (formula::Formula{:portela1999})(layer::EarthLayer{T}, f::T) where {T <: Real}
    β = convert(T, formula.options.beta)
    α = convert(T, formula.options.exponent)
    σ_low = inv(layer.rho)
    ω = convert(T, 2) * convert(T, π) * f
    σ = σ_low + β * convert(T, 1e-6) * ω^α
    ε_r = β * convert(T, 1e-6) * tan(convert(T, π) * α / convert(T, 2)) *
          ω^(α - one(T)) / convert(T, ε₀)
    ρ = inv(σ)
    return EarthLayer{T}(ρ, ε_r, layer.mu_r, layer.thickness)
end

"""
$(TYPEDSIGNATURES)

Check the options of the law.

# Arguments

- `options`: options of `Formula{:portela1999}`.

# Returns

- `options`, unchanged.

# Errors

- Throws `ArgumentError` for an odd integer exponent.
"""
function validate(options::FormulationOptions{Formula{:portela1999}})
    exponent = options.exponent
    isinteger(exponent) && isodd(exponent) &&
        throw(ArgumentError(
            "Portela exponent must not be an odd integer (tangent pole)"))
    return options
end
