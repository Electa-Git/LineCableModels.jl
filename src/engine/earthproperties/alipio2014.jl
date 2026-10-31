function defaults(::Type{Formula{:alipio2014}})
    (exponent = 0.54, epsilon_high = 12.0, scale = 1.26, sigma_exponent = -0.73)
end

validate(formula::Formula{:alipio2014}, ::Type{EarthLayer}) = formula

function (formula::Formula{:alipio2014})(layer::EarthLayer{T}, f::T) where {T <: Real}
    γ = convert(T, formula.options.exponent)
    ε_high = convert(T, formula.options.epsilon_high)
    a = convert(T, formula.options.scale)
    b = convert(T, formula.options.sigma_exponent)
    σ_low = convert(T, 1000) / layer.rho
    D = a * σ_low^b
    ε_r = ε_high +
          tan(convert(T, π) * γ / convert(T, 2)) *
          convert(T, 1e-3) * σ_low * D * f^(γ - one(T)) /
          (convert(T, 2) * convert(T, π) * convert(T, vacuum_permittivity) *
           convert(T, 10)^(convert(T, 6) * γ))
    σ = σ_low * (one(T) + D * (f / convert(T, 1e6))^γ)
    ρ = convert(T, 1000) / σ
    return EarthLayer{T}(ρ, ε_r, layer.mu_r, layer.thickness)
end

"""
$(TYPEDSIGNATURES)

Check the options of the law.

# Arguments

- `options`: options of `Formula{:alipio2014}`.

# Returns

- `options`, unchanged.

# Errors

- Throws `ArgumentError` for an odd integer exponent.
"""
function validate(options::FormulationOptions{Formula{:alipio2014}})
    exponent = options.exponent
    isinteger(exponent) && isodd(exponent) &&
        throw(ArgumentError(
            "Alipio–Visacro exponent must not be an odd integer (tangent pole)"))
    return options
end

"""
$(TYPEDSIGNATURES)

Warn once when the scan of `problem` leaves the range of the law.

# Arguments

- `problem`: frequencies of the scan.
- `formula`: law `Formula{:alipio2014}` that consumes the problem.

# Returns

- `problem`, unchanged.
"""
function validate(problem::LineParametersProblem, formula::Formula{:alipio2014})
    frequencies = problem.frequencies
    if first(frequencies) < 100 || last(frequencies) > 4e6
        @warn "the law :alipio2014 holds from 100 Hz to 4 MHz; the scan runs from $(first(frequencies)) Hz to $(last(frequencies)) Hz"
    end
    return problem
end
