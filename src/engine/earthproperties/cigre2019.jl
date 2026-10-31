function defaults(::Type{Formula{:cigre2019}})
    (epsilon_high = 12.0, epsilon_scale = 9.5e4,
        epsilon_sigma_exponent = 0.27, epsilon_frequency_exponent = -0.46,
        sigma_scale = 4.7e-6, sigma_frequency_exponent = 0.54)
end

validate(formula::Formula{:cigre2019}, ::Type{EarthLayer}) = formula

function (formula::Formula{:cigre2019})(layer::EarthLayer{T}, f::T) where {T <: Real}
    ε_high = convert(T, formula.options.epsilon_high)
    a = convert(T, formula.options.epsilon_scale)
    b = convert(T, formula.options.epsilon_sigma_exponent)
    c = convert(T, formula.options.epsilon_frequency_exponent)
    d = convert(T, formula.options.sigma_scale)
    e = convert(T, formula.options.sigma_frequency_exponent)
    σ_low = inv(layer.rho)
    ε_r = ε_high + a * σ_low^b * f^c
    σ = σ_low + d * σ_low^b * f^e
    ρ = inv(σ)
    return EarthLayer{T}(ρ, ε_r, layer.mu_r, layer.thickness)
end

"""
$(TYPEDSIGNATURES)

Warn once when the scan of `problem` leaves the range of the law.

# Arguments

- `problem`: frequencies of the scan.
- `formula`: law `Formula{:cigre2019}` that consumes the problem.

# Returns

- `problem`, unchanged.
"""
function validate(problem::LineParametersProblem, formula::Formula{:cigre2019})
    frequencies = problem.frequencies
    if first(frequencies) < 100 || last(frequencies) > 4e6
        @warn "the law :cigre2019 holds from 100 Hz to 4 MHz; the scan runs from $(first(frequencies)) Hz to $(last(frequencies)) Hz"
    end
    return problem
end
