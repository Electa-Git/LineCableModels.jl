defaults(::Type{Formula{:datsios2019}}) = (epsilon_dry = 3.5,)

validate(formula::Formula{:datsios2019}, ::Type{EarthLayer}) = formula

function (formula::Formula{:datsios2019})(layer::EarthLayer{T}, f::T) where {T <: Real}
    ε_dry = convert(T, formula.options.epsilon_dry)
    σ_low = convert(T, 1e4) / layer.rho
    p = convert(T, 0.537) * σ_low^convert(T, 0.16)
    ε_high = convert(T, 1.24) * σ_low^convert(T, 0.415) * ε_dry
    ε_3k = convert(T, 4) * σ_low^convert(T, 0.463) *
           (convert(T, 2.9) * ε_dry - convert(T, 3.8))
    ε_r = ε_high + (convert(T, 3000) / max(f, convert(T, 3000)))^p * (ε_3k - ε_high)
    σ_high = σ_low * (one(T) + convert(T, 0.65) / σ_low^convert(T, 0.57))
    σ = σ_high * f * convert(T, 1e-6) +
        (convert(T, 42) - convert(T, 42) * (f - convert(T, 42)) * convert(T, 1e-6)) *
        (σ_low / convert(T, 42) - σ_high * convert(T, 1e-6))
    ρ = convert(T, 1e4) / σ
    return EarthLayer{T}(ρ, ε_r, layer.mu_r, layer.thickness)
end

"""
$(TYPEDSIGNATURES)

Warn once when the scan of `problem` leaves the range of the law.

# Arguments

- `problem`: frequencies of the scan.
- `formula`: law `Formula{:datsios2019}` that consumes the problem.

# Returns

- `problem`, unchanged.
"""
function validate(problem::LineParametersProblem, formula::Formula{:datsios2019})
    frequencies = problem.frequencies
    if first(frequencies) < 42 || last(frequencies) > 1e6
        @warn "the law :datsios2019 holds from 42 Hz to 1 MHz; the scan runs from $(first(frequencies)) Hz to $(last(frequencies)) Hz"
    end
    return problem
end
