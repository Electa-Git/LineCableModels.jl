function defaults(::Type{Formula{:longmire1975}})
    (epsilon_high = 5.0, corner_scale = 125.0, corner_exponent = 0.8312)
end

validate(formula::Formula{:longmire1975}, ::Type{EarthLayer}) = formula

function (formula::Formula{:longmire1975})(layer::EarthLayer{T}, f::T) where {T <: Real}
    ε_high = convert(T, formula.options.epsilon_high)
    a = convert(T, formula.options.corner_scale)
    b = convert(T, formula.options.corner_exponent)
    σ_low = inv(layer.rho)
    f_c = (a * σ_low)^b
    Σε = zero(T)
    Σσ = zero(T)
    decade = one(T)
    for a_n in (convert(T, 3.4e6), convert(T, 2.74e5), convert(T, 2.58e4),
        convert(T, 3.38e3), convert(T, 5.26e2), convert(T, 1.33e2),
        convert(T, 2.72e1), convert(T, 1.25e1), convert(T, 4.8),
        convert(T, 2.17), convert(T, 0.98), convert(T, 0.392), convert(T, 0.173))
        f_n = f_c * decade
        Σε += a_n / (one(T) + (f / f_n)^convert(T, 2))
        Σσ += a_n * (f / f_n) / (one(T) + (f / f_n)^convert(T, 2))
        decade *= convert(T, 10)
    end
    ε_r = ε_high + Σε
    σ = σ_low + convert(T, 2) * convert(T, π) * f * convert(T, ε₀) * Σσ
    ρ = inv(σ)
    return EarthLayer{T}(ρ, ε_r, layer.mu_r, layer.thickness)
end

"""
$(TYPEDSIGNATURES)

Check the options of the law.

# Arguments

- `options`: options of `Formula{:longmire1975}`.

# Returns

- `options`, unchanged.

# Errors

- Throws `ArgumentError` for a nonpositive `corner_scale`.
"""
function validate(options::FormulationOptions{Formula{:longmire1975}})
    options.corner_scale > 0 || throw(ArgumentError(
        "Longmire–Smith corner_scale must be positive"))
    return options
end

"""
$(TYPEDSIGNATURES)

Warn once when the scan of `problem` leaves the range of the law.

# Arguments

- `problem`: frequencies of the scan.
- `formula`: law `Formula{:longmire1975}` that consumes the problem.

# Returns

- `problem`, unchanged.
"""
function validate(problem::LineParametersProblem, formula::Formula{:longmire1975})
    frequencies = problem.frequencies
    if first(frequencies) < 5 || last(frequencies) > 3e10
        @warn "the law :longmire1975 holds from 5 Hz to 30 GHz; the scan runs from $(first(frequencies)) Hz to $(last(frequencies)) Hz"
    end
    return problem
end
