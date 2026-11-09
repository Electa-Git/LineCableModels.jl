defaults(::Type{Formula{:visacro1987}}) = (normalization_frequency = 100.0,)
domain(::Formula{:visacro1987}) = (frequencies = 40.0 .. 2.0e6,)

validate(formula::Formula{:visacro1987}, ::Type{EarthLayer}) = formula

function (formula::Formula{:visacro1987})(layer::EarthLayer{T}, f::T) where {T <: Real}
    f_0 = convert(T, formula.options.normalization_frequency)
    σ_low = inv(layer.rho)
    ε_r = convert(T, 2.34e6) * σ_low^convert(T, 0.535) * f^convert(T, -0.597)
    σ = σ_low * (f / f_0)^convert(T, 0.072)
    ρ = inv(σ)
    return EarthLayer{T}(ρ, ε_r, layer.mu_r, layer.thickness)
end

"""
$(TYPEDSIGNATURES)

Check the options of the law.

# Arguments

- `options`: options of `Formula{:visacro1987}`.

# Returns

- `options`, unchanged.

# Errors

- Throws `ArgumentError` for a nonpositive `normalization_frequency`.
"""
function validate(options::FormulationOptions{Formula{:visacro1987}})
    options.normalization_frequency > 0 || throw(ArgumentError(
        "Visacro–Portela normalization_frequency must be positive [Hz]"))
    return options
end
