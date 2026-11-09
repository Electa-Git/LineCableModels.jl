defaults(::Type{Formula{:visacro2012}}) = (frequency_boundary = 100.0,)
domain(::Formula{:visacro2012}) = (frequencies = 100.0 .. 4.0e6,)

validate(formula::Formula{:visacro2012}, ::Type{EarthLayer}) = formula

function (formula::Formula{:visacro2012})(layer::EarthLayer{T}, f::T) where {T <: Real}
    f_0 = convert(T, formula.options.frequency_boundary)
    f_e = max(f, f_0)
    σ_low = inv(layer.rho)
    ε_r = convert(T, 1.3) + convert(T, 7.6e3) * f_e^convert(T, -0.4)
    σ = σ_low + convert(T, 1.2e-6) * σ_low^convert(T, 0.27) *
                (f_e - f_0)^convert(T, 0.65)
    ρ = inv(σ)
    return EarthLayer{T}(ρ, ε_r, layer.mu_r, layer.thickness)
end

"""
$(TYPEDSIGNATURES)

Check the options of the law.

# Arguments

- `options`: options of `Formula{:visacro2012}`.

# Returns

- `options`, unchanged.

# Errors

- Throws `ArgumentError` for a nonpositive `frequency_boundary`.
"""
function validate(options::FormulationOptions{Formula{:visacro2012}})
    options.frequency_boundary > 0 || throw(ArgumentError(
        "Visacro–Alipio frequency_boundary must be positive [Hz]"))
    return options
end
