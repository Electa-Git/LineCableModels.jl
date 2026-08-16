
struct Lossless <: InsulationImpedanceFormulation end
get_description(::Lossless) = "Lossless insulation (ideal dielectric)"

@inline function (f::Lossless)(
        r_in::T,
        r_ex::T,
        mur_i::T,
        jω::Complex{T}
) where {T <: REALSCALAR}
    if isapprox(r_in, 0.0, atol = eps(T)) || isapprox(r_in, r_ex, atol = eps(T))
        # TODO: Implement consistent handling of admittance for bare conductors
        # Issue URL: https://github.com/Electa-Git/LineCableModels.jl/issues/18
        return zero(Complex{T})
    end

    # Constants
    mu_i = _typed_μ₀(T) * mur_i

    return Complex{T}(jω * mu_i * log(r_ex / r_in) / (T(2) * T(π)))
end
