
struct Lossless <: InsulationImpedanceFormulation end
description(::Lossless) = "Lossless insulation (ideal dielectric)"

@inline function (f::Lossless)(
        r_in::T,
        r_ex::T,
        mur_i::T,
        jω::Complex{T}
) where {T <: Real}
    if isapprox(r_in, 0.0, atol = eps(T)) || isapprox(r_in, r_ex, atol = eps(T))
        # TODO: Implement consistent handling of admittance for bare conductors
        # Issue URL: https://github.com/Electa-Git/LineCableModels.jl/issues/18
        return zero(Complex{T})
    end

    mu0 = one(r_in) * 4 * (one(r_in) * π) * (one(r_in) * 10)^(-7)
    mu_i = mu0 * mur_i

    return Complex{T}(jω * mu_i * log(r_ex / r_in) / (2 * (one(r_in) * π)))
end
