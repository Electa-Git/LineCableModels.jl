function fem_lossless_semicon(
        material::Material{T}, frequency, temperature, assumptions
) where {T <: Real}
    ε₀ = one(T) * 88541878128 * (one(T) * 10)^(-22)
    ω = 2 * (one(T) * π) * convert(T, frequency)
    return complex(zero(T), ω) * ε₀ * material.eps_r
end

const FEM_LOSSLESS_SEMICON = SEMICON_ADMITTANCE.Formula(
    :FEMLossless,
    fem_lossless_semicon
)
