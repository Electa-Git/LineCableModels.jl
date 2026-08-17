struct Lossless <: InsulationAdmittanceFormulation end
description(::Lossless) = "Lossless insulation (ideal dielectric)"

@inline function (f::Lossless)(
        r_in::T,
        r_ex::T,
        epsr_i::T,
        jω::Complex{T},
        loss_factor::T
) where {T <: Real}
    if isapprox(r_in, 0.0, atol = eps(T)) || isapprox(r_in, r_ex, atol = eps(T))
        # TODO: Implement consistent handling of admittance for bare conductors
        # Issue URL: https://github.com/Electa-Git/LineCableModels.jl/issues/17
        return zero(Complex{T})
    end

    eps0 = one(r_in) * 88541878128 * (one(r_in) * 10)^(-22)
    eps_i = eps0 * epsr_i

    return Complex{T}(log(r_ex / r_in) / (2 * (one(r_in) * π) * eps_i))
end

@inline function potential_coefficient(
        f::Lossless,
        ws,
        component_idx::Int,
        jω::Complex{T}
) where {T <: Real}
    return f(
        ws.r_ins_in[component_idx],
        ws.r_ins_ext[component_idx],
        ws.eps_ins[component_idx],
        jω,
        ws.tan_ins[component_idx]
    )
end
