struct Lossless <: InsulationAdmittanceFormulation end
get_description(::Lossless) = "Lossless insulation (ideal dielectric)"

@inline function (f::Lossless)(
        r_in::T,
        r_ex::T,
        epsr_i::T,
        jω::Complex{T},
        loss_factor::T
) where {T <: REALSCALAR}
    if isapprox(r_in, 0.0, atol = eps(T)) || isapprox(r_in, r_ex, atol = eps(T))
        # TODO: Implement consistent handling of admittance for bare conductors
        # Issue URL: https://github.com/Electa-Git/LineCableModels.jl/issues/17
        return zero(Complex{T})
    end

    # Constants
    eps_i = T(ε₀) * epsr_i

    return Complex{T}(log(r_ex / r_in) / (2π * eps_i))
end

@inline function potential_coefficient(
        f::Lossless,
        ws,
        component_idx::Int,
        jω::Complex{T}
) where {T <: REALSCALAR}
    return f(
        ws.r_ins_in[component_idx],
        ws.r_ins_ext[component_idx],
        ws.eps_ins[component_idx],
        jω,
        ws.tan_ins[component_idx]
    )
end
