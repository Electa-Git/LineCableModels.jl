defaults(::Type{Formula{:schelkunoff1934}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the inner, outer, and transfer surface impedances of a tubular conductor by the
formula `Formula{:schelkunoff1934}`, whose equations are in the section "Formulas" of
[`InternalImpedance`](@ref).

# Arguments

- `Formula{:schelkunoff1934}`: cylindrical surface-impedance formula.
- `r_in`: inner radius \\[m\\].
- `r_ex`: outer radius, in \\[m\\].
- `rho`: resistivity \\[Ω·m\\].
- `mu_r`: relative permeability \\[dimensionless\\].
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].

# Returns

- `inner`: impedance of the inner surface \\[Ω/m\\].
- `outer`: impedance of the outer surface, in \\[Ω/m\\].
- `transfer`: transfer impedance between the two surfaces \\[Ω/m\\].

# Errors

- Throws `DomainError` unless the radii are finite with `0 ≤ r_in < r_ex`, `rho` and `mu_r`
  are positive and finite, and `s` is finite and nonzero.
"""
function impedance(::Formula{:schelkunoff1934}, r_in::T, r_ex::T, rho::T, mu_r::T,
        s::Complex{T}) where {T <: Real}
    isfinite(r_in) && isfinite(r_ex) && zero(T) <= r_in < r_ex ||
        throw(DomainError((r_in, r_ex), "conductor radii must satisfy 0 ≤ r_in < r_ex [m]"))
    isfinite(rho) && rho > zero(T) && isfinite(mu_r) && mu_r > zero(T) ||
        throw(DomainError((rho, mu_r), "conductor resistivity and permeability must be positive and finite"))
    isfinite(s) && !iszero(s) || throw(DomainError(s, "jω must be finite and nonzero"))
    mu = vacuum_permeability * mu_r
    sigma = inv(rho)
    m = sqrt(s * mu * sigma)
    w_ex = m * r_ex
    if iszero(r_in)
        numerator = besselix(0, w_ex)
        denominator = besselix(1, w_ex)
        outer = (s * mu / 2π) * (1 / w_ex) * (numerator / denominator)
        return (inner = zero(outer), outer, transfer = zero(outer))
    end

    w_in = m * r_in
    sc_in = exp(abs(real(w_in)) - w_ex)
    sc_ex = exp(abs(real(w_ex)) - w_in)
    sc = sc_in / sc_ex
    i0_in = besselix(0, w_in)
    i1_in = besselix(1, w_in)
    k0_in = besselkx(0, w_in)
    k1_in = besselkx(1, w_in)
    i0_ex = besselix(0, w_ex)
    i1_ex = besselix(1, w_ex)
    k0_ex = besselkx(0, w_ex)
    k1_ex = besselkx(1, w_ex)
    numerator = k0_in * i1_ex + sc * i0_in * k1_ex
    denominator = k1_in * i1_ex - sc * i1_in * k1_ex
    inner = (s * mu / 2π) * (1 / w_in) * (numerator / denominator)
    numerator = i0_ex * k1_in + sc * k0_ex * i1_in
    denominator = i1_ex * k1_in - sc * k1_ex * i1_in
    outer = (s * mu / 2π) * (1 / w_ex) * (numerator / denominator)
    numerator = one(sc_ex) / sc_ex
    denominator = i1_ex * k1_in - sc * i1_in * k1_ex
    transfer = (1 / (2π * r_in * r_ex * sigma)) * (numerator / denominator)
    return (; inner, outer, transfer)
end
