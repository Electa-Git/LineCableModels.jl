defaults(::Type{Formula{:schelkunoff1934}}) = (;)

function (::Formula{:schelkunoff1934})(
        r_in::T, r_ex::T, ρ::T, μ_r::T, s::Complex{T}) where {T <: Real}
    isfinite(r_in) && isfinite(r_ex) && zero(T) <= r_in < r_ex ||
        throw(DomainError((r_in, r_ex), "conductor radii must satisfy 0 ≤ r_in < r_ex [m]"))
    isfinite(ρ) && ρ > zero(T) && isfinite(μ_r) && μ_r > zero(T) ||
        throw(DomainError((ρ, μ_r), "conductor resistivity and permeability must be positive and finite"))
    isfinite(s) && !iszero(s) || throw(DomainError(s, "jω must be finite and nonzero"))
    mu = vacuum_permeability * μ_r
    sigma = inv(ρ)
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
