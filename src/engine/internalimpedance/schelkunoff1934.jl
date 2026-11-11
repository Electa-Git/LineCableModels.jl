defaults(::Type{Formula{:schelkunoff1934}}) = (;)

function (::Formula{:schelkunoff1934})(
        r_in::T, r_ex::T, ρ::T, μ_r::T, s::Complex{T}) where {T <: Real}
    isfinite(r_in) && isfinite(r_ex) && zero(T) <= r_in < r_ex ||
        throw(DomainError((r_in, r_ex), "conductor radii must satisfy 0 ≤ r_in < r_ex [m]"))
    isfinite(ρ) && ρ > zero(T) && isfinite(μ_r) && μ_r > zero(T) ||
        throw(DomainError((ρ, μ_r), "conductor resistivity and permeability must be positive and finite"))
    isfinite(s) && !iszero(s) || throw(DomainError(s, "jω must be finite and nonzero"))
    mu = μ₀ * μ_r
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

function (f::Formula{:schelkunoff1934})(
        n::Int, r_in::T, r_ex::T, ρ::T, μ_r::T, s::Complex{T}) where {T <: Real}
    n >= 0 || throw(DomainError(n, "the harmonic order must be nonnegative"))
    iszero(n) && return f(r_in, r_ex, ρ, μ_r, s)
    isfinite(r_in) && isfinite(r_ex) && zero(T) <= r_in < r_ex ||
        throw(DomainError((r_in, r_ex), "conductor radii must satisfy 0 ≤ r_in < r_ex [m]"))
    isfinite(ρ) && ρ > zero(T) && isfinite(μ_r) && μ_r > zero(T) ||
        throw(DomainError((ρ, μ_r), "conductor resistivity and permeability must be positive and finite"))
    isfinite(s) && !iszero(s) || throw(DomainError(s, "jω must be finite and nonzero"))
    μ = μ₀ * μ_r
    σ = inv(ρ)
    m = sqrt(s * μ * σ)
    w_in = m * r_in
    w_ex = m * r_ex
    N = n + ceil(Int, nominal(abs(w_ex))) + 80
    α_in = α_ex = complex(T(2N))
    for k in (N - 1):-1:n
        α_in = 2k + w_in^2 / α_in
        α_ex = 2k + w_ex^2 / α_ex
    end
    if iszero(r_in)
        outer = (s * μ / 2π) / (α_ex - n)
        return (inner = zero(outer), outer, transfer = zero(outer))
    end
    λ = besselix(0, w_in) * besselkx(0, w_ex) / (besselix(0, w_ex) * besselkx(0, w_in)) *
        exp(abs(real(w_in)) - abs(real(w_ex)) + w_in - w_ex)
    Π = besselix(0, w_ex) * besselkx(0, w_in)
    β_in = w_in * besselkx(0, w_in) / besselkx(1, w_in)
    β_ex = w_ex * besselkx(0, w_ex) / besselkx(1, w_ex)
    λ *= α_ex * β_in / (α_in * β_ex)
    Π *= w_ex * w_in / (α_ex * β_in)
    αk_in, αk_ex = α_in, α_ex
    for k in (n - 1):-1:1
        αk_in = 2k + w_in^2 / αk_in
        αk_ex = 2k + w_ex^2 / αk_ex
        β_in = w_in^2 / (β_in + 2(n - k))
        β_ex = w_ex^2 / (β_ex + 2(n - k))
        λ *= αk_ex * β_in / (αk_in * β_ex)
        Π *= w_ex * w_in / (αk_ex * β_in)
    end
    den = (β_in + n) * (α_ex - n) - λ * (α_in - n) * (β_ex + n)
    inner = (s * μ / 2π) * ((α_ex - n) + λ * (β_ex + n)) / den
    outer = (s * μ / 2π) * ((β_in + n) + λ * (α_in - n)) / den
    transfer = (s * μ / 2π) * exp(w_in - abs(real(w_ex))) / (Π * den)
    return (; inner, outer, transfer)
end
