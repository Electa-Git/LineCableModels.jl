defaults(::Type{Formula{:schelkunoff1934}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the surface impedances of a tubular conductor by Schelkunoff's cylindrical surface
terms, with the modified Bessel functions ``I_\\nu`` and ``K_\\nu``:

```math
\\begin{aligned}
z_\\mathrm{inner} &= \\frac{s \\mu}{2 \\pi w_\\mathrm{in}}
    \\frac{I_0(w_\\mathrm{in}) K_1(w_\\mathrm{ex}) + K_0(w_\\mathrm{in}) I_1(w_\\mathrm{ex})}{D}, \\\\
z_\\mathrm{outer} &= \\frac{s \\mu}{2 \\pi w_\\mathrm{ex}}
    \\frac{I_0(w_\\mathrm{ex}) K_1(w_\\mathrm{in}) + K_0(w_\\mathrm{ex}) I_1(w_\\mathrm{in})}{D}, \\\\
z_\\mathrm{transfer} &= \\frac{1}{2 \\pi r_\\mathrm{in} r_\\mathrm{ex} \\sigma D},
\\end{aligned}
```

where ``D = I_1(w_\\mathrm{ex}) K_1(w_\\mathrm{in}) - I_1(w_\\mathrm{in}) K_1(w_\\mathrm{ex})``,
``w = m r``, ``m = \\sqrt{s \\mu \\sigma}``, ``\\mu = \\mu_0 \\mu_r`` and ``\\sigma = 1 / \\rho``.
The computation uses exponentially scaled Bessel functions. A solid conductor, with
`r_in = 0`, has no inner surface: its inner and transfer impedances are zero, and its outer
impedance is ``s \\mu I_0(w_\\mathrm{ex}) / (2 \\pi w_\\mathrm{ex} I_1(w_\\mathrm{ex}))``.

# Arguments

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
