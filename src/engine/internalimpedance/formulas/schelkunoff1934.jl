
"""
$(TYPEDSIGNATURES)

**Identification.** Exact cylindrical surface impedances for a solid or
hollow round conductor.

**Expression.**

```math
\\begin{aligned}
Z_{is}&=\\frac{\\rho m}{2\\pi aD}
[I_0(ma)K_1(mb)+K_0(ma)I_1(mb)],\\\\
Z_{os}&=\\frac{\\rho m}{2\\pi bD}
[I_0(mb)K_1(ma)+K_0(mb)I_1(ma)],\\\\
Z_{ms}&=\\frac{\\rho m}{2\\pi abD},\\\\
D&=I_1(mb)K_1(ma)-K_1(mb)I_1(ma).
\\end{aligned}
```

Here ``a`` and ``b`` are the inner and outer conductor radii in meters,
``ρ`` is resistivity in Ω·m, and ``m=\\sqrt{jωμ/ρ}`` is in m⁻¹, with
``μ=μ_0μ_r`` in H/m. ``I_ν`` and ``K_ν`` are modified Bessel functions.
Each surface impedance is in Ω/m.

For ``a=0``, ``Z_{int}=\\rho mI_0(mb)/(2\\pi bI_1(mb))``.

Schelkunoff's surface terms were later recovered by Ametani to assemble the
complete core-sheath-armor impedance matrix. The Engine applies that outward
assembly recursively to any number of concentric conductive terminals.

**Reference.** S. A. Schelkunoff, “The Electromagnetic Theory of Coaxial
Transmission Lines and Cylindrical Shields,” *Bell System Technical Journal*,
13, 532–579, 1934. A. Ametani, “A General Formulation of Impedance and
Admittance of Cables,” *IEEE Transactions on Power Apparatus and Systems*,
PAS-99(3), 902–910, 1980. DOI: 10.1109/TPAS.1980.319718.

"""
function description(::Type{<:Formula{:schelkunoff1934}}; compact::Bool=false)
    compact ? "Schelkunoff" : "Schelkunoff exact round-conductor surface impedances (1934)"
end

"""
$(TYPEDSIGNATURES)

Build the Functor of one solid or hollow circular conductor at one frequency for the exact
Schelkunoff surface impedances:

```math
Z_{is}=\\frac{\\rho m}{2\\pi aD}
\\left[I_0(ma)K_1(mb)+K_0(ma)I_1(mb)\\right],
\\qquad
Z_{ms}=\\frac{\\rho m}{2\\pi abD},
```

```math
Z_{os}=\\frac{\\rho m}{2\\pi bD}
\\left[I_0(mb)K_1(ma)+K_0(mb)I_1(ma)\\right],
```

where

```math
D=I_1(mb)K_1(ma)-K_1(mb)I_1(ma),
\\qquad
m=\\sqrt{j\\omega\\mu/\\rho}.
```

Here ``μ=μ_0μ_r`` is absolute permeability in H/m, and ``I_ν`` and ``K_ν``
are modified Bessel functions. For ``a=0``, the outer term is evaluated from the solid-cylinder limit
``Z_{int}=\\rho m I_0(mb)/(2\\pi b I_1(mb))``.

# Input

- `r_in`: inner conductor radius ``a`` \\[m\\].
- `r_ex`: outer conductor radius ``b`` \\[m\\].
- `rho`: conductor resistivity ``\\rho`` \\[Ω·m\\].
- `mu_r`: relative conductor permeability \\[dimensionless\\].
- `jω`: complex angular frequency ``j\\omega`` \\[rad/s\\].

# Returns

- The Functor of the conductor at that frequency. Its state stores the scaled modified
  Bessel functions ``I_0``, ``I_1``, ``K_0`` and ``K_1`` at ``ma`` and ``mb``, evaluated once
  for the inner, outer and transfer surfaces. A solid conductor evaluates only ``I_0(mb)``
  and ``I_1(mb)``.

# Notes

Implements Schelkunoff's cylindrical surface terms. Ametani (1980) recovered
these terms for the complete core-sheath-armor impedance assembly performed
recursively by the Engine.
"""
function Functor(formula::Formula{:schelkunoff1934}, input::NamedTuple; workspace = nothing)
    (; r_in, r_ex, rho, mu_r, jω) = input
    T = typeof(r_in)
    isfinite(r_in) && isfinite(r_ex) && zero(T) <= r_in < r_ex ||
        throw(DomainError((r_in, r_ex), "conductor radii must satisfy 0 ≤ r_in < r_ex [m]"))
    isfinite(rho) && rho > zero(T) && isfinite(mu_r) && mu_r > zero(T) ||
        throw(DomainError((rho, mu_r), "conductor resistivity and permeability must be positive and finite"))
    isfinite(jω) && !iszero(jω) || throw(DomainError(jω, "jω must be finite and nonzero"))
    mu_c = vacuum_permeability(T) * mu_r
    sigma_c = conductivity(rho)
    m = sqrt(jω * mu_c * sigma_c)
    w_ex = m * r_ex
    w_in = m * r_in
    i0_ex = special_besselix(0, w_ex)
    i1_ex = special_besselix(1, w_ex)
    # A solid conductor has no inner surface, and K at zero would be infinite. Its values
    # that the outer surface does not read are zero, so the state type does not depend on
    # the geometry.
    absent = zero(i0_ex)
    sc_ex, sc, i0_in, i1_in, k0_in, k1_in, k0_ex, k1_ex =
        if isapprox(r_in, zero(T); atol = eps(T))
            ntuple(_ -> absent, 8)
        else
            sc_in = exp(abs(real(w_in)) - w_ex)
            sc_out = exp(abs(real(w_ex)) - w_in)
            (sc_out, sc_in / sc_out, special_besselix(0, w_in), special_besselix(1, w_in),
                special_besselkx(0, w_in), special_besselkx(1, w_in),
                special_besselkx(0, w_ex), special_besselkx(1, w_ex))
        end
    state = (; mu_c, sigma_c, m, w_in, w_ex, sc_ex, sc,
        i0_in, i1_in, k0_in, k1_in, i0_ex, i1_ex, k0_ex, k1_ex)
    return Functor(formula, input, state)
end

@inline function internal_impedance(
        ::Formula{:schelkunoff1934},
        ::Val{:inner},
        functor, workspace
)
    input, state = functor.input, functor.state
    T = typeof(input.r_in)
    if isapprox(input.r_in, zero(T); atol = eps(T))
        return zero(Complex{T})
    end

    numerator = state.k0_in * state.i1_ex + state.sc * state.i0_in * state.k1_ex
    denominator = state.k1_in * state.i1_ex - state.sc * state.i1_in * state.k1_ex
    return Complex{T}(
        (input.jω * state.mu_c / 2π) * (1 / state.w_in) * (numerator / denominator)
    )
end

@inline function internal_impedance(
        ::Formula{:schelkunoff1934},
        ::Val{:outer},
        functor, workspace
)
    input, state = functor.input, functor.state
    T = typeof(input.r_in)
    if isapprox(input.r_in, zero(T); atol = eps(T))
        numerator = state.i0_ex
        denominator = state.i1_ex
    else
        numerator = state.i0_ex * state.k1_in + state.sc * state.k0_ex * state.i1_in
        denominator = state.i1_ex * state.k1_in - state.sc * state.k1_ex * state.i1_in
    end
    return Complex{T}(
        (input.jω * state.mu_c / 2π) * (1 / state.w_ex) *
        (numerator / denominator)
    )
end

@inline function internal_impedance(
        ::Formula{:schelkunoff1934},
        ::Val{:transfer},
        functor, workspace
)
    input, state = functor.input, functor.state
    T = typeof(input.r_in)
    if isapprox(input.r_in, zero(T); atol = eps(T))
        return zero(Complex{T})
    end

    numerator = one(state.sc_ex) / state.sc_ex
    denominator = state.i1_ex * state.k1_in - state.sc * state.i1_in * state.k1_ex
    return Complex{T}(
        (1 / (2π * input.r_in * input.r_ex * state.sigma_c)) *
        (numerator / denominator)
    )
end

formulation_options(::Expression{<:Formula{:schelkunoff1934}, typeof(internal_impedance)}) = FormulationOptions()

:schelkunoff1934
