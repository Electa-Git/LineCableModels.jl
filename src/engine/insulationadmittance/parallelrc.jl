"""
$(TYPEDEF)

Represent each concentric insulation layer as a frequency-independent shunt
conductance in parallel with its capacitance.

For a layer with inner radius ``r_i`` \\[m\\], outer radius ``r_o`` \\[m\\],
resistivity ``\\rho`` \\[Ω·m\\], and relative permittivity
``\\varepsilon_r`` \\[dimensionless\\], the per-unit-length layer admittance is

```math
y(s) = G + sC,
\\qquad
C = \\frac{2\\pi\\varepsilon_0\\varepsilon_r}{\\ln(r_o/r_i)},
\\qquad
G = \\frac{2\\pi}{\\rho\\ln(r_o/r_i)}.
```

The EMT solver assembles potential coefficients, so this formulation returns
``p(s)=s/y(s)`` for every layer and adds series layer coefficients before matrix
assembly. Material uncertainties therefore propagate through the complete
admittance calculation.
"""
struct ParallelRC <: InsulationAdmittanceFormulation end

function get_description(::ParallelRC)
    "Parallel-RC insulation (constant conductivity and permittivity)"
end

"""
$(TYPEDSIGNATURES)

Calculate the potential coefficient of one concentric lossy dielectric layer:

```math
p(s) = \\frac{s}{G+sC}
     = \\frac{\\ln(r_o/r_i)}{2\\pi}
       \\frac{s}{1/\\rho+s\\varepsilon_0\\varepsilon_r}.
```

At infinite resistivity, the result reduces exactly to the lossless potential
coefficient ``1/C``.

# Arguments

- `r_in`: Inner layer radius \\[m\\].
- `r_ex`: Outer layer radius \\[m\\].
- `rho`: Dielectric resistivity \\[Ω·m\\].
- `eps_r`: Relative dielectric permittivity \\[dimensionless\\].
- `s`: Complex angular frequency \\[rad/s\\].

# Returns

- Complex potential coefficient per unit length \\[m/F\\].
"""
@inline function (f::ParallelRC)(
        r_in::T,
        r_ex::T,
        rho::T,
        eps_r::T,
        s::Complex{T}
) where {T <: REALSCALAR}
    if isapprox(r_in, 0.0, atol = eps(T)) || isapprox(r_in, r_ex, atol = eps(T))
        # Keep bare-conductor handling consistent with Lossless.
        return zero(Complex{T})
    end

    log_ratio = log(r_ex / r_in)
    capacitance = T(2π * ε₀) * eps_r / log_ratio
    conductivity = _to_σ(rho)
    conductance = T(2π) * conductivity / log_ratio

    return s / (conductance + s * capacitance)
end

@inline function potential_coefficient(
        f::ParallelRC,
        ws,
        component_idx::Int,
        s::Complex{T}
) where {T <: REALSCALAR}
    p = zero(Complex{T})
    @inbounds for layer_idx in ws.insulator_layer_ranges[component_idx]
        p += f(
            ws.r_ins_layer_in[layer_idx],
            ws.r_ins_layer_ext[layer_idx],
            ws.rho_ins_layer[layer_idx],
            ws.eps_ins_layer[layer_idx],
            s
        )
    end
    return p
end
