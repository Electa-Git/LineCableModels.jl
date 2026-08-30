assumptions(::Val{:ParallelRC}) = (;)

description(::Formula{:ParallelRC}) =
    "Parallel-RC insulation (constant conductivity and permittivity)"

"""
$(TYPEDSIGNATURES)

Calculate the potential coefficient of one concentric lossy dielectric layer:

```math
p(s) = \\frac{s}{G+sC}
     = \\frac{\\ln(r_o/r_i)}{2\\pi}
       \\frac{s}{1/\\rho+s\\varepsilon_0\\varepsilon_r}.
```

# Arguments

- `r_in`: Inner layer radius \\[m\\].
- `r_ex`: Outer layer radius \\[m\\].
- `rho`: Dielectric resistivity \\[Ω·m\\].
- `eps_r`: Relative dielectric permittivity \\[dimensionless\\].
- `s`: Complex angular frequency \\[rad/s\\].
- `values`: Formula assumptions.

# Returns

- Complex potential coefficient per unit length \\[m/F\\].
"""
@inline function parallelrc(
        r_in::T,
        r_ex::T,
        rho::T,
        eps_r::T,
        s::Complex{T},
        values::NamedTuple
) where {T <: Real}
    if isapprox(r_in, zero(T); atol = eps(T)) ||
       isapprox(r_in, r_ex; atol = eps(T))
        return zero(Complex{T})
    end

    log_ratio = log(r_ex / r_in)
    eps0 = one(r_in) * 88541878128 * (one(r_in) * 10)^(-22)
    circumference = 2 * (one(r_in) * π)
    capacitance = circumference * eps0 * eps_r / log_ratio
    conductance = circumference * conductivity(rho) / log_ratio
    return s / (conductance + s * capacitance)
end

@inline function (::Functor{:ParallelRC})(
        r_in::T,
        r_ex::T,
        rho::T,
        eps_r::T,
        s::Complex{T},
        values::NamedTuple
) where {T <: Real}
    return parallelrc(r_in, r_ex, rho, eps_r, s, values)
end

:ParallelRC
