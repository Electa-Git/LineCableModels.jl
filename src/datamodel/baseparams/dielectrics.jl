@inline function _epsilon0(value)
    return one(value) * 88541878128 * (one(value) * 10)^(-22)
end

"""
$(TYPEDSIGNATURES)

Calculate coaxial shunt capacitance
``C=2\\pi\\varepsilon_0\\varepsilon_r/\\log(r_{ex}/r_{in})`` \\[F/m\\].
"""
function shunt_capacitance(r_in::Real, r_ex::Real, eps_r::Real)
    rin, rex, permittivity = promote(float(r_in), float(r_ex), float(eps_r))
    rin > zero(rin) && rex > rin || throw(DomainError(
        (rin, rex), "shunt capacitance requires 0 < r_in < r_ex"
    ))
    permittivity >= zero(permittivity) || throw(DomainError(
        permittivity, "relative permittivity must be nonnegative"
    ))
    return 2 * (one(rin) * π) * _epsilon0(rin) * permittivity / log(rex / rin)
end

"""
$(TYPEDSIGNATURES)

Calculate coaxial shunt conductance
``G=2\\pi/(\\rho\\log(r_{ex}/r_{in}))`` \\[S/m\\].
"""
function shunt_conductance(r_in::Real, r_ex::Real, rho::Real)
    rin, rex, resistivity = promote(float(r_in), float(r_ex), float(rho))
    rin > zero(rin) && rex > rin || throw(DomainError(
        (rin, rex), "shunt conductance requires 0 < r_in < r_ex"
    ))
    resistivity > zero(resistivity) || throw(DomainError(
        resistivity, "resistivity must be positive"
    ))
    return 2 * (one(rin) * π) / resistivity / log(rex / rin)
end

"""
$(TYPEDSIGNATURES)

Combine two radial dielectric layers connected in series. Each layer is
represented by the shunt admittance
``Y_i=G_i+\\mathrm{j}\\omega C_i``; the equivalent is

```math
Y_{eq}=\\frac{Y_1Y_2}{Y_1+Y_2}.
```

# Arguments

- `conductance1`: Shunt conductance of the accumulated dielectric \\[S/m\\].
- `capacitance1`: Shunt capacitance of the accumulated dielectric \\[F/m\\].
- `conductance2`: Shunt conductance of the added dielectric \\[S/m\\].
- `capacitance2`: Shunt capacitance of the added dielectric \\[F/m\\].
- `omega`: Angular reference frequency \\[rad/s\\].

# Returns

- Named tuple containing the equivalent `conductance` \\[S/m\\] and
  `capacitance` \\[F/m\\].
"""
function series_shunt_admittance(
        conductance1::Real,
        capacitance1::Real,
        conductance2::Real,
        capacitance2::Real,
        omega::Real
)
    G1, C1, G2, C2, angular_frequency = promote(
        float(conductance1), float(capacitance1),
        float(conductance2), float(capacitance2), float(omega)
    )
    angular_frequency > zero(angular_frequency) || throw(DomainError(
        angular_frequency, "angular frequency must be positive"
    ))
    equivalent = parallel(
        complex(G1, angular_frequency * C1),
        complex(G2, angular_frequency * C2)
    )
    return (
        conductance = real(equivalent),
        capacitance = imag(equivalent) / angular_frequency
    )
end

"""
$(TYPEDSIGNATURES)

Recover equivalent relative permittivity from coaxial capacitance:
``\\varepsilon_{eq}=C\\log(r_{ex}/r_{in})/(2\\pi\\varepsilon_0)``.
"""
function equivalent_eps(capacitance::Real, r_ex::Real, r_in::Real)
    value, rex, rin = promote(float(capacitance), float(r_ex), float(r_in))
    value >= zero(value) || throw(DomainError(value, "capacitance must be nonnegative"))
    rin > zero(rin) && rex > rin || throw(DomainError(
        (rin, rex), "equivalent permittivity requires 0 < r_in < r_ex"
    ))
    return value * log(rex / rin) / (2 * (one(rin) * π) * _epsilon0(rin))
end

"""
$(TYPEDSIGNATURES)

Calculate dielectric loss tangent ``\\tan\\delta=G/(\\omega C)``.
"""
function loss_tangent(conductance::Real, capacitance::Real, omega::Real)
    G, C, frequency = promote(
        float(conductance), float(capacitance), float(omega)
    )
    G >= zero(G) || throw(DomainError(G, "conductance must be nonnegative"))
    C > zero(C) || throw(DomainError(C, "capacitance must be positive"))
    frequency > zero(frequency) || throw(DomainError(
        frequency, "angular frequency must be positive"
    ))
    return G / (frequency * C)
end

"""
$(TYPEDSIGNATURES)

Recover equivalent conductivity from coaxial conductance:
``\\sigma_{eq}=G\\log(r_{ex}/r_{in})/(2\\pi)``.
"""
function equivalent_conductivity(conductance::Real, r_in::Real, r_ex::Real)
    value, rin, rex = promote(float(conductance), float(r_in), float(r_ex))
    value >= zero(value) || throw(DomainError(value, "conductance must be nonnegative"))
    rin > zero(rin) && rex > rin || throw(DomainError(
        (rin, rex), "equivalent conductivity requires 0 < r_in < r_ex"
    ))
    return value * log(rex / rin) / (2 * (one(rin) * π))
end
