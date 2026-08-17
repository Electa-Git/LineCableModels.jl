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
