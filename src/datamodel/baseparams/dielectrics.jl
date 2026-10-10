"""
$(TYPEDSIGNATURES)

Recover equivalent relative permittivity from coaxial capacitance:
``\\varepsilon_{eq}=C\\log(r_{ex}/r_{in})/(2\\pi\\varepsilon_0)``.

`capacitance` is nonnegative capacitance per unit length in F/m.
`r_ex` and `r_in` are the outer and inner dielectric radii in meters,
with ``0<r_{in}<r_{ex}``. ``ε_0`` is vacuum permittivity in F/m.
The returned relative permittivity is dimensionless.
"""
function equivalent_eps(capacitance::Real, r_ex::Real, r_in::Real)
    value, rex, rin = promote(float(capacitance), float(r_ex), float(r_in))
    value >= zero(value) || throw(DomainError(value, "capacitance must be nonnegative"))
    rin > zero(rin) && rex > rin || throw(DomainError(
        (rin, rex), "equivalent permittivity requires 0 < r_in < r_ex"
    ))
    return value * log(rex / rin) / (2 * (one(rin) * π) * vacuum_permittivity(typeof(rin)))
end

"""
$(TYPEDSIGNATURES)

Recover equivalent conductivity from coaxial conductance:
``\\sigma_{eq}=G\\log(r_{ex}/r_{in})/(2\\pi)``.

`conductance` is nonnegative conductance per unit length in S/m.
`r_in` and `r_ex` are the inner and outer dielectric radii in meters,
with ``0<r_{in}<r_{ex}``. The returned conductivity is in S/m.
"""
function equivalent_conductivity(conductance::Real, r_in::Real, r_ex::Real)
    value, rin, rex = promote(float(conductance), float(r_in), float(r_ex))
    value >= zero(value) || throw(DomainError(value, "conductance must be nonnegative"))
    rin > zero(rin) && rex > rin || throw(DomainError(
        (rin, rex), "equivalent conductivity requires 0 < r_in < r_ex"
    ))
    return value * log(rex / rin) / (2 * (one(rin) * π))
end
