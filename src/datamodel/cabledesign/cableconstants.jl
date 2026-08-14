"""
    CableConstants{T}

Canonical per-unit-length cable constants.

The fields `R`, `L`, and `C` are stored in Ω/m, H/m, and F/m respectively.
Display conversions belong to `UnitHandler` and presentation adapters.
"""
struct CableConstants{T}
    "Series resistance per unit length \\[Ω/m\\]."
    R::T
    "Series inductance per unit length \\[H/m\\]."
    L::T
    "Shunt capacitance per unit length \\[F/m\\]."
    C::T
end

function Base.:(==)(left::CableConstants, right::CableConstants)
    left.R == right.R && left.L == right.L && left.C == right.C
end

function CableConstants(R::Real, L::Real, C::Real)
    values = promote(R, L, C)
    return CableConstants{typeof(first(values))}(values...)
end

"""
    CableConstants(design::CableDesign; S=nothing, rho_e=100.0)

Compute the scalar cable constants represented by `design`.

# Arguments

- `design`: Cable design whose core, shield, and outer geometry are used.
- `S`: Cable separation in metres. The outer diameter is used when omitted.
- `rho_e`: Earth resistivity in Ω·m.

# Returns

A [`CableConstants`](@ref) value storing `R` in Ω/m, `L` in H/m, and `C` in
F/m.

# Notes

The implemented expressions are the same tubular-resistance, trefoil-
inductance, and coaxial-capacitance expressions historically used by
`DataFrame(design, :baseparams)`. This constructor changes their storage units,
not their physical calculation.
"""
function CableConstants(
        design::CableDesign;
        S::Union{Nothing, Number} = nothing,
        rho_e::Number = 100.0
)
    length(design.components) >= 2 || throw(
        ArgumentError("at least two cable components are required"),
    )

    cable_core = design.components[1]
    cable_shield = design.components[2]
    cable_outer = design.components[end]
    separation = if S === nothing
        if isnan(cable_outer.insulator_group.r_ex)
            2 * cable_outer.conductor_group.r_ex
        else
            2 * cable_outer.insulator_group.r_ex
        end
    else
        S
    end

    resistance_value = calc_tubular_resistance(
        cable_core.conductor_group.r_in,
        cable_core.conductor_group.r_ex,
        cable_core.conductor_props.rho,
        0.0,
        20.0,
        20.0
    )
    inductance_value = calc_inductance_trifoil(
        cable_core.conductor_group.r_in,
        cable_core.conductor_group.r_ex,
        cable_core.conductor_props.rho,
        cable_core.conductor_props.mu_r,
        cable_shield.conductor_group.r_in,
        cable_shield.conductor_group.r_ex,
        cable_shield.conductor_props.rho,
        cable_shield.conductor_props.mu_r,
        separation;
        rho_e = rho_e
    )
    capacitance_value = calc_shunt_capacitance(
        cable_core.conductor_group.r_ex,
        cable_core.insulator_group.r_ex,
        cable_core.insulator_props.eps_r
    )
    return CableConstants(resistance_value, inductance_value, capacitance_value)
end

R(constants::CableConstants) = constants.R
L(constants::CableConstants) = constants.L
C(constants::CableConstants) = constants.C
basis(::CableConstants) = :per_length
resistance(constants::CableConstants) = R(constants)
inductance(constants::CableConstants) = L(constants)
capacitance(constants::CableConstants) = C(constants)
