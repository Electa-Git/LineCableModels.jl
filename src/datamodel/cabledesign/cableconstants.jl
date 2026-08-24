"""
$(TYPEDEF)

Store cable constants per unit length.

The fields `R`, `L`, and `C` are stored in Ω/m, H/m, and F/m respectively.
Display conversions belong to [`LineCableModels.Units`](@ref) and presentation
adapters.

$(TYPEDFIELDS)
"""
struct CableConstants{T} <: AbstractProblemResult
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
$(TYPEDSIGNATURES)

Compute the scalar cable constants represented by `design`.

The calculation uses tubular resistance, trefoil inductance, and coaxial
capacitance. Presentation is handled separately after the result exists.

# Arguments

- `design`: Cable design whose core, shield, and outer geometry are used.
- `S`: Cable separation in metres. The outer diameter is used when omitted.
- `rho_e`: Earth resistivity in Ω·m.

# Returns

A [`CableConstants`](@ref) value storing `R` in Ω/m, `L` in H/m, and `C` in
F/m.
"""
function _compute_cable_constants(
        design::CableDesign{T};
        S::Union{Nothing, Number} = nothing,
        rho_e::Number = oftype(first(design.components).conductor_props.rho, 100)
) where {T}
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

    reference = cable_core.conductor_props.T0
    resistance_value = tubular_resistance(
        cable_core.conductor_group.r_in,
        cable_core.conductor_group.r_ex,
        cable_core.conductor_props.rho,
        cable_core.conductor_props.alpha,
        reference,
        reference
    )
    inductance_value = trifoil_inductance(
        cable_core.conductor_group.r_in,
        cable_core.conductor_group.r_ex,
        cable_core.conductor_props.rho,
        cable_core.conductor_props.mu_r,
        cable_shield.conductor_group.r_in,
        cable_shield.conductor_group.r_ex,
        cable_shield.conductor_props.rho,
        cable_shield.conductor_props.mu_r,
        separation;
        rho_e,
        f = oftype(separation, 50)
    )
    capacitance_value = shunt_capacitance(
        cable_core.conductor_group.r_ex,
        cable_core.insulator_group.r_ex,
        cable_core.insulator_props.eps_r
    )
    return CableConstants(resistance_value, inductance_value, capacitance_value)
end

observe(constants::CableConstants, ::typeof(R)) = constants.R
observe(constants::CableConstants, ::typeof(L)) = constants.L
observe(constants::CableConstants, ::typeof(C)) = constants.C

R(constants::CableConstants) = observe(constants, R)
L(constants::CableConstants) = observe(constants, L)
C(constants::CableConstants) = observe(constants, C)
basis(::CableConstants) = :per_length
resistance(constants::CableConstants) = R(constants)
inductance(constants::CableConstants) = L(constants)
capacitance(constants::CableConstants) = C(constants)

observables(::Type{<:CableConstants}) = (R, L, C)

function Base.show(io::IO, constants::CableConstants)
    selectors = (R, L, C)
    print(io, "CableConstants(")
    for (index, selector) in enumerate(selectors)
        index > 1 && print(io, ", ")
        quantity = Units.quantity(selector)
        unit = Units.native_unit(quantity, basis(constants))
        print(
            io,
            Units.symbol(quantity),
            "=",
            observe(constants, selector),
            " ",
            Units.label(unit)
        )
    end
    print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", constants::CableConstants) = show(io, constants)
