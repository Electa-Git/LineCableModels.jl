"""
$(TYPEDEF)

Represent one conductor/insulation pair and its eagerly flattened equivalent
materials at their common reference temperature.

$(TYPEDFIELDS)
"""
struct CableComponent{T <: Real}
    id::String
    conductor_group::ConductorGroup{T}
    conductor_props::Material{T}
    insulator_group::InsulatorGroup{T}
    insulator_props::Material{T}
end

function _component_reference_temperature(component::CableComponent)
    return _reference_temperature((
        component.conductor_group.layers...,
        component.insulator_group.layers...
    ))
end

function _check_cable_component(component::CableComponent)
    isapprox(component.conductor_group.r_ex, component.insulator_group.r_in) ||
        throw(DomainError(
            (component.conductor_group.r_ex, component.insulator_group.r_in),
            "conductor and insulator boundaries must coincide"
        ))
    _component_reference_temperature(component)
    return nothing
end

function Validation.rules(::Type{<:CableComponent})
    (Validation.OwnerRule(:cable_component_boundaries, _check_cable_component),)
end

function CableComponent(
        id::AbstractString,
        conductor_group::ConductorGroup,
        insulator_group::InsulatorGroup
)
    T = promote_type(eltype(conductor_group), eltype(insulator_group))
    conductors = convert(ConductorGroup{T}, conductor_group)
    insulators = convert(InsulatorGroup{T}, insulator_group)
    isapprox(conductors.r_ex, insulators.r_in) || throw(DomainError(
        (conductors.r_ex, insulators.r_in),
        "conductor and insulator boundaries must coincide"
    ))
    reference_temperature = _reference_temperature((
        conductors.layers...,
        insulators.layers...
    ))

    conductor_rho = equivalent_rho(
        conductors.resistance,
        conductors.r_ex,
        conductors.r_in
    )
    conductor_mu = equivalent_mu(
        conductors.gmr,
        conductors.r_ex,
        conductors.r_in
    )
    conductor_props = Material(
        :conductor,
        conductor_rho,
        zero(T),
        conductor_mu,
        reference_temperature,
        conductors.alpha
    )

    insulator_eps = equivalent_eps(
        insulators.shunt_capacitance,
        insulators.r_ex,
        insulators.r_in
    )
    conductivity = equivalent_conductivity(
        insulators.shunt_conductance,
        insulators.r_in,
        insulators.r_ex
    )
    insulator_mu = equivalent_mu(insulators) * solenoid_factor(
        conductors.num_turns,
        conductors.r_ex,
        insulators.r_ex
    )
    insulator_props = Material(
        :insulator,
        inv(conductivity),
        insulator_eps,
        insulator_mu,
        reference_temperature,
        zero(T)
    )
    return validate(CableComponent{T}(
        String(id), conductors, conductor_props, insulators, insulator_props
    ))
end

function Tubular(component::CableComponent)
    group = component.conductor_group
    return Tubular(group.r_in, group.r_ex, component.conductor_props)
end

function Insulator(component::CableComponent)
    group = component.insulator_group
    return Insulator(group.r_in, group.r_ex, component.insulator_props)
end

function ConductorGroup(component::CableComponent{T}) where {T}
    original = component.conductor_group
    group = ConductorGroup(Tubular(component))
    group.num_wires = original.num_wires
    group.num_turns = original.num_turns
    return validate(group)
end

function InsulatorGroup(component::CableComponent)
    return InsulatorGroup(
        Insulator(component);
        reference_frequency = component.insulator_group.reference_frequency
    )
end

function Base.convert(::Type{CableComponent{T}}, component::CableComponent) where {T <:
                                                                                   Real}
    return CableComponent(
        component.id,
        convert(ConductorGroup{T}, component.conductor_group),
        convert(InsulatorGroup{T}, component.insulator_group)
    )
end

function Base.convert(::Type{CableComponent{T}}, component::CableComponent{T}) where {T <:
                                                                                      Real}
    component
end

include("base.jl")
