"""
$(TYPEDEF)

Store a materialized cable design at one common material reference temperature.

$(TYPEDFIELDS)
"""
mutable struct CableDesign{T <: Real}
    cable_id::String
    nominal_data::Union{Nothing, NominalData{T}}
    components::Vector{CableComponent{T}}
end

function _design_reference_temperature(components)
    isempty(components) && throw(ArgumentError("a cable design requires one component"))
    reference = _component_reference_temperature(first(components))
    all(component -> isapprox(_component_reference_temperature(component), reference),
        components) || throw(ArgumentError(
        "all cable materials must share one reference temperature",
    ))
    return reference
end

function validate(design::CableDesign)
    isempty(design.cable_id) && throw(ArgumentError("cable_id cannot be empty"))
    _design_reference_temperature(design.components)
    ids = getproperty.(design.components, :id)
    allunique(ids) || throw(ArgumentError("component identifiers must be unique"))
    for (left, right) in zip(design.components, Iterators.drop(design.components, 1))
        isapprox(left.insulator_group.r_ex, right.conductor_group.r_in) ||
            throw(DomainError(
                (left.insulator_group.r_ex, right.conductor_group.r_in),
                "adjacent cable-component boundaries must coincide"
            ))
    end
    return design
end

function CableDesign(
        cable_id::AbstractString,
        component::CableComponent{T};
        nominal_data::Union{Nothing, NominalData} = nothing
) where {T <: Real}
    nominal = nominal_data === nothing ? nothing : convert(NominalData{T}, nominal_data)
    return validate(CableDesign{T}(
        String(cable_id), nominal, CableComponent{T}[component]
    ))
end

function CableDesign(
        cable_id::AbstractString,
        components::AbstractVector{<:CableComponent};
        nominal_data::Union{Nothing, NominalData} = nothing
)
    isempty(components) && throw(ArgumentError("a cable design requires one component"))
    T = promote_type(eltype.(components)...)
    converted = CableComponent{T}[convert(CableComponent{T}, item) for item in components]
    nominal = nominal_data === nothing ? nothing : convert(NominalData{T}, nominal_data)
    return validate(CableDesign{T}(String(cable_id), nominal, converted))
end

function CableDesign(
        cable_id::AbstractString,
        conductor_group::ConductorGroup,
        insulator_group::InsulatorGroup;
        component_id::AbstractString = "component1",
        nominal_data::Union{Nothing, NominalData} = nothing
)
    return CableDesign(
        cable_id,
        CableComponent(component_id, conductor_group, insulator_group);
        nominal_data
    )
end

function add!(design::CableDesign{T}, component::CableComponent{T}) where {T}
    any(item -> item.id == component.id, design.components) && throw(ArgumentError(
        "component '$(component.id)' already exists",
    ))
    reference = _design_reference_temperature(design.components)
    isapprox(_component_reference_temperature(component), reference) ||
        throw(ArgumentError("component materials have a different reference temperature"))
    isempty(design.components) ||
        isapprox(
            last(design.components).insulator_group.r_ex,
            component.conductor_group.r_in
        ) ||
        throw(DomainError(
            component.conductor_group.r_in,
            "component must start at the current cable outer radius"
        ))
    push!(design.components, component)
    return validate(design)
end

function add!(design::CableDesign{T}, component::CableComponent{U}) where {T, U}
    throw(ArgumentError(
        "cannot add CableComponent{$U} to CableDesign{$T}; explicitly convert the " *
        "complete design before mutation",
    ))
end

function add!(
        design::CableDesign,
        component_id::AbstractString,
        conductor_group::ConductorGroup,
        insulator_group::InsulatorGroup
)
    return add!(design, CableComponent(component_id, conductor_group, insulator_group))
end

function equivalent(original::CableDesign; new_id::AbstractString = "")
    target = isempty(new_id) ? original.cable_id * "_equivalent" : String(new_id)
    components = CableComponent[CableComponent(
                                    component.id,
                                    ConductorGroup(component),
                                    InsulatorGroup(component)
                                ) for component in original.components]
    return CableDesign(target, components; nominal_data = original.nominal_data)
end

function nonsensify(original::CableDesign; new_id::AbstractString = "")
    target = isempty(new_id) ? original.cable_id * "_nonsense" : String(new_id)
    components = CableComponent[]
    for component in original.components
        conductors = component.conductor_group
        insulators = component.insulator_group
        conductor = Tubular(
            first(conductors.layers).r_in,
            last(conductors.layers).r_ex,
            first(conductors.layers).material_props
        )
        dielectric_index = something(
            findfirst(layer -> layer isa Insulator, insulators.layers),
            1
        )
        dielectric = Insulator(
            conductor.r_ex,
            insulators.r_ex,
            insulators.layers[dielectric_index].material_props
        )
        push!(components,
            CableComponent(
                component.id,
                ConductorGroup(conductor),
                InsulatorGroup(
                    dielectric;
                    reference_frequency = insulators.reference_frequency
                )
            ))
    end
    return CableDesign(target, components; nominal_data = original.nominal_data)
end

function Base.convert(::Type{CableDesign{T}}, design::CableDesign) where {T <: Real}
    return CableDesign(
        design.cable_id,
        CableComponent{T}[convert(CableComponent{T}, item) for item in design.components];
        nominal_data = design.nominal_data === nothing ? nothing :
                       convert(NominalData{T}, design.nominal_data)
    )
end

Base.convert(::Type{CableDesign{T}}, design::CableDesign{T}) where {T <: Real} = design

include("cabledesign/base.jl")
include("cabledesign/cableconstants.jl")
include("cabledesign/dataframe.jl")
