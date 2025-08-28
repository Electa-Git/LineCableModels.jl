import ..LineCableModels: add!

"""
$(TYPEDEF)

Represents the design of a cable, including its unique identifier, nominal data, and components.

$(TYPEDFIELDS)
"""
mutable struct CableDesign{T<:REALSCALAR}
    "Unique identifier for the cable design."
    cable_id::String
    "Informative reference data."
    nominal_data::Union{Nothing,NominalData{T}}
    "Vector of cable components."
    components::Vector{CableComponent{T}}

    @doc """
    $(TYPEDSIGNATURES)

    **Strict numeric kernel**: constructs a `CableDesign{T}` from one component
    (typed) and optional nominal data (typed or `nothing`). Assumes all inputs
    are already at scalar type `T`.

    # Arguments

    - `cable_id`: Unique identifier for the cable design.
    - `component`: Initial [`CableComponent`](@ref) for the design.
    - `nominal_data`: Reference data for the cable design. Default: `NominalData()`.

    # Returns

    - A [`CableDesign`](@ref) object with the specified properties.

    # Examples

    ```julia
    conductor_group = ConductorGroup(central_conductor)
    insulator_group = InsulatorGroup(main_insulator)
    component = CableComponent(conductor_group, insulator_group)
    design = $(FUNCTIONNAME)("example", component)
    ```

    # See also

    - [`CableComponent`](@ref)
    - [`ConductorGroup`](@ref)
    - [`InsulatorGroup`](@ref)
    """
    @inline function CableDesign{T}(
        cable_id::String,
        component::CableComponent{T};
        nominal_data::Union{Nothing,NominalData{T}}=nothing,
    ) where {T<:REALSCALAR}
        new{T}(cable_id, nominal_data, CableComponent{T}[component])
    end

    @inline function CableDesign{T}(
        cable_id::String,
        components::Vector{CableComponent{T}};
        nominal_data::Union{Nothing,NominalData{T}}=nothing,
    ) where {T<:REALSCALAR}
        new{T}(cable_id, nominal_data, components)
    end
end

"""
$(TYPEDSIGNATURES)

**Weakly-typed constructor** that infers the scalar type from the `component` (and nominal data if present), coerces values to that type, and calls the typed kernel.
"""
function CableDesign(
    cable_id::String,
    component::CableComponent;
    nominal_data::NominalData=NominalData(),
)
    # Resolve T from component and nominal_data (ignoring `nothing` fields in the latter)
    T = resolve_T(component, nominal_data)

    compT = coerce_to_T(component, T)
    ndT = coerce_to_T(nominal_data, T)  # identity if already T

    return CableDesign{T}(cable_id, compT; nominal_data=ndT)
end

"""
$(TYPEDSIGNATURES)

Constructs a [`CableDesign`](@ref) instance **from conductor and insulator groups**.
Convenience wrapper that builds the component with reduced boilerplate.
"""
function CableDesign(
    cable_id::String,
    conductor_group::ConductorGroup,
    insulator_group::InsulatorGroup;
    component_id::String="component1",
    nominal_data::NominalData=NominalData(),
)
    component = CableComponent(component_id, conductor_group, insulator_group)
    return CableDesign(cable_id, component; nominal_data)
end

function add!(design::CableDesign{T}, component::CableComponent) where {T}
    Tnew = resolve_T(design, component)

    if Tnew === T
        compT = coerce_to_T(component, T)
        if (idx = findfirst(c -> c.id == compT.id, design.components)) !== nothing
            @warn "Component with ID '$(compT.id)' already exists and will be overwritten."
            design.components[idx] = compT
        else
            push!(design.components, compT)
        end
        return design
    else
        @warn """
        Adding a `$Tnew` component to a `CableDesign{$T}` returns a **promoted** design.
        Capture the result:  design = add!(design, component)
        """
        # promote whole design, then insert coerced component
        promoted = coerce_to_T(design, Tnew)
        compT = coerce_to_T(component, Tnew)
        if (idx = findfirst(c -> c.id == compT.id, promoted.components)) !== nothing
            promoted.components[idx] = compT
        else
            push!(promoted.components, compT)
        end
        return promoted
    end
end

# --- add!(design, by groups): wraps the above ---
function add!(
    design::CableDesign{T},
    component_id::String,
    conductor_group::ConductorGroup,
    insulator_group::InsulatorGroup,
) where {T}
    comp = CableComponent(component_id, conductor_group, insulator_group)
    add!(design, comp)  # may return the same or a promoted design
end

include("cabledesign/base.jl")
include("cabledesign/dataframe.jl")
