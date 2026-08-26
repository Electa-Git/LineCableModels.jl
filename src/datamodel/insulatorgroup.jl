"""
$(TYPEDEF)

Store concentric dielectric layers and their eager equivalent parallel-RC
value at `reference_frequency`.

$(TYPEDFIELDS)
"""
mutable struct InsulatorGroup{T <: Real} <: AbstractInsulatorPart{T}
    r_in::T
    r_ex::T
    cross_section::T
    shunt_capacitance::T
    shunt_conductance::T
    reference_frequency::T
    layers::Vector{AbstractInsulatorPart{T}}
end

_preview_layer_name(::InsulatorGroup) = "insulator group"

function preview_materials(group::InsulatorGroup)
    return Tuple(material
    for layer in group.layers
    for material in preview_materials(layer))
end

function preview_shapes(group::InsulatorGroup, context)
    shapes = PreviewPolygon[]
    for (index, layer) in enumerate(group.layers)
        nested = merge(
            context,
            (; include_label = index == 1 && context.include_label)
        )
        append!(shapes, preview_shapes(layer, nested))
    end
    return shapes
end

Base.eltype(::InsulatorGroup{T}) where {T} = T
Base.eltype(::Type{InsulatorGroup{T}}) where {T} = T

"""
$(TYPEDSIGNATURES)

Calculate the log-radius-weighted relative magnetic permeability of concentric
insulation layers:

```math
\\mu_{r,\\mathrm{eq}} =
\\frac{\\sum_k \\mu_{r,k}\\log(r_{k,\\mathrm{ex}}/r_{k,\\mathrm{in}})}
{\\log(r_\\mathrm{ex}/r_\\mathrm{in})}.
```
"""
function BaseParams.equivalent_mu(group::InsulatorGroup)
    group.r_in > zero(group.r_in) || throw(DomainError(
        group.r_in,
        "equivalent insulation permeability requires a positive inner radius"
    ))
    weighted = zero(group.r_in)
    for layer in group.layers
        weighted += layer.material_props.mu_r * log(layer.r_ex / layer.r_in)
    end
    return weighted / log(group.r_ex / group.r_in)
end

function _reference_temperature(layers::AbstractVector{<:AbstractInsulatorPart})
    isempty(layers) && throw(ArgumentError("an insulator group requires one layer"))
    reference = first(layers).material_props.T0
    all(layer -> isapprox(layer.material_props.T0, reference), layers) ||
        throw(ArgumentError("all insulator materials must share one reference temperature"))
    return reference
end

function _check_insulator_group(group::InsulatorGroup)
    isempty(group.layers) && throw(ArgumentError("an insulator group cannot be empty"))
    _reference_temperature(group.layers)
    group.reference_frequency > zero(group.reference_frequency) || throw(DomainError(
        group.reference_frequency,
        "reference frequency must be positive"
    ))
    group.r_in == first(group.layers).r_in ||
        throw(DomainError(group.r_in, "group inner radius differs from its first layer"))
    group.r_ex == last(group.layers).r_ex ||
        throw(DomainError(group.r_ex, "group outer radius differs from its last layer"))
    return nothing
end

function Validation.rules(::Type{<:InsulatorGroup})
    (Validation.OwnerRule(:insulator_group_structure, _check_insulator_group),)
end

function InsulatorGroup(
        layer::AbstractInsulatorPart{T};
        reference_frequency::Real = oftype(layer.r_in, 50)
) where {T <: Real}
    frequency = convert(T, reference_frequency)
    return validate(InsulatorGroup{T}(
        layer.r_in,
        layer.r_ex,
        layer.cross_section,
        layer.shunt_capacitance,
        layer.shunt_conductance,
        frequency,
        AbstractInsulatorPart{T}[layer]
    ))
end

Base.convert(::Type{InsulatorGroup{T}}, group::InsulatorGroup{T}) where {T <: Real} = group

"""
$(TYPEDSIGNATURES)

Append one concentric dielectric layer using the group's fixed reference
frequency. The layer must use the group's scalar type and reference temperature.
"""
function add!(group::InsulatorGroup{T}, layer::AbstractInsulatorPart{T}) where {T}
    isapprox(layer.r_in, group.r_ex) || throw(DomainError(
        layer.r_in,
        "new insulator layer must start at radius $(group.r_ex)"
    ))
    isapprox(layer.material_props.T0, _reference_temperature(group.layers)) ||
        throw(ArgumentError("new insulator material has a different reference temperature"))

    angular_frequency = (2 * (one(group.reference_frequency) * π)) *
                        group.reference_frequency
    group_admittance = complex(
        group.shunt_conductance,
        angular_frequency * group.shunt_capacitance
    )
    layer_admittance = complex(
        layer.shunt_conductance,
        angular_frequency * layer.shunt_capacitance
    )
    equivalent = parallel(group_admittance, layer_admittance)
    next_conductance = real(equivalent)
    next_capacitance = imag(equivalent) / angular_frequency
    next_cross_section = group.cross_section + layer.cross_section

    candidate = validate(InsulatorGroup{T}(
        group.r_in,
        layer.r_ex,
        next_cross_section,
        next_capacitance,
        next_conductance,
        group.reference_frequency,
        AbstractInsulatorPart{T}[group.layers; layer]
    ))
    group.r_ex = candidate.r_ex
    group.cross_section = candidate.cross_section
    group.shunt_capacitance = candidate.shunt_capacitance
    group.shunt_conductance = candidate.shunt_conductance
    group.layers = candidate.layers
    return group
end

function add!(group::InsulatorGroup{T}, layer::AbstractInsulatorPart{U}) where {T, U}
    throw(ArgumentError(
        "cannot add $(typeof(layer)) to InsulatorGroup{$T}; explicitly convert the " *
        "complete design before mutation",
    ))
end

function Base.convert(::Type{InsulatorGroup{T}}, group::InsulatorGroup) where {T <: Real}
    layers = AbstractInsulatorPart{T}[convert(AbstractInsulatorPart{T}, layer)
                                      for layer in group.layers]
    return validate(InsulatorGroup{T}(
        convert(T, group.r_in), convert(T, group.r_ex),
        convert(T, group.cross_section), convert(T, group.shunt_capacitance),
        convert(T, group.shunt_conductance), convert(T, group.reference_frequency),
        layers
    ))
end
