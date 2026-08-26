"""
$(TYPEDEF)

Store concentric conductive layers and their eagerly calculated equivalent
properties at one material reference temperature.

$(TYPEDFIELDS)
"""
mutable struct ConductorGroup{T <: Real} <: AbstractConductorPart{T}
    r_in::T
    r_ex::T
    cross_section::T
    num_wires::Int
    num_turns::T
    resistance::T
    alpha::T
    gmr::T
    layers::Vector{AbstractConductorPart{T}}
end

_preview_layer_name(::ConductorGroup) = "conductor group"

function preview_materials(group::ConductorGroup)
    return Tuple(material
    for layer in group.layers
    for material in preview_materials(layer))
end

function preview_shapes(group::ConductorGroup, context)
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

Base.eltype(::ConductorGroup{T}) where {T} = T
Base.eltype(::Type{ConductorGroup{T}}) where {T} = T

BaseParams.gmd_elements(group::ConductorGroup) = BaseParams.gmd_elements(last(group.layers))

_wire_count(::AbstractConductorPart) = 0
_wire_count(layer::CircStrands) = layer.num_wires
_wire_count(::Strip) = 1

_turn_count(layer::AbstractConductorPart{T}) where {T} = zero(T)
function _turn_count(layer::CircStrands)
    iszero(layer.pitch_length) ?
    zero(eltype(layer)) : inv(layer.pitch_length)
end
function _turn_count(layer::Strip)
    iszero(layer.pitch_length) ?
    zero(eltype(layer)) : inv(layer.pitch_length)
end

function _reference_temperature(layers)
    isempty(layers) && throw(ArgumentError("a conductor group requires one layer"))
    reference = first(layers).material_props.T0
    all(layer -> isapprox(layer.material_props.T0, reference), layers) ||
        throw(ArgumentError("all conductor materials must share one reference temperature"))
    return reference
end

function _check_conductor_group(group::ConductorGroup)
    isempty(group.layers) && throw(ArgumentError("a conductor group cannot be empty"))
    _reference_temperature(group.layers)
    group.r_in == first(group.layers).r_in ||
        throw(DomainError(group.r_in, "group inner radius differs from its first layer"))
    group.r_ex == last(group.layers).r_ex ||
        throw(DomainError(group.r_ex, "group outer radius differs from its last layer"))
    return nothing
end

function Validation.rules(::Type{<:ConductorGroup})
    (Validation.OwnerRule(:conductor_group_structure, _check_conductor_group),)
end

function ConductorGroup(layer::AbstractConductorPart{T}) where {T <: Real}
    wires = _wire_count(layer)
    return validate(ConductorGroup{T}(
        layer.r_in,
        layer.r_ex,
        layer.cross_section,
        wires,
        _turn_count(layer),
        layer.resistance,
        layer.material_props.alpha,
        layer.gmr,
        AbstractConductorPart{T}[layer]
    ))
end

"""
$(TYPEDSIGNATURES)

Append one concentric conductor layer. The layer must use the group's scalar
type and material reference temperature.
"""
function add!(group::ConductorGroup{T}, layer::AbstractConductorPart{T}) where {T}
    isapprox(layer.r_in, group.r_ex) || throw(DomainError(
        layer.r_in,
        "new conductor layer must start at radius $(group.r_ex)"
    ))
    isapprox(layer.material_props.T0, _reference_temperature(group.layers)) ||
        throw(ArgumentError("new conductor material has a different reference temperature"))

    next_gmr = equivalent_gmr(group, layer)
    next_alpha = equivalent_alpha(
        group.alpha,
        group.resistance,
        layer.material_props.alpha,
        layer.resistance
    )
    next_resistance = parallel(group.resistance, layer.resistance)
    next_cross_section = group.cross_section + layer.cross_section
    added_wires = _wire_count(layer)
    next_wires = group.num_wires + added_wires
    next_turns = iszero(added_wires) ? group.num_turns :
                 (group.num_wires * group.num_turns + added_wires * _turn_count(layer)) /
                 next_wires

    candidate = validate(ConductorGroup{T}(
        group.r_in,
        layer.r_ex,
        next_cross_section,
        next_wires,
        next_turns,
        next_resistance,
        next_alpha,
        next_gmr,
        AbstractConductorPart{T}[group.layers; layer]
    ))
    group.r_ex = candidate.r_ex
    group.cross_section = candidate.cross_section
    group.num_wires = candidate.num_wires
    group.num_turns = candidate.num_turns
    group.resistance = candidate.resistance
    group.alpha = candidate.alpha
    group.gmr = candidate.gmr
    group.layers = candidate.layers
    return group
end

function add!(group::ConductorGroup{T}, layer::AbstractConductorPart{U}) where {T, U}
    throw(ArgumentError(
        "cannot add $(typeof(layer)) to ConductorGroup{$T}; explicitly convert the " *
        "complete design before mutation",
    ))
end

function Base.convert(::Type{ConductorGroup{T}}, group::ConductorGroup) where {T <: Real}
    layers = AbstractConductorPart{T}[convert(AbstractConductorPart{T}, layer)
                                      for layer in group.layers]
    return validate(ConductorGroup{T}(
        convert(T, group.r_in), convert(T, group.r_ex),
        convert(T, group.cross_section), group.num_wires,
        convert(T, group.num_turns), convert(T, group.resistance),
        convert(T, group.alpha), convert(T, group.gmr), layers
    ))
end

Base.convert(::Type{ConductorGroup{T}}, group::ConductorGroup{T}) where {T <: Real} = group
