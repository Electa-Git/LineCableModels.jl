"""
    LineCableModels.ParametricBuilder.Conductor

Declare solid, tubular, wire, strip, and stranded conductor parts for
[`LineCableModels.ParametricBuilder.CableBuilder`](@ref). An explicit `Grid` or
`Gridspace` input returns a finite space of part builders.
"""
module Conductor

using ..ParametricBuilder:
                           AbstractGrid, Grid, Gridspace, PartBuilder,
                           _radial_declaration
import ...DataModel

function Solid(
        component::Symbol;
        radius,
        material,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    values = (
        Val(:conductor), Val(DataModel.Tubular), Val(:solid), component, 1,
        radius, material
    )
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PartBuilder}(PartBuilder, grids; combine)
    end
    return PartBuilder(values...)
end

function Tubular(
        component::Symbol;
        layers = 1,
        radius = nothing,
        thickness = nothing,
        material,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    mode, dimension = _radial_declaration(radius, thickness)
    values = (
        Val(:conductor), Val(DataModel.Tubular), mode, component, layers,
        dimension, material
    )
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PartBuilder}(PartBuilder, grids; combine)
    end
    return PartBuilder(values...)
end

function Wires(
        component::Symbol;
        layers = 1,
        wire_radius,
        num_wires,
        lay_ratio = 11.0,
        material,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    values = (
        Val(:conductor), Val(DataModel.CircStrands), Val(:wire_radius),
        component, layers, wire_radius, num_wires, lay_ratio, material
    )
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PartBuilder}(PartBuilder, grids; combine)
    end
    return PartBuilder(values...)
end

function Strip(
        component::Symbol;
        layers = 1,
        radius = nothing,
        thickness = nothing,
        width,
        lay_ratio = 0.0,
        material,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    mode, dimension = _radial_declaration(radius, thickness)
    values = (
        Val(:conductor), Val(DataModel.Strip), mode, component, layers,
        dimension, width, lay_ratio, material
    )
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PartBuilder}(PartBuilder, grids; combine)
    end
    return PartBuilder(values...)
end

function Stranded(
        component::Symbol;
        layers::Int,
        wire_radius,
        num_wires::Int = 6,
        lay_ratio = 11.0,
        material,
        combine::Symbol = :product
)
    layers >= 1 || throw(ArgumentError("layers must be at least one"))
    central = Wires(
        component;
        layers = 1,
        wire_radius,
        num_wires = 1,
        lay_ratio = 0.0,
        material,
        combine
    )
    layers == 1 && return (central,)
    rings = Wires(
        component;
        layers = layers - 1,
        wire_radius,
        num_wires,
        lay_ratio,
        material,
        combine
    )
    return (central, rings)
end

end
