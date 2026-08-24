"""
    LineCableModels.ParametricBuilder.Insulator

Declare tubular insulation and semiconducting parts for
[`LineCableModels.ParametricBuilder.CableBuilder`](@ref). An explicit `Grid` or
`Gridspace` input returns a finite space of part builders.
"""
module Insulator

using ..ParametricBuilder:
                           AbstractGrid, Grid, Gridspace, PartBuilder,
                           _radial_declaration
import ...DataModel

function _annular(
        part,
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
        Val(:insulator), Val(part), mode, component, layers, dimension, material
    )
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PartBuilder}(PartBuilder, grids; combine)
    end
    return PartBuilder(values...)
end

Tubular(component::Symbol; kwargs...) = _annular(DataModel.Insulator, component; kwargs...)
Semicon(component::Symbol; kwargs...) = _annular(DataModel.Semicon, component; kwargs...)

end
