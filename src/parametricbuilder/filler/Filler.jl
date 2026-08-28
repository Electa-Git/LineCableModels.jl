"""
    LineCableModels.ParametricBuilder.Filler

Construct explicit filling Regions without imposing a material class.
"""
module Filler

using ..ParametricBuilder: AbstractGrid, Grid, Gridspace
import ...DataModel
import ...Materials

function Region(tag::Symbol, material; primitive, combine::Symbol = :product)
    values = (tag, primitive, material)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{DataModel.Region}(DataModel.Region, grids; combine)
    end
    material isa Materials.Material ||
        throw(ArgumentError("filler material must resolve to Material"))
    return DataModel.Region(values...)
end

end
