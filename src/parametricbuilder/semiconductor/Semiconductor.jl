"""
    LineCableModels.ParametricBuilder.Semiconductor

Construct semiconducting Regions while enforcing `Material.kind == :semicon`.
"""
module Semiconductor

using ..ParametricBuilder: AbstractGrid, Grid, Gridspace
import ...DataModel
import ...Materials

function _region(tag, primitive, material)
    material isa Materials.Material ||
        throw(ArgumentError("semiconductor material must resolve to Material"))
    material.kind === :semicon || throw(ArgumentError(
        "semiconductor material must have kind :semicon; got :$(material.kind)"
    ))
    return DataModel.Region(tag, primitive, material)
end

function Shell(tag::Symbol, material; t, combine::Symbol = :product)
    values = (tag, t, material)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        build = (resolved_tag, thickness, resolved_material) ->
                _region(resolved_tag, DataModel.Shell(thickness), resolved_material)
        return Gridspace{DataModel.Region}(build, grids; combine)
    end
    return _region(tag, DataModel.Shell(t), material)
end

end
