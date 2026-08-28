"""
    LineCableModels.ParametricBuilder.Insulator

Construct insulating Regions while enforcing `Material.kind == :insulator`.
"""
module Insulator

using ..ParametricBuilder: AbstractGrid, Grid, Gridspace
import ...DataModel
import ...Materials

function _region(tag, primitive, material)
    material isa Materials.Material ||
        throw(ArgumentError("insulator material must resolve to Material"))
    material.kind === :insulator || throw(ArgumentError(
        "insulator material must have kind :insulator; got :$(material.kind)"
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
                _region(resolved_tag, DataModel.ShellDefinition(thickness), resolved_material)
        return Gridspace{DataModel.Region}(build, grids; combine)
    end
    return _region(tag, DataModel.ShellDefinition(t), material)
end

end
