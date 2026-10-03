"""
    LineCableModels.ParametricBuilder.Insulator

Construct insulating Regions while enforcing `Material.kind == :insulator`.
"""
module Insulator

import ..ParametricBuilder: insulation

function Shell(tag::Symbol, material; t, combine::Symbol = :product)
    return insulation(material; t, tag, combine)
end

end
