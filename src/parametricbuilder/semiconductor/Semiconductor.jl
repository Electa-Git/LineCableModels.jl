"""
    LineCableModels.ParametricBuilder.Semiconductor

Construct semiconducting Regions while enforcing `Material.kind == :semicon`.
"""
module Semiconductor

import ..ParametricBuilder: screen

function Shell(tag::Symbol, material; t, combine::Symbol = :product)
    return screen(material; t, tag, combine)
end

end
