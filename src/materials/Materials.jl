"""
    LineCableModels.Materials

Define electromagnetic material records and an in-memory material library.

# Public actions

- Construct and validate [`Material`](@ref) values.
- Add, remove, and retrieve materials in a [`MaterialsLibrary`](@ref).
- Present material data through the Base display protocol.

JSON persistence belongs to `LineCableModels.ImportExport`.

# Dependencies

$(IMPORTS)
"""
module Materials

export AbstractMaterial, Material, MaterialsLibrary, add!

#! explicit-imports: off
# IMPORTS is expanded in the module docstring rather than called as Julia code.
using DocStringExtensions: IMPORTS
#! explicit-imports: on
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES, FUNCTIONNAME
using RequiredInterfaces: @required
import ..LineCableModels: add!, validate
import ..TextDisplay

include("material.jl")
include("materialslibrary.jl")
include("base.jl")

@required AbstractMaterial begin
    validate(::AbstractMaterial)
end

end # module Materials
