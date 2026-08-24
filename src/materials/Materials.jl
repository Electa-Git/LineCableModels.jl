"""
    LineCableModels.Materials

Define electromagnetic material records and an in-memory material library.

# Public actions

- Construct and validate [`Material`](@ref) values.
- Add, remove, and retrieve materials in a [`MaterialsLibrary`](@ref).
- Present material data through Base and DataFrames protocols.

JSON persistence belongs to `LineCableModels.ImportExport`.

# Dependencies

$(IMPORTS)
"""
module Materials

export Material, MaterialsLibrary, add!

using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES,
                           FUNCTIONNAME
import ..LineCableModels: add!, validate
import ..Validation

include("material.jl")
include("materialslibrary.jl")
include("dataframe.jl")
include("base.jl")

end # module Materials
