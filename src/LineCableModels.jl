module LineCableModels

# Package-wide definitions
include("globals.jl")
include("utils/logging.jl")

# Submodule `Utils`
include("Utils.jl")
@force using .Utils

# Submodule `Validation`
include("Validation.jl")
@force using .Validation

# Submodule `Materials`
include("Materials.jl")
@force using .Materials
@reexport using .Materials

# Submodule `EarthProps`
include("EarthProps.jl")
@force using .EarthProps
@reexport using .EarthProps

# Submodule `DataModel`
include("DataModel.jl")
@force using .DataModel
@reexport using .DataModel

# Submodule `Engine`
include("Engine.jl")
@force using .Engine
@reexport using .Engine

# Submodule `ImportExport`
include("ImportExport.jl")
@force using .ImportExport
@reexport using .ImportExport

end