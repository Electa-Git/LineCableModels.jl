"""
	LineCableModels.ImportExport

The [`ImportExport`](@ref) module provides methods for serializing and deserializing data structures in [`LineCableModels.jl`](index.md), and data exchange with external programs.

# Overview

This module provides functionality for:

- Saving and loading cable designs and material libraries to/from JSON and other formats.
- Exporting cable system models to PSCAD format.
- Serializing custom types with special handling for measurements and complex numbers.

The module implements a generic serialization framework with automatic type reconstruction
and proper handling of Julia-specific types like `Measurement` objects and `Inf`/`NaN` values.

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module ImportExport

# Export public API
export export_data
export read_data
export save
export load!
export get
export delete!

# Load common dependencies
using ..LineCableModels
include("commondeps.jl")

# Module-specific dependencies
using Measurements
using EzXML # For PSCAD export
using Dates # For PSCAD export
using Printf # For ATP export
using JSON3
using Serialization # For .jls format
import Base: get
using ..Utils
using ..Materials
using ..EarthProps
using ..DataModel
import ..LineCableModels: add!, load!, export_data, save, _is_headless, _display_path, _CLEANMETHODLIST


include("ImportExport/serialize.jl")
include("ImportExport/deserialize.jl")
include("ImportExport/cableslibrary.jl")
include("ImportExport/materialslibrary.jl")
include("ImportExport/pscad.jl")
include("ImportExport/atp.jl")

end # module ImportExport
