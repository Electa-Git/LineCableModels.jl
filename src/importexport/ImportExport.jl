"""
	LineCableModels.ImportExport

The [`ImportExport`](@ref) module provides methods for serializing and deserializing data structures in [`LineCableModels.jl`](index.md), and data exchange with external programs.

# Overview

This module provides functionality for:

- Saving and loading cable designs and material libraries to/from JSON and other formats.
- Exporting cable system models to PSCAD and ATP formats.
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


# Module-specific dependencies
using ..Commons
using ..Utils: display_path, to_nominal, resolve_T, coerce_to_T
using ..Materials: Material, MaterialsLibrary
using ..EarthProps: EarthModel
using ..DataModel: CablesLibrary, CableDesign, CableComponent, ConductorGroup,
	InsulatorGroup, WireArray, Strip, Tubular, Semicon, Insulator, LineCableSystem,
	NominalData
import ..Engine: LineParameters, SeriesImpedance, ShuntAdmittance
using Measurements
using EzXML
using Dates
using Printf # For ATP export
using JSON3
using Serialization # For .jls format
using LinearAlgebra


"""
$(TYPEDSIGNATURES)

Export [`LineCableModels`](@ref) data for use in different EMT-type programs.

# Methods

$(METHODLIST)
"""
# function export_data end
export_data(backend::Symbol, args...; kwargs...) =
	export_data(Val(backend), args...; kwargs...)

include("serialize.jl")
include("deserialize.jl")
include("cableslibrary.jl")
include("materialslibrary.jl")
include("pscad.jl")
include("atp.jl")
include("tralin.jl")

end # module ImportExport
