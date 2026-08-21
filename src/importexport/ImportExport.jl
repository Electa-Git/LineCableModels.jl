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

"""
module ImportExport

# Export public API
export export_data
export import_data
export read_data
export save
export load!

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDSIGNATURES, METHODLIST, FUNCTIONNAME,
                           TYPEDEF, TYPEDFIELDS
import ..LineCableModels: add!, validate, nominal, standard_uncertainty
import ..LineCableModels: retired_legacy_json
using ..Materials: Material, MaterialsLibrary
using ..EarthProps: EarthModel
using ..DataModel: CablesLibrary, CableDesign, CableComponent, ConductorGroup,
                   InsulatorGroup, CircStrands, RectStrands, Strip, Tubular, Semicon,
                   Insulator, CablePosition,
                   LineCableSystem, NominalData
import ..Engine: LineParameters, SeriesImpedance, ShuntAdmittance
using EzXML
using Dates
using Printf # For ATP export
using JSON3
using Serialization # For .jls format
using LinearAlgebra
using XLSX
using Tables
using DataFrames

_display_path(path::AbstractString) =
    try
        relpath(abspath(path), pwd())
    catch
        basename(path)
    end

function _json_path(file_name::AbstractString)
    path = isabspath(file_name) ? String(file_name) : abspath(file_name)
    extension = lowercase(splitext(path)[2])
    isempty(extension) && return path * ".json"
    extension == ".json" || throw(ArgumentError(
        "JSON output requires a .json extension; got '$extension'",
    ))
    return path
end

"""
$(TYPEDSIGNATURES)

Export [`LineCableModels`](@ref) data for use in different EMT-type programs.

# Methods

$(METHODLIST)
"""
# function export_data end
function export_data(backend::Symbol, args...; kwargs...)
    export_data(Val(backend), args...; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Import external data through a format-specific backend.

# Arguments

- `backend`: Registered backend name.
- `args`: Backend-specific input arguments.

# Keywords

Backend-specific options.

# Returns

Backend-specific materialized objects.

# Methods

$(METHODLIST)
"""
function import_data(backend::Symbol, args...; kwargs...)
    import_data(Val(backend), args...; kwargs...)
end

include("serialize.jl")
include("deserialize.jl")
include("cableslibrary.jl")
include("materialslibrary.jl")
include("pscad.jl")
include("atp.jl")
include("xlsx.jl")
include("tralin.jl")

end # module ImportExport
