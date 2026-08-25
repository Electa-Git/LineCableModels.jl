"""
    LineCableModels.ImportExport

Read and write LineCableModels data.

# Overview

- Save and load material and cable libraries.
- Encode package-owned model types in the versioned JSON schema.
- Export cable systems and calculated matrices to PSCAD, ATPDraw, and TRALIN
  formats.
- Import supported PSCAD and TRALIN results.

# Dependencies

$(IMPORTS)

"""
module ImportExport

export export_data
export import_data
export save
export load!

using DocStringExtensions: IMPORTS, TYPEDSIGNATURES, METHODLIST, FUNCTIONNAME,
                           TYPEDEF, TYPEDFIELDS
import ..LineCableModels: add!, validate, nominal, standard_uncertainty
import ..LineCableModels: retired_legacy_json
import ..Grammar: observe
import ..ReportBuilder
using ..Materials: Material, MaterialsLibrary
using ..EarthProps: EarthModel
using ..DataModel: CablesLibrary, CableDesign, CableComponent, ConductorGroup,
                   InsulatorGroup, CircStrands, RectStrands, Strip, Tubular, Semicon,
                   Insulator, CablePosition,
                   LineCableSystem, NominalData
import ..Engine: LineParameters, SeriesImpedance, ShuntAdmittance,
                 frequencies, Z, Y, C
using EzXML
using Dates
using Printf # For ATP export
using JSON3
using Serialization # Read and write the .jls format.
using LinearAlgebra

include("interfaces.jl")
include("paths.jl")
include("serialize.jl")
include("deserialize.jl")
include("cableslibrary.jl")
include("materialslibrary.jl")
include("pscad/pscad.jl")
include("atp.jl")
include("tralin.jl")

end # module ImportExport
