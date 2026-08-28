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

#! explicit-imports: off
# IMPORTS is expanded in the module docstring rather than called as Julia code.
using DocStringExtensions: IMPORTS
#! explicit-imports: on
using DocStringExtensions: TYPEDSIGNATURES, METHODLIST
import ..LineCableModels: build, validate, nominal
import ..Grammar: observe
import ..ReportBuilder
using ..Materials: Material, MaterialsLibrary
using ..EarthProps: EarthLayer, EarthModel
import ..DataModel
using ..DataModel: CablesLibrary, CableDesign, LineCableSystem,
                   AbstractCablePart, Region, Stack, Group, Assembly, Enclosure,
                   DiskDefinition, RectangleDefinition, EllipseDefinition,
                   SectorDefinition, AnnulusDefinition, ShellDefinition,
                   PolygonDefinition, Pose2,
                   Ring, Polar, Fill, Lattice, capacity,
                   FillFactor, DiameterFactor, TabulatedCompaction,
                   AffineCompaction,
                   LayRatio, Pitch, LayAngle, Helix
using ..ParametricBuilder: AbstractGrid, DeterministicGrid, RelativeGrid,
                          AbsoluteGrid, Grid, Gridspace, AbsoluteError
import ..Engine
import ..Engine: LineParameters, SeriesImpedance, ShuntAdmittance,
                 frequencies, Z, Y, C
import EzXML
using EzXML: ElementNode, XMLDocument, addelement!, nodename, prettyprint,
             readxml, root, setroot!
using Printf: @printf, @sprintf
import JSON3
import Serialization
using LinearAlgebra: tril

include("interfaces.jl")
include("paths.jl")
include("serialize.jl")
include("deserialize.jl")
include("cableslibrary.jl")
include("materialslibrary.jl")
include("pscad/pscad.jl")
include("atp.jl")
include("tralin.jl")

public serialize_value, deserialize_value, deserialize_extension

end # module ImportExport
