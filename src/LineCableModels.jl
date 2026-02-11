module LineCableModels

## Public API
# -------------------------------------------------------------------------
# Core generics:
export add!, set_verbosity!, set_backend!

# Materials:
export Material, MaterialsLibrary

# Data model (design + system):
export Thickness, Diameter, WireArray, Strip, Tubular, Semicon, Insulator, Sector, SectorParams, SectorInsulator
export ConductorGroup, InsulatorGroup
export CableComponent, CableDesign, NominalData
export CablesLibrary
export CablePosition, LineCableSystem
export trifoil_formation, flat_formation, preview, equivalent

# Earth properties:
export EarthModel

# Engine:
export LineParametersProblem,
	FormulationSet,
	compute!, SeriesImpedance, ShuntAdmittance, per_km, per_m, kronify

# Import/Export:
export export_data, save, load!
# -------------------------------------------------------------------------

import DocStringExtensions: DocStringExtensions

# Submodule `Commons`
include("commons/Commons.jl")
using .Commons: IMPORTS, EXPORTS, add!

# Submodule `UncertainBessels`
include("uncertainbessels/UncertainBessels.jl")

# Submodule `Utils`
include("utils/Utils.jl")
using .Utils: set_verbosity!

# Submodule `BackendHandler`
include("backendhandler/BackendHandler.jl")
using .BackendHandler: set_backend!

# Submodule `PlotUIComponents`
include("plotuicomponents/PlotUIComponents.jl")

# Submodule `Validation`
include("validation/Validation.jl")

# Submodule `Materials`
include("materials/Materials.jl")
using .Materials: Material, MaterialsLibrary

# Submodule `EarthProps`
include("earthprops/EarthProps.jl")
using .EarthProps: EarthModel

# Submodule `DataModel`
include("datamodel/DataModel.jl")
using .DataModel: Thickness, Diameter, WireArray, Strip, Tubular, Semicon, Insulator, Sector, SectorParams, SectorInsulator,
	ConductorGroup, InsulatorGroup, CableComponent, CableDesign, NominalData, CablesLibrary,
	CablePosition, LineCableSystem, trifoil_formation, flat_formation, preview, equivalent

# Submodule `Engine`
include("engine/Engine.jl")
using .Engine: LineParametersProblem, compute!, SeriesImpedance, ShuntAdmittance, per_km,
	per_m, kronify, FormulationSet

# Submodule `ImportExport`
include("importexport/ImportExport.jl")
using .ImportExport: export_data, load!, save

end