module LineCableModels

## Public API
# -------------------------------------------------------------------------
# Core generics:
export add!, set_verbosity!, set_backend!

# Materials:
export Material, MaterialsLibrary

# Data model (design + system):
export Thickness, Diameter, CircStrands, RectStrands, Strip, Tubular, Semicon, Insulator
export ConductorGroup, InsulatorGroup
export CableComponent, CableDesign, NominalData
export CablesLibrary
export CablePosition, LineCableSystem
export trifoil_formation, flat_formation, preview, equivalent, MaxFill

# Earth properties:
export EarthModel

# Engine:
export LineParametersProblem,
	FormulationSet, compute!, SeriesImpedance, ShuntAdmittance, per_km, per_m, kronify,
	LineParameters, PhaseDomain, ModalDomain

# Parametric builder:
# export make_stranded, make_screened
# export conductor, insulator

# Import/Export:
export export_data, save, load!
# -------------------------------------------------------------------------

import DocStringExtensions: DocStringExtensions

# Submodule `Commons`
include("commons/Commons.jl")
using .Commons: IMPORTS, EXPORTS, add!, PhaseDomain, ModalDomain, domain
# Submodule `UncertainBessels`
include("uncertainbessels/UncertainBessels.jl")

# Submodule `Utils`
include("utils/Utils.jl")
using .Utils: set_verbosity!

# Submodule `UnitHandler`
include("unithandler/UnitHandler.jl")

# Submodule `Validation`
include("validation/Validation.jl")

# Submodule `PlotBuilder`
include("plotbuilder/PlotBuilder.jl")
using .PlotBuilder.BackendHandler: set_backend!

# Submodule `Materials`
include("materials/Materials.jl")
using .Materials: Material, MaterialsLibrary

# Submodule `EarthProps`
include("earthprops/EarthProps.jl")
using .EarthProps: EarthModel

# Submodule `DataModel`
include("datamodel/DataModel.jl")
using .DataModel: Thickness, Diameter, CircStrands, RectStrands, Strip, Tubular, Semicon,
	Insulator, ConductorGroup, InsulatorGroup, CableComponent, CableDesign, NominalData,
	CablesLibrary, CablePosition, LineCableSystem, trifoil_formation, flat_formation,
	preview, equivalent, MaxFill

# Submodule `Engine`
include("engine/Engine.jl")
using .Engine: LineParametersProblem, compute!, LineParameters, SeriesImpedance,
	ShuntAdmittance, per_km, per_m, kronify, FormulationSet

# Submodule `ParametricBuilder`
include("parametricbuilder/ParametricBuilder.jl")

# Submodule `UQ`
include("uq/UQ.jl")

# Submodule `ImportExport`
include("importexport/ImportExport.jl")
using .ImportExport: export_data, load!, save

# Aliases for backward compatibility
const WireArray = CircStrands  # alias for now
export WireArray  # export aliases

end