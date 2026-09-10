"""
    PSCAD

PSCAD model exchange and native line-parameter computation. Loading this module
defines the backend without contacting a station or launching PSCAD.
"""
module PSCAD

using Base64: base64encode
using SHA: sha256
import TOML
using LineCableModels
using LineCableModels.DataModel: LineCableSystem
using LineCableModels.Earth: EarthModel
using LineCableModels.Engine
using LineCableModels.ImportExport
import LineCableModels: description, parameterize, computation_details, validate,
                        FormulaMethod
import LineCableModels.Engine.EarthImpedance: earth_impedance
import LineCableModels.Engine.EarthAdmittance: earth_potential_coefficient
import LineCableModels.Engine.InternalImpedance: internal_impedance
import LineCableModels.Engine.InsulationImpedance: insulation_impedance
import LineCableModels.Engine: Formulation,
                               LineParametersProblem,
                               verbosity
import LineCableModels.Grammar: AbstractFormulation, ComputationOptions,
                                FormulationOptions, computation_options, compute,
                                formulation_options
using DocStringExtensions: TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
import LineCableModels: constitutive
import LineCableModels.DataModel
import LineCableModels.Engine
import LineCableModels.ImportExport: import_data, export_data
import EzXML
using EzXML: ElementNode, XMLDocument, addelement!, nodename,
             readxml, root, setroot!

export PSCADFormulation, RemoteConfig,
       read_pscad_result, remote_command, run_remote_pscad, formulas, identify
public pscad_setting

include("formulations.jl")
include("importexport/pscad.jl")
include("results.jl")
include("validate.jl")
include("remote/configuration.jl")
include("remote/files.jl")
include("remote/remote.jl")
include("compute.jl")

end
