"""
    PSCAD

PSCAD model exchange and native line-parameter computation. Loading this module
defines the backend without contacting PSCAD. Native runs connect to PSCAD's automation
server through the `LineCableModelsPSCADExt` extension, which `using PythonCall` loads.
"""
module PSCAD

using SHA: sha256
import TOML
using LineCableModels
using LineCableModels.DataModel: LineCableSystem
using LineCableModels.Earth: EarthModel
using LineCableModels.Commons: vacuum_permittivity
using LineCableModels.Engine
using LineCableModels.ImportExport
import LineCableModels: description, parameterize, computation_details, validate,
                        Expression
import LineCableModels.Engine.EarthImpedance: earth_impedance
import LineCableModels.Engine.EarthAdmittance: earth_potential_coefficient
import LineCableModels.Engine.InternalImpedance: internal_impedance
import LineCableModels.Engine.InsulationImpedance: insulation_impedance
import LineCableModels.Engine: Formulation,
                               LineParametersProblem
import LineCableModels.Commons: AbstractFormulation, ComputationOptions, ComputationDetails,
                                FormulationOptions, computation_options, compute,
                                formulation_options, gridpoint_id, resolve
using LineCableModels.Commons: Resolved
using DocStringExtensions: TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
import LineCableModels: constitutive
import LineCableModels.DataModel
using LineCableModels.DataModel.BaseParams: equivalent_eps
import LineCableModels.Engine
import LineCableModels.Units
import LineCableModels.ImportExport: import_data, export_data
import EzXML
using EzXML: ElementNode, XMLDocument, addelement!, nodename,
             readxml, root, setroot!

export PSCADFormulation, RemoteConfig,
       read_pscad_result, run_remote_pscad, identify
public pscad_setting, PSCAD_VERSION

include("formulations.jl")
include("importexport/pscad.jl")
include("results.jl")
include("validate.jl")
include("remote/configuration.jl")
include("remote/remote.jl")
include("compute.jl")

end
