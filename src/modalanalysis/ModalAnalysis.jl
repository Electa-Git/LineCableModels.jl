"""
    LineCableModels.ModalAnalysis

Transform fully coupled line-parameter matrices between phase and modal
coordinate domains independently of the backend that calculated them.

# Dependencies

$(IMPORTS)
"""
module ModalAnalysis

export ModalAnalysisProblem, ModalAnalysisFormulation
export LineCableModelsModal, ModalOperators, Formula
export operators, Tv, Ti, gamma, alpha, beta, velocity, Zc, Yc, H, PropagationParameters, transform
export formula_id, formulas

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: FormulaMethod, nominal, FormulaDefinition, formula, parameterize, validate
import ..LineCableModels: line_length
import ..LineCableModels
import ..Grammar: AbstractProblemDefinition, AbstractFormulation,
                  FormulationOptions, ComputationOptions, ComputationDetails,
                  compute, computation_options, computation_details, formulation_options, details
import ..Grammar: observe, observables, request_identity, request_indices, observation_indices
import ..Engine: LineParameters, LineParametersFormulation, PhaseDomain, ModalDomain,
                 SeriesImpedance, ShuntAdmittance, basis, frequencies,
                 description, formula_id, selectdomain, selectdetails, initialize_buffers
using LinearAlgebra: Diagonal, I, checksquare, cond, diag, dot, eigen,
                     issuccess, ldiv!, lu!, mul!, norm, rdiv!
import ..Grammar: AbstractCoreResult
import ..Engine
import ..Grammar
import ..Units
import ..TextDisplay
#! explicit-imports: on

include("interfaces.jl")
include("problems.jl")

include("formulations.jl")
public decompose!
include("eigensystems.jl")
include("compute.jl")

#! explicit-imports: off
const FORMULAS = (
    include("formulas/chrysochos2014.jl"),
    include("formulas/default.jl"),
)
#! explicit-imports: on

"""
Return the built-in modal-transformation formula identifiers.
"""
formulas() = FORMULAS

include("quantities.jl")
include("propagation.jl")
include("observations.jl")
include("textdisplay.jl")

end # module ModalAnalysis
