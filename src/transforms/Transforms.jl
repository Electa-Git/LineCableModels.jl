"""
    LineCableModels.Transforms

Transform fully coupled line-parameter matrices between phase and modal
coordinate domains independently of the backend that calculated them.

# Dependencies

$(IMPORTS)
"""
module Transforms

export ModalTransformationProblem, ModalTransformationFormulation
export LineCableModelsModal, ModalOperators, Formula
export operators, formula_id, formulas, gamma, modal_quantities

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: FormulaMethod, nominal, FormulaDefinition, formula, parameterize, validate
import ..Grammar: AbstractProblemDefinition, AbstractFormulation,
                  FormulationOptions, ComputationOptions, ComputationDetails,
                  compute, computation_options, computation_details, formulation_options, details
import ..Engine: LineParameters, LineParametersFormulation, PhaseDomain, ModalDomain,
                 SeriesImpedance, ShuntAdmittance,
                 description, formula_id, selectdomain,
                 offdiagonal_ratio
using LinearAlgebra: Diagonal, I, checksquare, cond, diag, dot, eigen,
                     issuccess, ldiv!, lu!, mul!, norm, rdiv!
#! explicit-imports: on

include("interfaces.jl")
include("problems.jl")

include("formulations.jl")
public modal_operators
include("eigensystems.jl")

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

include("compute.jl")
include("quantities.jl")

end # module Transforms
