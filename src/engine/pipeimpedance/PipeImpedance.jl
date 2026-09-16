"""
    LineCableModels.Engine.PipeImpedance

Own pipe-type formula selections and backend applicability. No analytical
pipe-type implementation is supplied yet.
"""
module PipeImpedance
import ...Grammar: FormulationOptions
import ...LineCableModels: FormulaDefinition

export Formula, formula_id, formulas

using DocStringExtensions: TYPEDEF, TYPEDSIGNATURES
import ..Engine: PipeImpedanceFormulation, Formulation
import ...LineCableModels: formula_id
#! explicit-imports: off
# These bindings are consumed by the formula definitions below.
import ..Engine: LineCableModelsCoaxial
import ...LineCableModels: description
#! explicit-imports: on
import ...DataModel
using ...DataModel: CableDesign

include("interface.jl")

#! explicit-imports: off
const FORMULAS = (
    include("formulas/none.jl"),
    include("formulas/default.jl"),
)
#! explicit-imports: on

"""
Return registered pipe selections, including the explicit default policy.
"""
formulas() = FORMULAS

end
