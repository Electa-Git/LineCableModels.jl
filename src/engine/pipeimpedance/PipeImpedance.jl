"""
    LineCableModels.Engine.PipeImpedance

Own pipe-type formula selections and backend applicability. No analytical
pipe-type implementation is supplied yet.
"""
module PipeImpedance
import ...Commons: FormulationOptions, formulas, formulation_options
import ...LineCableModels: FormulaDefinition

export Formula, formula_id, formulas

using DocStringExtensions: TYPEDEF
#! explicit-imports: off
# Expanded in the docstrings of the included formula files.
using DocStringExtensions: TYPEDSIGNATURES
#! explicit-imports: on
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
Return registered pipe selections, including the explicit default selection.
"""
formulas(::Type{<:Formula}) = FORMULAS

end
