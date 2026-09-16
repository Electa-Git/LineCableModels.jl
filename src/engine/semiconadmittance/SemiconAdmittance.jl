"""
    LineCableModels.Engine.SemiconAdmittance

Define registered constitutive relations for semiconducting-screen admittance.

# Dependencies

$(IMPORTS)

"""
module SemiconAdmittance
import ...Grammar: FormulationOptions
import ...Grammar: formulation_options

export Formula, formula_id, formulas

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: SemiconAdmittanceFormulation, formula_id, validate
import ...LineCableModels: FormulaDefinition, FormulaMethod
using ...Materials: Material
#! explicit-imports: off
import ..Engine: description, conductivity
#! explicit-imports: on

include("interface.jl")

public semicon_material

#! explicit-imports: off
const FORMULAS = (
    include("formulas/default.jl"),
    include("formulas/lossless.jl"),
    include("formulas/lossy.jl"),
)
#! explicit-imports: on

"""
Return the built-in semicon-admittance formula identifiers.
"""
formulas() = FORMULAS

end # module SemiconAdmittance
