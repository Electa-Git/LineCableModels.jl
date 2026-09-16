"""
    LineCableModels.Engine.InsulationImpedance

Define registered series-impedance formulas for cable insulation.

# Dependencies

$(IMPORTS)

"""
module InsulationImpedance
import ...Grammar: FormulationOptions
import ...Grammar: formulation_options

# Export public API
export Formula, formula_id, formulas

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: InsulationImpedanceFormulation, formula_id
import ...LineCableModels: FormulaDefinition, FormulaMethod
#! explicit-imports: off
import ..Engine: description
#! explicit-imports: on

include("interface.jl")

public insulation_impedance

#! explicit-imports: off
const FORMULAS = (
    include("formulas/ametani1980.jl"),
    include("formulas/default.jl"),
)
#! explicit-imports: on

"""
Return the built-in insulation-impedance formula identifiers.
"""
formulas() = FORMULAS

end # module InsulationImpedance
