"""
    LineCableModels.Engine.InsulationAdmittance

Define registered constitutive relations for cable-insulation admittance.
`:lossy` retains material conduction and displacement current; `:lossless`
explicitly selects the lossless approximation; and `:default` routes to
`:lossless`.

# Dependencies

$(IMPORTS)

"""
module InsulationAdmittance
import ...Grammar: computation_options

# Export public API
export Formula, formula_id, formulas

# Module-specific dependencies
#! explicit-imports: off
# IMPORTS is expanded in this module docstring rather than called as Julia code.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: InsulationAdmittanceFormulation, formula_id
import ...LineCableModels: FormulaDefinition, FormulaMethod
using ...Materials: Material
#! explicit-imports: off
import ..Engine: description, conductivity
#! explicit-imports: on

include("interface.jl")

#! explicit-imports: off
const FORMULAS = (
    include("formulas/default.jl"),
    include("formulas/lossless.jl"),
    include("formulas/lossy.jl"),
)
#! explicit-imports: on

"""
Return the built-in insulation-admittance formula identifiers.
"""
formulas() = FORMULAS

end # module InsulationAdmittance
