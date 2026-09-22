"""
    LineCableModels.Engine.InternalImpedance

Define conductor internal-impedance recipes and their electromagnetic
interaction formulas.

# Dependencies

$(IMPORTS)

"""
module InternalImpedance
import ...Grammar: FormulationOptions

# Export public API
export Formula, formula_id, formulas, internal_impedance, surface_impedances

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: InternalImpedanceFormulation, Formulation, formula_id, formulation_options
#! explicit-imports: off
import ...LineCableModels: FormulaDefinition, FormulaMethod
import ..Engine: description, conductivity
import ..Engine: special_besselix, special_besselkx
#! explicit-imports: on

vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("interface.jl")

#! explicit-imports: off
const FORMULAS = (
    include("formulas/default.jl"),
    include("formulas/schelkunoff1934.jl"),
    include("formulas/wedepohl1973.jl"),
)
#! explicit-imports: on

"""
Return the built-in internal-impedance formula identifiers.
"""
formulas() = FORMULAS

end # module InternalImpedance
