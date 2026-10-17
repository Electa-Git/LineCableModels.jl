"""
    LineCableModels.Engine.InsulationImpedance

Define registered series-impedance formulas for cable insulation.

# Dependencies

$(IMPORTS)

"""
module InsulationImpedance
import ...Commons: FormulationOptions, formulas
using ...Commons: Functor
import ...Commons: formulation_options

# Export public API
export Formula, formula_id, formulas

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: InsulationImpedanceFormulation, formula_id
import ...LineCableModels: FormulaDefinition, Expression
#! explicit-imports: off
import ..Engine: description
using ...Commons: vacuum_permeability
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
formulas(::Type{<:Formula}) = FORMULAS

end # module InsulationImpedance
