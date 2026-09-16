"""
    LineCableModels.Materials.TemperatureDependent

Evaluate electrical resistivity at a prescribed temperature from a material's
reference calibration and a selected constitutive equation.

# Dependencies

$(IMPORTS)
"""
module TemperatureDependent
import ...Grammar: FormulationOptions

export Formula, formula_id, formulas
public temperature_resistivity

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
using ..Materials: Material
import ...TextDisplay
import ...Grammar: AbstractFormulation, formulation_options
import ...LineCableModels: FormulaDefinition, FormulaMethod, constitutive,
                          formula_id, validate
#! explicit-imports: off
# Used by the formula files included below.
import ...LineCableModels: description
#! explicit-imports: on

include("interface.jl")

public TemperatureDependentFormulation

#! explicit-imports: off
const FORMULAS = (
    include("formulas/linear.jl"),
    include("formulas/default.jl"),
)
#! explicit-imports: on

"""Return the registered electrical-resistivity temperature laws."""
formulas() = FORMULAS

end # module TemperatureDependent
