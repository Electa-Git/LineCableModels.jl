"""
    LineCableModels.Engine.EarthAdmittance

Define earth-return admittance recipes, numerical primitives, and
formula-owned frequency functors.

# Dependencies

$(IMPORTS)

"""
module EarthAdmittance
import ...Grammar: FormulationOptions

# Export public API
export Formula, formula_id, earth_potential_coefficient, assumptions, formulas

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ...LineCableModels: validate
import ..Engine: EarthPair, earth_parameters, layer_index
import ...Earth: EquivalentHomogeneous
import ..Engine: EarthAdmittanceFormulation, formula_id
#! explicit-imports: off
# Explicitly included equations share these physical and numerical operations.
import ..Engine: earth_bindings, initialize_buffers, earth!, same_physical_state
import ..Engine: computation_type
import ..EarthImpedance
import ...LineCableModels: FormulaDefinition, FormulaMethod, nominal
import ..Engine: description, conductivity, media
import ..Engine: formulation_options
#! explicit-imports: on

include("interface.jl")

#! explicit-imports: off
const FORMULAS = (
    include("formulas/unified.jl"),
    include("formulas/default.jl"),
    include("formulas/pollaczek1926.jl"),
    include("formulas/wise1948.jl"),
    include("formulas/xue2018.jl")
)
#! explicit-imports: on

"""
Return registered earth-admittance identities, including unimplemented equations.
"""
formulas() = FORMULAS

end # module EarthAdmittance
