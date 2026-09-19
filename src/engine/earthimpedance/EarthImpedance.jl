"""
    LineCableModels.Engine.EarthImpedance

Define earth-return impedance recipes, numerical primitives, and formula-owned
frequency functors.

# Dependencies

$(IMPORTS)

"""
module EarthImpedance
import ...Grammar: FormulationOptions

# Export public API
export Formula, formula_id, earth_impedance, assumptions, formulas

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ...LineCableModels: validate
import ..Engine: EarthPair, earth_parameters, layer_index
import ...Earth: EquivalentHomogeneous
import ..Engine: EarthImpedanceFormulation, formula_id
#! explicit-imports: off
# Explicitly included equations share these physical and numerical operations.
import ..Engine: earth_bindings, initialize_buffers, earth!,
                 same_physical_state, numerical_magnitude,
                 special_besselix, special_besselkx,
                 special_besseljx
using LinearAlgebra: lu!, ldiv!
import ...LineCableModels: FormulaDefinition, FormulaMethod, nominal
import ..Engine: SpectralIntegral, integrate
import ..Engine: description, conductivity, media, special_besselk
import ..Engine: formulation_options
#! explicit-imports: on

include("interface.jl")

#! explicit-imports: off
const FORMULAS = (
    include("formulas/unified.jl"),
    include("formulas/ametani2009.jl"),
    include("formulas/carson1926.jl"),
    include("formulas/default.jl"),
    include("formulas/gary1976.jl"),
    include("formulas/lucca1994.jl"),
    include("formulas/pollaczek1926.jl"),
    include("formulas/saad1996.jl"),
    include("formulas/wedepohl1973.jl"),
    include("formulas/wise1934.jl"),
    include("formulas/xue2018.jl")
)
#! explicit-imports: on

"""
Return registered earth-impedance identities, including unimplemented equations.
"""
formulas() = FORMULAS

end # module EarthImpedance
