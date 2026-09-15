"""
    LineCableModels.Engine.EarthImpedance

Define earth-return impedance recipes, numerical primitives, and formula-owned
frequency functors.

# Dependencies

$(IMPORTS)

"""
module EarthImpedance

# Export public API
export Formula, formula_id, earth_impedance, assumptions, propagation, formulas, Γ

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ...LineCableModels: validate
import ..Engine: EarthPair, earth_parameters
import ...LineCableModels: constitutive
import ...Earth: EquivalentHomogeneous
import ..Engine: EarthImpedanceFormulation, formula_id
#! explicit-imports: off
# Explicitly included equations share these physical and numerical operations.
import ..Engine: system_earth, unified_entry, retained_earth_features
import ...LineCableModels: FormulaDefinition, FormulaMethod, nominal
import ..Engine: SpectralIntegral, integrate
import ..Engine: description, conductivity, media, special_besselk
import ..Engine: computation_options
#! explicit-imports: on

vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("interface.jl")
include("homogeneous.jl")

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
    include("formulas/xue2018.jl"),
)
#! explicit-imports: on

"""
Return numerical earth-impedance identifiers.
"""
formulas() = FORMULAS

end # module EarthImpedance
