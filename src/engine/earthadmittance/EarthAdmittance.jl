"""
    LineCableModels.Engine.EarthAdmittance

Define earth-return admittance recipes, numerical primitives, and
formula-owned frequency functors.

# Dependencies

$(IMPORTS)

"""
module EarthAdmittance

# Export public API
export Formula, formula_id, earth_potential_coefficient, assumptions, propagation, formulas,
       Γ

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ...LineCableModels: validate
import ..Engine: EarthPair, hooks, earth_parameters
import ...Earth: EquivalentHomogeneous
import ..Engine: EarthAdmittanceFormulation, formula_id
#! explicit-imports: off
# Formula files are discovered dynamically; their imports are verified by the
# equation-ownership tests because the static import scanner cannot follow them.
import ..Engine: system_earth, unified_entry, retained_earth_features
import ...LineCableModels: FormulaDefinition, FormulaMethod, nominal
import ..Engine: Formulation, SpectralIntegral, integrate
import ..Engine: description, conductivity, media, special_besselk
import ..Engine: computation_options, LineCableModelsCoaxial
#! explicit-imports: on

vacuum_permittivity(value) = one(value) * 88541878128 * (one(value) * 10)^(-22)
vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("interface.jl")
include("homogeneous.jl")

#! explicit-imports: off
const FORMULAS = (
    include("formulas/default.jl"),
    include("formulas/pollaczek1926.jl"),
    include("formulas/wise1948.jl"),
    include("formulas/xue2018.jl"),
)
#! explicit-imports: on

"""
Return numerical earth-admittance identifiers.
"""
formulas() = FORMULAS

end # module EarthAdmittance
