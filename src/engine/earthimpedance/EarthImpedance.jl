"""
    LineCableModels.Engine.EarthImpedance

Define earth-return impedance recipes, numerical primitives, and formula-owned
frequency functors.

# Dependencies

$(IMPORTS)

"""
module EarthImpedance
import ...Commons: FormulationOptions, bindings, formulas

# Export public API
export Formula, formula_id, earth_impedance, formulas

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ...LineCableModels: validate
import ..Engine: EarthPair, layer_index
import ...Earth: EquivalentHomogeneous
import ..Engine: EarthImpedanceFormulation, formula_id
#! explicit-imports: off
# Explicitly included equations share these physical and numerical operations.
import ..Engine: earth!
import ...LineCableModels: FormulaDefinition, FormulaMethod
import ..Engine: description, conductivity, special_besselk
import ..Engine: formulation_options
import ..Engine: earth_spectral_term, earth_direct
using ...Commons: vacuum_permeability
#! explicit-imports: on

public Functor, axial_field_coefficient

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
formulas(::Type{<:Formula}) = FORMULAS

end # module EarthImpedance
