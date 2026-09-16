"""
    LineCableModels.Earth.FrequencyDependent

Define measured and material-physics relations that map one soil material and
one frequency to its frequency-dependent electromagnetic properties.

The explicit `:constant` formula is the frequency-independent pass-through.
The `:default` formula is a routing alias for `:constant`; the remaining
registered identifiers implement literature-based frequency-dispersive laws.

# Dependencies

$(IMPORTS)
"""
module FrequencyDependent
import ...Grammar: FormulationOptions
import ...Grammar: formulation_options
import ...LineCableModels: FormulaDefinition

export Formula, formula_id, formulas, assumptions
public FrequencyDependentFormulation, earth_material

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
using ..Earth: EarthMaterial
import ...Grammar: AbstractFormulation
import ...LineCableModels: FormulaMethod, constitutive, formula_id, validate
#! explicit-imports: off
import ...LineCableModels: description
#! explicit-imports: on

include("interface.jl")

#! explicit-imports: off
const FORMULAS = (
    include("formulas/alipio2014.jl"),
    include("formulas/cigre2019.jl"),
    include("formulas/constant.jl"),
    include("formulas/datsios2019.jl"),
    include("formulas/default.jl"),
    include("formulas/longmire1975.jl"),
    include("formulas/messier1985.jl"),
    include("formulas/portela1999.jl"),
    include("formulas/scott1967.jl"),
    include("formulas/visacro1987.jl"),
    include("formulas/visacro2012.jl"),
)
#! explicit-imports: on

"""
Return the built-in frequency-dependent earth-material formula identifiers.
"""
formulas() = FORMULAS

end # module FrequencyDependent
