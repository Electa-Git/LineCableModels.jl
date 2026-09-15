"""
    LineCableModels.Earth.EquivalentHomogeneous

Define equivalent homogeneous-earth rules used when a homogeneous earth-return
formulation consumes a horizontally layered [`EarthModel`](@ref).

Material frequency dependence and equivalent-earth reduction remain separate
operations. [`AfterFD`](@ref) and [`BeforeFD`](@ref) select their composition
order through dispatch.

# Dependencies

$(IMPORTS)
"""
module EquivalentHomogeneous
import ...Grammar: computation_options
import ...LineCableModels: validate

export Formula, AfterFD, BeforeFD
export formula_id, formulas, rule

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
using ..Earth: EarthMaterial, EarthModel
import ...Grammar: AbstractFormulation
import ...LineCableModels: FormulaDefinition, FormulaMethod, formula_id
#! explicit-imports: off
import ...LineCableModels: description
#! explicit-imports: on

include("interface.jl")

#! explicit-imports: off
const FORMULAS = (
    include("formulas/default.jl"),
)
#! explicit-imports: on

"""
Return the built-in equivalent homogeneous-earth formula identifiers.
"""
formulas() = FORMULAS

end # module EquivalentHomogeneous
