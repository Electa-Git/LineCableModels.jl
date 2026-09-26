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
import ...Grammar: FormulationOptions
import ...Grammar: formulation_options
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

public equivalent_material, AbstractRule, AbstractSequence

#! explicit-imports: off
const FORMULAS = (
    include("formulas/bottommost.jl"),
    include("formulas/default.jl"),
)
#! explicit-imports: on

"""
Return the built-in equivalent homogeneous-earth formula identifiers.
"""
formulas() = FORMULAS

end # module EquivalentHomogeneous
