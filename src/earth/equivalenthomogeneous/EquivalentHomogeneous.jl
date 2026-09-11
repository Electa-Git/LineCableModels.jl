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
const FORMULAS = let
    directory = joinpath(@__DIR__, "formulas")
    Base.include_dependency(directory)
    files = sort!(filter(
        path -> endswith(path, ".jl"),
        readdir(directory; join = true)
    ))
    identifiers = Symbol[]
    for file in files
        identifier = include(file)
        identifier isa Symbol || error(
            "EquivalentHomogeneous formula file $(basename(file)) must return its Symbol identifier"
        )
        identifier in identifiers && error(
            "duplicate EquivalentHomogeneous formula identifier :$identifier"
        )
        push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"Return the built-in equivalent homogeneous-earth formula identifiers."
formulas() = FORMULAS

end # module EquivalentHomogeneous
