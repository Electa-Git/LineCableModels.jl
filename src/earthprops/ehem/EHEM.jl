"""
    LineCableModels.EarthProps.EHEM

Define equivalent homogeneous-earth rules used when a homogeneous earth-return
formulation consumes a horizontally layered [`EarthModel`](@ref).

Material frequency dependence and equivalent-earth reduction remain separate
operations. [`AfterFD`](@ref) and [`BeforeFD`](@ref) select their composition
order through dispatch.

# Dependencies

$(IMPORTS)
"""
module EHEM

export Formula, Layer, AfterFD, BeforeFD
export formula_id, assumptions, formulas, rule

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
using ..EarthProps: EarthMaterial, EarthModel
import ...Grammar: AbstractFormulation
import ...LineCableModels: FormulaMethod, formula_id
#! explicit-imports: off
import ...LineCableModels: description
#! explicit-imports: on

"Last-layer reconstruction policy selected by `:default`."
const DEFAULT = :Layer

include("interface.jl")
include("layer.jl")

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
            "EHEM formula file $(basename(file)) must return its Symbol identifier"
        )
        identifier in identifiers && error(
            "duplicate EHEM formula identifier :$identifier"
        )
        push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"Return the built-in equivalent homogeneous-earth formula identifiers."
formulas() = FORMULAS

end # module EHEM
