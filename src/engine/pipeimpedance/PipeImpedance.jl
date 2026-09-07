"""
    LineCableModels.Engine.PipeImpedance

Evaluate circular-pipe cavity and inner-wall impedance coefficients.
Individual conductor skin impedances and the pipe's outer-surface/transfer
assembly are separate contributions.

# Dependencies

$(IMPORTS)
"""
module PipeImpedance

export Formula, Pair, formula_id, formulas, routes, assumptions
#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: PipeImpedanceFormulation, formula_id
import ..Engine: special_besselix, special_besselkx
#! explicit-imports: off
import ...LineCableModels: FormulaMethod
import ..Engine: description, InternalImpedance
#! explicit-imports: on
# SpecialFunctions exposes numerical failures through this unexported type.
# Keep the precision retry restricted to that failure, not arbitrary exceptions.
#! explicit-imports: off
import SpecialFunctions: AmosException
#! explicit-imports: on

const DEFAULT = :DaSilva2006
include("interface.jl")
include("circular.jl")
include("proximity.jl")

#! explicit-imports: off
const FORMULAS = let
    directory = joinpath(@__DIR__, "formulas")
    Base.include_dependency(directory)
    identifiers = Symbol[]
    for path in sort(filter(endswith(".jl"), readdir(directory; join=true)))
        id = include(path)
        id isa Symbol || error("pipe formula file must return its Symbol identifier")
        id in identifiers && error("duplicate pipe formula identifier :$id")
        push!(identifiers, id)
    end
    Tuple(identifiers)
end

#! explicit-imports: on

"Return the built-in circular-pipe formula identifiers."
formulas() = FORMULAS

end
