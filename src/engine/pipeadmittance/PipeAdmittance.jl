"""
    LineCableModels.Engine.PipeAdmittance

Evaluate pipe-interior potential coefficients before complete matrix assembly
and inversion. The individual scalar values are not nodal admittances.

# Dependencies

$(IMPORTS)
"""
module PipeAdmittance

export Formula, formula_id, formulas, routes, assumptions
#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: PipeAdmittanceFormulation, formula_id, PipeImpedance
#! explicit-imports: off
import ...LineCableModels: FormulaMethod
import ..Engine: description
#! explicit-imports: on

const DEFAULT=:Kane1995
include("interface.jl")

#! explicit-imports: off
const FORMULAS=let
    directory=joinpath(@__DIR__,"formulas")
    Base.include_dependency(directory)
    identifiers=Symbol[]
    for path in sort(filter(endswith(".jl"),readdir(directory;join=true)))
        id=include(path)
        id isa Symbol || error("pipe potential formula must return its Symbol identifier")
        id in identifiers && error("duplicate pipe potential formula :$id")
        push!(identifiers,id)
    end
    Tuple(identifiers)
end

#! explicit-imports: on

"Return the registered pipe-interior potential formula identifiers."
formulas()=FORMULAS

end
