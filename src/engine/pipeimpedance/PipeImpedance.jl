"""
    LineCableModels.Engine.PipeImpedance

Own pipe-type formula selections and backend applicability. No analytical
pipe-type implementation is supplied yet.
"""
module PipeImpedance
import ...LineCableModels: FormulaDefinition

export Formula, formula_id, formulas

using DocStringExtensions: TYPEDEF, TYPEDSIGNATURES
import ..Engine: PipeImpedanceFormulation, Formulation
import ...LineCableModels: formula_id
#! explicit-imports: off
# These bindings are consumed by the dynamically included formula definitions.
import ..Engine: LineCableModelsCoaxial, LineCableModelsFEM
import ...LineCableModels: description
#! explicit-imports: on
import ...DataModel
using ...DataModel: CableDesign

include("interface.jl")

#! explicit-imports: off
# Formula files are discovered and included dynamically, as in the other owners.
const FORMULAS = let
    directory = joinpath(@__DIR__, "formulas")
    Base.include_dependency(directory)
    identifiers = Symbol[]
    for path in sort!(filter(endswith(".jl"), readdir(directory; join = true)))
        identifier = include(path)
        identifier isa Symbol ||
            error("pipe formula file must return its Symbol identifier: $path")
        identifier in identifiers &&
            error("duplicate pipe-impedance identifier :$identifier")
        push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"Return registered pipe selections, including the explicit default policy."
formulas() = FORMULAS

end
