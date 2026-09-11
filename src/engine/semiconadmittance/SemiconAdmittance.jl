"""
    LineCableModels.Engine.SemiconAdmittance

Define registered constitutive relations for semiconducting-screen admittance.

# Dependencies

$(IMPORTS)

"""
module SemiconAdmittance
import ...Grammar: computation_options

export Formula, formula_id, formulas

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: SemiconAdmittanceFormulation, formula_id
import ...LineCableModels: FormulaDefinition, FormulaMethod
using ...Materials: Material
#! explicit-imports: off
import ..Engine: description, conductivity
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
            "semicon-admittance formula file $(basename(file)) must return its Symbol identifier"
        )
        identifier in identifiers && error(
            "duplicate semicon-admittance formula identifier :$identifier"
        )
        push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"Return the built-in semicon-admittance formula identifiers."
formulas() = FORMULAS

end # module SemiconAdmittance
