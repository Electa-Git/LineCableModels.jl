"""
    LineCableModels.Engine.InternalImpedance

Define conductor internal-impedance recipes and their electromagnetic
interaction formulas.

# Dependencies

$(IMPORTS)

"""
module InternalImpedance

# Export public API
export Formula, formula_id, routes, assumptions, formulas

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: InternalImpedanceFormulation, formula_id
#! explicit-imports: off
import ...LineCableModels: FormulaMethod
import ..Engine: description, conductivity
import ..Engine: special_besselix, special_besselkx
#! explicit-imports: on

vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

"Registered internal-impedance formula selected by `:default`."
const DEFAULT = :Schelkunoff1934

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
            "internal-impedance formula file $(basename(file)) must return its Symbol identifier"
        )
        identifier in identifiers && error(
            "duplicate internal-impedance formula identifier :$identifier"
        )
        push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"Return the built-in internal-impedance formula identifiers."
formulas() = FORMULAS

end # module InternalImpedance
