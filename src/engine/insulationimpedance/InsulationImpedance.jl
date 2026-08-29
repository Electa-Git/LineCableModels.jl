"""
    LineCableModels.Engine.InsulationImpedance

Define registered series-impedance formulas for cable insulation.

# Dependencies

$(IMPORTS)

"""
module InsulationImpedance

# Export public API
export Formula, formula_id, assumptions, formulas

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: InsulationImpedanceFormulation, formula_id
#! explicit-imports: off
import ..Engine: description
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
            "insulation-impedance formula file $(basename(file)) must return its Symbol identifier"
        )
        identifier in identifiers && error(
            "duplicate insulation-impedance formula identifier :$identifier"
        )
        push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"Return the built-in insulation-impedance formula identifiers."
formulas() = FORMULAS

end # module InsulationImpedance
