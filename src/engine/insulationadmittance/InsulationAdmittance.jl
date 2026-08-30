"""
    LineCableModels.Engine.InsulationAdmittance

Define registered constitutive relations for cable-insulation admittance.
`:Marti2001` retains parallel conductance and capacitance;
`:Gustavsen2013` selects the conventional lossless approximation.

# Dependencies

$(IMPORTS)

"""
module InsulationAdmittance

# Export public API
export Formula, formula_id, assumptions, formulas

# Module-specific dependencies
#! explicit-imports: off
# IMPORTS is expanded in this module docstring rather than called as Julia code.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: InsulationAdmittanceFormulation, formula_id
using ...Materials: Material
import ...LineCableModels: constitutive
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
            "insulation-admittance formula file $(basename(file)) must return its Symbol identifier"
        )
        identifier in identifiers && error(
            "duplicate insulation-admittance formula identifier :$identifier"
        )
        push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"Return the built-in insulation-admittance formula identifiers."
formulas() = FORMULAS

end # module InsulationAdmittance
