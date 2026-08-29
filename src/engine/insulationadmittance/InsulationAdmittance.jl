"""
    LineCableModels.Engine.InsulationAdmittance

Define registered shunt-admittance formulas for cable insulation. `:Lossless`
retains capacitance only. `:ParallelRC` retains capacitance and dielectric
conductance.

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
#! explicit-imports: off
import ..Engine: conductivity, description
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
