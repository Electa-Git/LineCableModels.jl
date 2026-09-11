"""
    LineCableModels.Materials.TemperatureDependent

Evaluate electrical resistivity at a prescribed temperature from a material's
reference calibration and a selected constitutive equation.

# Dependencies

$(IMPORTS)
"""
module TemperatureDependent

export Formula, formula_id, formulas
public temperature_resistivity

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
using ..Materials: Material
import ...TextDisplay
import ...Grammar: AbstractFormulation, computation_options
import ...LineCableModels: FormulaDefinition, FormulaMethod, constitutive,
                          formula_id, validate
#! explicit-imports: off
# Used by the formula files discovered and included below.
import ...LineCableModels: description
#! explicit-imports: on

include("interface.jl")

#! explicit-imports: off
const FORMULAS = let
    directory = joinpath(@__DIR__, "formulas")
    Base.include_dependency(directory)
    files = sort!(filter(path -> endswith(path, ".jl"), readdir(directory; join=true)))
    identifiers = Symbol[]
    for file in files
        identifier = include(file)
        identifier isa Symbol || error(
            "temperature formula file $(basename(file)) must return its Symbol identifier")
        identifier in identifiers && error("duplicate temperature formula :$identifier")
        push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"""Return the registered electrical-resistivity temperature laws."""
formulas() = FORMULAS

end # module TemperatureDependent
