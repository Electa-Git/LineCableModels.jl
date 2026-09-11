"""
    LineCableModels.Engine.EarthImpedance

Define earth-return impedance recipes, numerical primitives, and formula-owned
frequency functors.

# Dependencies

$(IMPORTS)

"""
module EarthImpedance

# Export public API
export Formula, formula_id, earth_impedance, assumptions, propagation, formulas, Γ

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ...LineCableModels: validate
import ..Engine: EarthPair, hooks, earth_parameters
import ...Earth: EquivalentHomogeneous
import ..Engine: EarthImpedanceFormulation, formula_id
#! explicit-imports: off
# Formula files are discovered dynamically; their imports are verified by the
# equation-ownership tests because the static import scanner cannot follow them.
import ..Engine: system_earth, unified_entry, retained_earth_features
import ...LineCableModels: FormulaDefinition, FormulaMethod, nominal
import ..Engine: Formulation, SpectralIntegral, integrate
import ..Engine: description, conductivity, media, special_besselk
import ..Engine: computation_options, LineCableModelsCoaxial
#! explicit-imports: on

vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("interface.jl")
include("homogeneous.jl")

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
            "earth-impedance formula file $(basename(file)) must return its Symbol identifier"
        )
        identifier in identifiers && error(
            "duplicate earth-impedance formula identifier :$identifier"
        )
        push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"Return numerical earth-impedance identifiers."
formulas() = FORMULAS

end # module EarthImpedance
