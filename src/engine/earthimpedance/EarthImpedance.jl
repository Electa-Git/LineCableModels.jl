"""
    LineCableModels.Engine.EarthImpedance

Define earth-return impedance recipes, numerical primitives, and formula-owned
frequency functors.

# Dependencies

$(IMPORTS)

"""
module EarthImpedance

# Export public API
export Formula, formula_id, routes, assumptions, propagation, formulas, Γ

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ...LineCableModels: nominal
import ...LineCableModels: validate
import ..Engine: EarthPair
import ..Engine: EarthImpedanceFormulation, formula_id
#! explicit-imports: off
import ...LineCableModels: FormulaMethod
import ..Engine: description, conductivity, media, special_besselk
using SpecialFunctions: bessely, hankelh1
#! explicit-imports: on
using QuadGK: quadgk

vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("interface.jl")
include("homogeneous.jl")

#! explicit-imports: off
const REGISTERED,
FORMULAS = let
    directory = joinpath(@__DIR__, "formulas")
    Base.include_dependency(directory)
    files = sort!(filter(
        path -> endswith(path, ".jl"),
        readdir(directory; join = true)
    ))
    identifiers = Symbol[]
    discovered = Symbol[]
    for file in files
        identifier = include(file)
        identifier isa Symbol || error(
            "earth-impedance formula file $(basename(file)) must return its Symbol identifier"
        )
        identifier in discovered && error(
            "duplicate earth-impedance formula identifier :$identifier"
        )
        push!(discovered, identifier)
        (identifier === :default || propagation(Val(identifier)) === Val(:backend)) ||
            push!(identifiers, identifier)
    end
    Tuple(discovered), Tuple(identifiers)
end
#! explicit-imports: on

"Return numerical earth-impedance identifiers; `:default` is a context selector."
formulas() = FORMULAS

end # module EarthImpedance
