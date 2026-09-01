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
import ..Engine: EarthImpedanceFormulation, formula_id, bessel_difference
#! explicit-imports: off
import ..Engine: description, conductivity, formula_method, media, special_besselk
using SpecialFunctions: bessely, hankelh1
#! explicit-imports: on
using QuadGK: quadgk

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
    discovered = Set{Symbol}()
    for file in files
        identifier = include(file)
        identifier isa Symbol || error(
            "earth-impedance formula file $(basename(file)) must return its Symbol identifier"
        )
        identifier in discovered && error(
            "duplicate earth-impedance formula identifier :$identifier"
        )
        push!(discovered, identifier)
        propagation(Val(identifier)) === Val(:backend) ||
            push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"Return the built-in earth-impedance formula identifiers."
formulas() = FORMULAS

end # module EarthImpedance
