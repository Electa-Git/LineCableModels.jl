"""
    LineCableModels.Engine.EarthAdmittance

Define earth-return admittance recipes, numerical primitives, and
formula-owned frequency functors.

# Dependencies

$(IMPORTS)

"""
module EarthAdmittance

# Export public API
export Formula, formula_id, routes, assumptions, propagation, formulas, Γ

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ...LineCableModels: nominal
import ..Engine: EarthAdmittanceFormulation, formula_id, bessel_difference
#! explicit-imports: off
import ..Engine: description, conductivity, media, special_besselk
#! explicit-imports: on
using QuadGK: quadgk

vacuum_permittivity(value) = one(value) * 88541878128 * (one(value) * 10)^(-22)
vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

"Registered earth-admittance formula selected by `:default`."
const DEFAULT = :Papadopoulos2010

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
            "earth-admittance formula file $(basename(file)) must return its Symbol identifier"
        )
        identifier in identifiers && error(
            "duplicate earth-admittance formula identifier :$identifier"
        )
        push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"Return the built-in earth-admittance formula identifiers."
formulas() = FORMULAS

end # module EarthAdmittance
