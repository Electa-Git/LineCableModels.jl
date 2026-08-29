"""
    LineCableModels.EarthProps.FD

Define measured and material-physics relations that map one soil material and
one frequency to its frequency-dependent electromagnetic properties.

# Dependencies

$(IMPORTS)
"""
module FD

export Formula, formula_id, assumptions, formulas

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
using ..EarthProps: EarthMaterial
import ...Grammar: AbstractFormulation
import ...LineCableModels: constitutive, formula_id
#! explicit-imports: off
import ...LineCableModels: description
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
            "earth-material formula file $(basename(file)) must return its Symbol identifier"
        )
        identifier in identifiers && error(
            "duplicate earth-material formula identifier :$identifier"
        )
        push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"Return the built-in frequency-dependent earth-material formula identifiers."
formulas() = FORMULAS

end # module FD
