"""
$(TYPEDEF)

Store cable designs by `cable_id`.

$(TYPEDFIELDS)
"""
mutable struct CablesLibrary
    "Cable designs indexed by `cable_id`."
    data::Dict{String, CableDesign}
    "Catalogue records indexed by `cable_id`."
    catalogues::Dict{String, NamedTuple}

    @doc """
     $(TYPEDSIGNATURES)

     Construct an empty cable-design library.

     # Returns

     - An empty [`CablesLibrary`](@ref).

     # Examples

     ```jldoctest
     library = $(FUNCTIONNAME)()
     isempty(library)
     # output
     true
     ```

     """
    function CablesLibrary()::CablesLibrary
        return new(Dict{String, CableDesign}(), Dict{String, NamedTuple}())
    end
end

"""
$(TYPEDSIGNATURES)

Add a cable design under its `cable_id`.

# Arguments

- `library`: Destination library.
- `design`: Validated cable design.

# Returns

- The modified `library`.

# Errors

- Throws `ArgumentError` when `library` already contains `design.cable_id`.
"""
function add!(
        library::CablesLibrary,
        design::CableDesign;
        catalogue::NamedTuple = design.nominal_data
)
    haskey(library, design.cable_id) && throw(ArgumentError(
        "cable design '$(design.cable_id)' already exists",
    ))
    library.data[design.cable_id] = validate(design)
    library.catalogues[design.cable_id] = catalogue
    return library
end

"Return the catalogue record associated with `cable_id`."
catalogue(library::CablesLibrary, cable_id::AbstractString) =
    library.catalogues[String(cable_id)]

include("base.jl")
include("dataframe.jl")
include("vdeparse.jl")
