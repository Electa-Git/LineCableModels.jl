"""
$(TYPEDEF)

Store cable designs by `cable_id`.

$(TYPEDFIELDS)
"""
mutable struct CablesLibrary
    "Cable designs indexed by `cable_id`."
    data::Dict{String, CableDesign}

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
        return new(Dict{String, CableDesign}())
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
function add!(library::CablesLibrary, design::CableDesign)
    haskey(library, design.cable_id) && throw(ArgumentError(
        "cable design '$(design.cable_id)' already exists",
    ))
    library.data[design.cable_id] = validate(design)
    return library
end

include("base.jl")
include("dataframe.jl")
include("vdeparse.jl")
