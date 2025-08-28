"""
$(TYPEDEF)

Represents a library of cable designs stored as a dictionary.

$(TYPEDFIELDS)
"""
mutable struct CablesLibrary
    "Dictionary mapping cable IDs to the respective CableDesign objects."
    data::Dict{String,CableDesign}

    @doc """
    $(TYPEDSIGNATURES)

    Constructs an empty [`CablesLibrary`](@ref) instance.

    # Arguments

    - None.

    # Returns

    - A [`CablesLibrary`](@ref) object with an empty dictionary of cable designs.

    # Examples

    ```julia
    # Create a new, empty library
    library = $(FUNCTIONNAME)()
    ```

    # See also

    - [`CableDesign`](@ref)
    - [`add!`](@ref)
    - [`delete!`](@ref)
    - [`LineCableModels.ImportExport.save`](@ref)
    - [`DataFrame`](@ref)
    """
    function CablesLibrary()::CablesLibrary
        library = new(Dict{String,CableDesign}())
        @info "Initializing empty cables database..."
        return library
    end
end

import ..LineCableModels: add!

"""
Stores a cable design in a [`CablesLibrary`](@ref) object.

# Arguments

- `library`: An instance of [`CablesLibrary`](@ref) to which the cable design will be added.
- `design`: A [`CableDesign`](@ref) object representing the cable design to be added. This object must have a `cable_id` field to uniquely identify it.

# Returns

- None. Modifies the `data` field of the [`CablesLibrary`](@ref) object in-place by adding the new cable design.

# Examples
```julia
library = CablesLibrary()
design = CableDesign("example", ...) # Initialize CableDesign with required fields
add!(library, design)
println(library) # Prints the updated dictionary containing the new cable design
```
# See also

- [`CablesLibrary`](@ref)
- [`CableDesign`](@ref)
- [`delete!`](@ref)
"""
function add!(library::CablesLibrary, design::CableDesign)
    library.data[design.cable_id] = design
    @info "Cable design with ID `$(design.cable_id)` added to the library."
    library
end

include("cableslibrary/base.jl")
include("cableslibrary/dataframe.jl")


