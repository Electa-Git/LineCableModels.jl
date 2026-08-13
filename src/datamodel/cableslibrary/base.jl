
# Implement the AbstractDict interface
Base.length(lib::CablesLibrary) = length(lib.data)
Base.setindex!(lib::CablesLibrary, value::CableDesign, key::String) = (lib.data[key] = value)
Base.iterate(lib::CablesLibrary, state...) = iterate(lib.data, state...)
Base.keys(lib::CablesLibrary) = keys(lib.data)
Base.values(lib::CablesLibrary) = values(lib.data)
Base.haskey(lib::CablesLibrary, key::String) = haskey(lib.data, key)
Base.getindex(lib::CablesLibrary, key::String) = getindex(lib.data, key)

"""
$(TYPEDSIGNATURES)

Retrieves a cable design from a [`CablesLibrary`](@ref) object by its ID.

# Arguments

- `library`: An instance of [`CablesLibrary`](@ref) from which the cable design will be retrieved.
- `cable_id`: The ID of the cable design to retrieve.

# Returns

- A [`CableDesign`](@ref) object corresponding to the given `cable_id` if found, otherwise `nothing`.

# Examples

```julia
library = CablesLibrary()
design = CableDesign("example", ...) # Initialize a CableDesign
add!(library, design)

# Retrieve the cable design
retrieved_design = $(FUNCTIONNAME)(library, "cable1")
println(retrieved_design.id)  # Prints "example"

# Attempt to retrieve a non-existent design
missing_design = $(FUNCTIONNAME)(library, "nonexistent_id")
println(missing_design === nothing)  # Prints true
```

# See also

- [`CablesLibrary`](@ref)
- [`CableDesign`](@ref)
- [`add!`](@ref)
- [`delete!`](@ref)
"""
function Base.get(library::CablesLibrary, cable_id::String, default=nothing)
    if haskey(library, cable_id)
        @info "Cable design with ID `$cable_id` loaded from the library."
        return library[cable_id]
    else
        @warn "Cable design with ID `$cable_id` not found in the library; returning default."
        return default
    end
end

"""
$(TYPEDSIGNATURES)

Removes a cable design from a [`CablesLibrary`](@ref) object by its ID.

# Arguments

- `library`: An instance of [`CablesLibrary`](@ref) from which the cable design will be removed.
- `cable_id`: The ID of the cable design to remove.

# Returns

- Nothing. Modifies the `data` field of the [`CablesLibrary`](@ref) object in-place by removing the specified cable design if it exists.

# Examples

```julia
library = CablesLibrary()
design = CableDesign("example", ...) # Initialize a CableDesign
add!(library, design)

# Remove the cable design
$(FUNCTIONNAME)(library, "example")
haskey(library, "example")  # Returns false
```

# See also

- [`CablesLibrary`](@ref)
- [`add!`](@ref)
"""
function Base.delete!(library::CablesLibrary, cable_id::String)
    if haskey(library, cable_id)
        delete!(library.data, cable_id)
        @info "Cable design with ID `$cable_id` removed from the library."
    else
        @error "Cable design with ID `$cable_id` not found in the library; cannot delete."
        throw(KeyError(cable_id))
    end
end