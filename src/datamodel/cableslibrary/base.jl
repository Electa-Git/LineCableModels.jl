
# Implement the AbstractDict interface
Base.length(lib::CablesLibrary) = length(lib.data)
function Base.setindex!(lib::CablesLibrary, value::CableDesign, key::String)
    lib.data[key] = value
    get!(lib.catalogues, key, (;))
    return value
end
Base.iterate(lib::CablesLibrary, state...) = iterate(lib.data, state...)
Base.keys(lib::CablesLibrary) = keys(lib.data)
Base.values(lib::CablesLibrary) = values(lib.data)
Base.haskey(lib::CablesLibrary, key::String) = haskey(lib.data, key)
Base.getindex(lib::CablesLibrary, key::String) = getindex(lib.data, key)

"""
$(TYPEDSIGNATURES)

Return the cable design stored under `cable_id`, or `default` when absent.

# Arguments

- `library`: Cable-design library.
- `cable_id`: Stored cable identifier.
- `default`: Value returned when `cable_id` is absent. Default: `nothing`.

# Returns

- The stored [`CableDesign`](@ref), or `default`.

"""
function Base.get(library::CablesLibrary, cable_id::String, default = nothing)
    return get(library.data, cable_id, default)
end

"""
$(TYPEDSIGNATURES)

Remove the cable design stored under `cable_id`.

# Arguments

- `library`: Cable-design library.
- `cable_id`: Stored cable identifier.

# Returns

- The modified `library`.

# Errors

- Throws `KeyError` when `cable_id` is absent.
"""
function Base.delete!(library::CablesLibrary, cable_id::String)
    haskey(library, cable_id) || throw(KeyError(cable_id))
    delete!(library.data, cable_id)
    delete!(library.catalogues, cable_id)
    return library
end
