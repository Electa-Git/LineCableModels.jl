
# Implement the AbstractDict interface
Base.length(lib::CablesLibrary) = length(lib.entries)
function Base.setindex!(lib::CablesLibrary, design::CableDesign, key)
    cable_id = convert(String, key)
    cable_id == design.cable_id || throw(ArgumentError(
        "cable key '$cable_id' differs from cable_id '$(design.cable_id)'",
    ))
    validate(design)
    record = DatasheetInfo(design.nominal_data)
    lib.entries[cable_id] = design
    lib.datasheets[cable_id] = record
    return lib
end
Base.iterate(lib::CablesLibrary, state...) = iterate(lib.entries, state...)
Base.keys(lib::CablesLibrary) = keys(lib.entries)
Base.values(lib::CablesLibrary) = values(lib.entries)
Base.haskey(lib::CablesLibrary, key) = haskey(lib.entries, key)
Base.getindex(lib::CablesLibrary, key) = getindex(lib.entries, key)

"""
$(TYPEDSIGNATURES)

Return the cable design stored under `cable_id`, or `default` when absent.

# Arguments

- `library`: cable-design library.
- `cable_id`: stored cable identifier.
- `default`: value returned when `cable_id` is absent.

# Returns

- The stored [`CableDesign`](@ref), or `default`.

"""
function Base.get(library::CablesLibrary, cable_id, default)
    return get(library.entries, cable_id, default)
end

function Base.get(
        default::Union{Function, Type}, library::CablesLibrary, cable_id
)
    return get(default, library.entries, cable_id)
end

"""
$(TYPEDSIGNATURES)

Remove the cable design stored under `cable_id`.

# Arguments

- `library`: cable-design library.
- `cable_id`: stored cable identifier.

# Returns

- The modified `library`.

"""
function Base.delete!(library::CablesLibrary, cable_id)
    delete!(library.entries, cable_id)
    delete!(library.datasheets, cable_id)
    return library
end

"""
Return an empty cable library.
"""
Base.empty(::CablesLibrary) = CablesLibrary()

"""
Remove every cable design and catalog record and return `library`.
"""
function Base.empty!(library::CablesLibrary)
    empty!(library.entries)
    empty!(library.datasheets)
    return library
end

"""
Return a shallow cable-library copy with independent dictionary storage.
"""
function Base.copy(library::CablesLibrary)
    copied = CablesLibrary()
    copied.entries = copy(library.entries)
    copied.datasheets = copy(library.datasheets)
    return copied
end
