
# Implement the AbstractDict interface
Base.length(lib::CablesLibrary) = length(lib.data)
function Base.setindex!(lib::CablesLibrary, value::CableDesign, key)
    cable_id = convert(String, key)
    cable_id == value.cable_id || throw(ArgumentError(
        "cable key '$cable_id' differs from cable_id '$(value.cable_id)'",
    ))
    candidate = validate(value)
    record = DatasheetInfo(candidate.nominal_data)
    lib.data[cable_id] = candidate
    lib.catalogues[cable_id] = record
    return lib
end
Base.iterate(lib::CablesLibrary, state...) = iterate(lib.data, state...)
Base.keys(lib::CablesLibrary) = keys(lib.data)
Base.values(lib::CablesLibrary) = values(lib.data)
Base.haskey(lib::CablesLibrary, key) = haskey(lib.data, key)
Base.getindex(lib::CablesLibrary, key) = getindex(lib.data, key)

"""
$(TYPEDSIGNATURES)

Return the cable design stored under `cable_id`, or `default` when absent.

# Arguments

- `library`: Cable-design library.
- `cable_id`: Stored cable identifier.
- `default`: Value returned when `cable_id` is absent.

# Returns

- The stored [`CableDesign`](@ref), or `default`.

"""
function Base.get(library::CablesLibrary, cable_id, default)
    return get(library.data, cable_id, default)
end

function Base.get(
        default::Union{Function, Type}, library::CablesLibrary, cable_id
)
    return get(default, library.data, cable_id)
end

"""
$(TYPEDSIGNATURES)

Remove the cable design stored under `cable_id`.

# Arguments

- `library`: Cable-design library.
- `cable_id`: Stored cable identifier.

# Returns

- The modified `library`.

"""
function Base.delete!(library::CablesLibrary, cable_id)
    delete!(library.data, cable_id)
    delete!(library.catalogues, cable_id)
    return library
end

"Return an empty cable library."
Base.empty(::CablesLibrary) = CablesLibrary()

"Remove every cable design and catalogue record and return `library`."
function Base.empty!(library::CablesLibrary)
    empty!(library.data)
    empty!(library.catalogues)
    return library
end

"Return a shallow cable-library copy with independent dictionary storage."
function Base.copy(library::CablesLibrary)
    copied = CablesLibrary()
    copied.data = copy(library.data)
    copied.catalogues = copy(library.catalogues)
    return copied
end
