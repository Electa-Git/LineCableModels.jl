Base.eltype(::Material{T}) where {T} = T
Base.eltype(::Type{Material{T}}) where {T} = T

# Implement the AbstractDict interface
Base.length(lib::MaterialsLibrary) = length(lib.data)
function Base.setindex!(lib::MaterialsLibrary, value::Material, key::String)
    (lib.data[key] = value)
end
Base.iterate(lib::MaterialsLibrary, state...) = iterate(lib.data, state...)
Base.keys(lib::MaterialsLibrary) = keys(lib.data)
Base.values(lib::MaterialsLibrary) = values(lib.data)
Base.haskey(lib::MaterialsLibrary, key::String) = haskey(lib.data, key)
Base.getindex(lib::MaterialsLibrary, key::String) = getindex(lib.data, key)

"""
$(TYPEDSIGNATURES)

Remove the material stored under `name`.

# Arguments

- `library`: Material library.
- `name`: Stored material name.

# Returns

- The modified `library`.

# Errors

- Throws `KeyError` when `name` is absent.

"""
function Base.delete!(library::MaterialsLibrary, name::String)
    haskey(library, name) || throw(KeyError(name))
    delete!(library.data, name)
    return library
end

"""
$(TYPEDSIGNATURES)

Return the material stored under `name`, or `default` when absent.

# Arguments

- `library`: Material library.
- `name`: Stored material name.
- `default`: Value returned when `name` is absent. Default: `nothing`.

# Returns

- The stored [`Material`](@ref), or `default`.

"""
function Base.get(library::MaterialsLibrary, name::String, default = nothing)
    return get(library.data, name, default)
end

"""
$(TYPEDSIGNATURES)

Write the plain-text representation of a material to `io`.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `material`: Material to display.

# Returns

- `nothing` after writing to `io`.
"""
function Base.show(io::IO, ::MIME"text/plain", material::Material)
    print(io, "Material with properties: [")

    # Select displayed fields.
    fields = [:rho, :eps_r, :mu_r, :T0, :alpha]

    # Print each field with four significant digits.
    for (i, field) in enumerate(fields)
        value = getproperty(material, field)
        # Separate adjacent fields with a comma.
        delimiter = i < length(fields) ? ", " : ""
        print(io, "$field=$(round(value, sigdigits=4))$delimiter")
    end

    print(io, "]")
end

"""
$(TYPEDSIGNATURES)

Write a summary of a material library to `io`.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `library`: Material library to display.

# Returns

- `nothing` after writing to `io`.
"""
function Base.show(io::IO, ::MIME"text/plain", library::MaterialsLibrary)
    num_materials = length(library)
    material_word = num_materials == 1 ? "material" : "materials"
    print(io, "MaterialsLibrary with $num_materials $material_word")

    if num_materials > 0
        print(io, ":")
        # Optional: list the first few materials
        shown_materials = min(5, num_materials)
        material_names = sort!(collect(keys(library)))[1:shown_materials]

        for (i, name) in enumerate(material_names)
            print(io, "\n$(i == shown_materials ? "└─" : "├─") $name")
        end

        # Report entries omitted from the display.
        if num_materials > shown_materials
            print(io, "\n└─ ... and $(num_materials - shown_materials) more")
        end
    end
end

"""
$(TYPEDSIGNATURES)

Write a summary of a material dictionary to `io`.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `dict`: Material dictionary to display.

# Returns

- `nothing` after writing to `io`.
"""
function Base.show(io::IO, ::MIME"text/plain", dict::Dict{String, Material})
    num_materials = length(dict)
    material_word = num_materials == 1 ? "material" : "materials"
    print(io, "Dict{String, Material} with $num_materials $material_word")

    if num_materials > 0
        print(io, ":")
        # List the first few materials
        shown_materials = min(5, num_materials)
        material_names = sort!(collect(keys(dict)))[1:shown_materials]

        for (i, name) in enumerate(material_names)
            print(io, "\n$(i == shown_materials ? "└─" : "├─") $name")
        end

        # Report entries omitted from the display.
        if num_materials > shown_materials
            print(io, "\n└─ ... and $(num_materials - shown_materials) more")
        end
    end
end

function Base.convert(::Type{Material{T}}, m::Material) where {T <: Real}
    Material{T}(convert(T, m.rho), convert(T, m.eps_r), convert(T, m.mu_r),
        convert(T, m.T0), convert(T, m.alpha))
end

Base.convert(::Type{Material{T}}, material::Material{T}) where {T <: Real} = material
