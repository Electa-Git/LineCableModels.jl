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

Removes a material from a [`MaterialsLibrary`](@ref).

# Arguments

- `library`: Instance of [`MaterialsLibrary`](@ref) from which the material will be removed.
- `name`: Name of the material to be removed.

# Returns

- The modified instance of [`MaterialsLibrary`](@ref) without the specified material.

# Errors

Throws an error if the material does not exist in the library.

# Examples

```julia
library = MaterialsLibrary()
$(FUNCTIONNAME)(library, "copper")
```

"""
function Base.delete!(library::MaterialsLibrary, name::String)
    haskey(library, name) || throw(KeyError(name))
    delete!(library.data, name)
    return library
end

"""
$(TYPEDSIGNATURES)

Retrieves a material from a [`MaterialsLibrary`](@ref) by name.

# Arguments

- `library`: Instance of [`MaterialsLibrary`](@ref) containing the materials.
- `name`: Name of the material to retrieve.

# Returns

- The requested [`Material`](@ref) if found, otherwise `nothing`.

# Examples

```julia
library = MaterialsLibrary()
material = $(FUNCTIONNAME)(library, "copper")
```

"""
function Base.get(library::MaterialsLibrary, name::String, default = nothing)
    return get(library.data, name, default)
end

"""
$(TYPEDSIGNATURES)

Defines the display representation of a [`Material`](@ref) object for REPL or text output.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `material`: The [`Material`](@ref) instance to be displayed.

# Returns

- Nothing. Modifies `io` by writing text representation of the material.
"""
function Base.show(io::IO, ::MIME"text/plain", material::Material)
    print(io, "Material with properties: [")

    # Define fields to display
    fields = [:rho, :eps_r, :mu_r, :T0, :alpha]

    # Print each field with proper formatting
    for (i, field) in enumerate(fields)
        value = getproperty(material, field)
        # Add comma only between items, not after the last one
        delimiter = i < length(fields) ? ", " : ""
        print(io, "$field=$(round(value, sigdigits=4))$delimiter")
    end

    print(io, "]")
end

"""
$(TYPEDSIGNATURES)

Defines the display representation of a [`MaterialsLibrary`](@ref) object for REPL or text output.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `library`: The [`MaterialsLibrary`](@ref) instance to be displayed.

# Returns

- Nothing. Modifies `io` by writing text representation of the library.
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

        # If there are more materials than we're showing
        if num_materials > shown_materials
            print(io, "\n└─ ... and $(num_materials - shown_materials) more")
        end
    end
end

"""
$(TYPEDSIGNATURES)

Defines the display representation of a [`MaterialsLibrary`](@ref) object for REPL or text output.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `dict`: The [`MaterialsLibrary`](@ref) contents to be displayed.

# Returns

- Nothing. Modifies `io` by writing text representation of the library.
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

        # If there are more materials than we're showing
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
