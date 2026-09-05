Base.eltype(::Material{T}) where {T} = T
Base.eltype(::Type{Material{T}}) where {T} = T

# Implement the AbstractDict interface
Base.length(lib::MaterialsLibrary) = length(lib.data)
function Base.setindex!(lib::MaterialsLibrary, value::Material, key)
    setindex!(lib.data, validate(value), key)
    return lib
end
Base.iterate(lib::MaterialsLibrary, state...) = iterate(lib.data, state...)
Base.keys(lib::MaterialsLibrary) = keys(lib.data)
Base.values(lib::MaterialsLibrary) = values(lib.data)
Base.haskey(lib::MaterialsLibrary, key) = haskey(lib.data, key)
Base.getindex(lib::MaterialsLibrary, key) = getindex(lib.data, key)

"""
$(TYPEDSIGNATURES)

Remove the material stored under `name`.

# Arguments

- `library`: Material library.
- `name`: Stored material name.

# Returns

- The modified `library`.

"""
function Base.delete!(library::MaterialsLibrary, name)
    delete!(library.data, name)
    return library
end

"""
$(TYPEDSIGNATURES)

Return the material stored under `name`, or `default` when absent.

# Arguments

- `library`: Material library.
- `name`: Stored material name.
- `default`: Value returned when `name` is absent.

# Returns

- The stored [`Material`](@ref), or `default`.

"""
function Base.get(library::MaterialsLibrary, name, default)
    return get(library.data, name, default)
end

function Base.get(
        default::Union{Function, Type}, library::MaterialsLibrary, name
)
    return get(default, library.data, name)
end

"Return an empty material library without repopulating built-in records."
Base.empty(::MaterialsLibrary) = MaterialsLibrary(; add_defaults = false)

"Remove every stored material and return `library`."
function Base.empty!(library::MaterialsLibrary)
    empty!(library.data)
    return library
end

"Return a shallow material-library copy with independent dictionary storage."
Base.copy(library::MaterialsLibrary) = MaterialsLibrary(copy(library.data))

TextDisplay.name(::Type{<:Material}) = "Material"

function Base.summary(io::IO, material::Material)
    print(io, "Material · ", material.kind)
end

function _material_fields(material::Material)
    return (
        ρ = TextDisplay.engineering(material.rho, :ohm_meter),
        εᵣ = TextDisplay.value(material.eps_r),
        μᵣ = TextDisplay.value(material.mu_r),
        T₀ = TextDisplay.engineering(material.T0, :celsius),
        α = iszero(material.alpha) ? nothing :
            TextDisplay.engineering(material.alpha, :kelvin_inverse),
        ρₜₕ = iszero(material.rho_thermal) ? nothing :
              string(TextDisplay.value(material.rho_thermal), " K·m/W"),
        θₘₐₓ = material.theta_max == oftype(material.theta_max, 90) ? nothing :
               TextDisplay.engineering(material.theta_max, :celsius),
        tanδ = iszero(material.tan_delta) ? nothing :
               TextDisplay.value(material.tan_delta),
        σₛ = iszero(material.sigma_solar) ? nothing :
             TextDisplay.value(material.sigma_solar)
    )
end

function Base.show(io::IO, material::Material)
    fields = _material_fields(material)
    print(io, "Material(:", material.kind, "; ρ=", fields.ρ,
        ", εᵣ=", fields.εᵣ, ", μᵣ=", fields.μᵣ)
    for key in (:ρₜₕ, :θₘₐₓ, :tanδ, :σₛ)
        displayed = getproperty(fields, key)
        displayed === nothing && continue
        print(io, ", ", key, "=", displayed)
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", material::Material)
    get(io, :compact, false) && return show(io, material)
    return TextDisplay.fields(
        io,
        "Material · $(material.kind)",
        _material_fields(material);
        multiline = true
    )
end

TextDisplay.name(::Type{<:MaterialsLibrary}) = "MaterialsLibrary"

function Base.summary(io::IO, library::MaterialsLibrary)
    count = length(library)
    print(io, "MaterialsLibrary with ", count, count == 1 ? " material" : " materials")
end

function Base.show(io::IO, library::MaterialsLibrary)
    count = length(library)
    print(io, "MaterialsLibrary(", count, count == 1 ? " material)" : " materials)")
end

function Base.show(io::IO, ::MIME"text/plain", library::MaterialsLibrary)
    get(io, :compact, false) && return show(io, library)
    count = length(library)
    header = "MaterialsLibrary · $count $(count == 1 ? "material" : "materials")"
    children = [
        string(name, " · ", sprint(summary, library[name]))
        for name in sort!(collect(keys(library)))
    ]
    return TextDisplay.tree(io, header, children; noun = "materials")
end

function Base.convert(::Type{Material{T}}, m::Material) where {T <: Real}
    Material{T}(
        m.kind,
        convert(T, m.rho),
        convert(T, m.eps_r),
        convert(T, m.mu_r),
        convert(T, m.T0),
        convert(T, m.alpha),
        convert(T, m.rho_thermal),
        convert(T, m.theta_max),
        convert(T, m.tan_delta),
        convert(T, m.sigma_solar)
    )
end

Base.convert(::Type{Material{T}}, material::Material{T}) where {T <: Real} = material
