"""
$(TYPEDEF)

Store materials by name.

$(TYPEDFIELDS)
"""
mutable struct MaterialsLibrary <: AbstractDict{String, Material}
    "Materials indexed by name."
    data::Dict{String, Material}
end

"""
$(TYPEDSIGNATURES)

Construct a material library.

# Keywords

- `add_defaults`: Add the package's built-in material records. Default: `true`.

# Returns

- A [`MaterialsLibrary`](@ref), populated with built-in records when
  `add_defaults` is `true`.

# Examples

```jldoctest
library = $(FUNCTIONNAME)(; add_defaults=false)
isempty(library)
# output
true
```

"""
function MaterialsLibrary(; add_defaults::Bool = true)::MaterialsLibrary
    library = MaterialsLibrary(Dict{String, Material}())

    if add_defaults
        _add_default_materials!(library)
    end

    return library
end

"""
$(TYPEDSIGNATURES)

Add the built-in material records to `library`.

# Arguments

- `library`: Destination material library.

# Returns

- The modified `library`.

"""
function _add_default_materials!(library::MaterialsLibrary)
    add!(library, "air", Material(:insulator, Inf, 1.0, 1.0, 20.0, 0.0))
    add!(library, "pec", Material(:conductor, eps(), 1.0, 1.0, 20.0, 0.0))
    add!(
        library,
        "copper",
        Material(:conductor, 1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
    )
    add!(
        library,
        "aluminum",
        Material(:conductor, 2.8264e-8, 1.0, 1.000022, 20.0, 0.00429)
    )
    add!(library, "xlpe", Material(:insulator, 1.97e14, 2.5, 1.0, 20.0, 0.0))
    add!(library, "pe", Material(:insulator, 1.97e14, 2.3, 1.0, 20.0, 0.0))
    add!(
        library,
        "semicon1",
        Material(:semicon, 1000.0, 1000.0, 1.0, 20.0, 0.0)
    )
    add!(
        library,
        "semicon2",
        Material(:semicon, 500.0, 1000.0, 1.0, 20.0, 0.0)
    )
    add!(
        library,
        "polyacrylate",
        Material(:semicon, 5.3e3, 32.3, 1.0, 20.0, 0.0)
    )
    add!(library, "lead", Material(:conductor, 21.4e-8, 1.0, 0.999983, 20.0, 0.00400)) # Lead or lead alloy
    add!(library, "steel", Material(:conductor, 13.8e-8, 1.0, 300.0, 20.0, 0.00450)) # Steel
    add!(library, "pp", Material(:insulator, 1e15, 2.8, 1.0, 20.0, 0.0)) # Laminated paper propylene
end

"""
$(TYPEDSIGNATURES)

Add `material` under `name`.

# Arguments

- `library`: Destination material library.
- `name`: Name of the material.
- `material`: Validated material record.

# Returns

- The modified `library`.

# Errors

- Throws `ArgumentError` when `name` already exists.
"""
function add!(
        library::MaterialsLibrary,
        name::Union{AbstractString, Symbol},
        material::Material
)
    key = String(name)
    haskey(library, key) && throw(ArgumentError("material '$key' already exists"))
    candidate = validate(material)
    library[key] = candidate
    return library
end
