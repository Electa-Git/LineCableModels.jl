"""
    LineCableModels.Engine.EarthProperties

Compute the resistivity, relative permittivity, and relative permeability of earth and air,
and the equivalent earth of the layers below the air.

# Formulas

## `Formula{:constant}`

Properties of a layer of the earth or of the air independent of the frequency ``f``:

```math
\\rho(f) = \\rho, \\qquad \\varepsilon_r(f) = \\varepsilon_r, \\qquad \\mu_r(f) = \\mu_r,
```

where ``\\rho`` is the resistivity of the layer \\[Ω·m\\], and ``\\varepsilon_r`` and
``\\mu_r`` are its relative permittivity and permeability \\[dimensionless\\].

### Callable

`formula(layer, f)` returns `layer`, whose properties hold at every frequency.

- `layer`: layer of the earth or of the air.
- `f`: frequency \\[Hz\\].

## `Formula{:bottommost}`

Equivalent earth of the layers below the air: the deepest described layer, the identity for
one layer.

### Callable

`formula(layers, f)` returns the last of `layers`.

- `layers`: the layers below the air, as the law describes them at `f`.
- `f`: frequency \\[Hz\\].
"""
module EarthProperties

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults, validate
using ...Commons: FormulationOptions
using ..Engine: AbstractEarthProperties
using ...Earth: EarthLayer, EarthModel

export Formula, FORMULAS

"Identifiers of the earth-properties formulas."
const FORMULAS = (:constant, :bottommost)

"""
$(TYPEDEF)

Represent a law of the properties of earth and air, or a model of the equivalent earth.
Its options form a record of type `FormulationOptions{Formula{ID}}`, with `ID` the formula
identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, O <: FormulationOptions} <: AbstractEarthProperties
    "Options of the formula."
    options::O
    function Formula{ID}(options::O) where {ID, O <: FormulationOptions{Formula{ID}}}
        return new{ID, O}(options)
    end
end

"""
$(TYPEDSIGNATURES)

Build the formula named by `ID` with its declared defaults and supplied options.

# Arguments

- `Val{ID}`: formula identifier from [`FORMULAS`](@ref).

# Keywords

- `kw...`: options declared by [`defaults`](@ref) for `Formula{ID}`.

# Returns

- A [`Formula`](@ref) with a typed option record.

# Errors

- Throws `ArgumentError` for an unknown identifier or option.
"""
Formula(::Val{ID}; kw...) where {ID} = Formula{ID}(FormulationOptions(Formula{ID}; kw...))

# An identifier without a formula has no defaults.
function defaults(::Type{Formula{ID}}) where {ID}
    throw(ArgumentError(
        "unknown earth-properties formula :$ID; the formulas are $FORMULAS"))
end

include("constant.jl")
include("bottommost.jl")

# A law describes a layer, and a model reduces the described layers of an earth model to one.
validate(formula::Formula{:constant}, ::Type{EarthLayer}) = formula
validate(formula::Formula{:bottommost}, ::Type{EarthModel}) = formula
function validate(formula::Formula{ID}, ::Type{EarthLayer}) where {ID}
    throw(DomainError(formula,
        "the earth-properties formula :$ID describes no layer; a model of the equivalent earth belongs in the slot equivalent_earth"))
end
function validate(formula::Formula{ID}, ::Type{EarthModel}) where {ID}
    throw(DomainError(formula,
        "the earth-properties formula :$ID reduces no earth model; a law belongs in the slot earth_properties"))
end

end
