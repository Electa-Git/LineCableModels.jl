"""
    LineCableModels.Engine.InsulationAdmittance

Compute the admittivity of insulator layers.

# Formulas

## `Formula{:lossless}`

Admittivity of a dielectric without conduction or polarization loss, the frequency-domain
constitutive relation of its displacement current:

```math
\\kappa = \\mathrm{j} 2 \\pi f \\varepsilon_0 \\varepsilon_r,
```

where ``f`` is the frequency \\[Hz\\], ``\\kappa`` is the admittivity \\[S/m\\] and
``\\varepsilon_r`` is the relative permittivity `eps_r` of the material \\[dimensionless\\].
The constant ``\\varepsilon_0`` is [`ε₀`](@ref).

### Callable

`formula(material, f)` returns ``\\kappa`` \\[S/m\\] of an insulator layer.

- `material`: material of the layer.
- `f`: frequency \\[Hz\\].

## `Formula{:lossy}`

Admittivity of a dielectric with conduction, polarization loss and displacement current:

```math
\\kappa = \\frac{1}{\\rho} + 2 \\pi f \\varepsilon_0 \\varepsilon_r \\tan\\delta
    + \\mathrm{j} 2 \\pi f \\varepsilon_0 \\varepsilon_r,
```

where ``\\rho`` is the resistivity `rho` of the material \\[Ω·m\\] and ``\\tan\\delta`` is its
loss tangent `tan_delta` \\[dimensionless\\]. The other symbols are as in
`Formula{:lossless}`.

### Callable

`formula(material, f)` returns ``\\kappa`` \\[S/m\\] of an insulator layer.

- `material`: material of the layer.
- `f`: frequency \\[Hz\\].
"""
module InsulationAdmittance

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults
using ...Commons: FormulationOptions
using ..Engine: AbstractInsulationAdmittance
using ..Engine: ε₀
using ...Materials: Material

export Formula, FORMULAS

"Identifiers of the insulation-admittance formulas."
const FORMULAS = (:lossless, :lossy)

"""
$(TYPEDEF)

Represent an admittivity law for an insulator layer.
Its options form a record of type `FormulationOptions{Formula{ID}}`, with `ID` the formula
identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, O <: FormulationOptions} <: AbstractInsulationAdmittance
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
        "unknown insulation-admittance formula :$ID; the formulas are $FORMULAS"))
end

include("lossless.jl")
include("lossy.jl")

end
