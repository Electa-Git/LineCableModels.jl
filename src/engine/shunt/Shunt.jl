"""
    LineCableModels.Engine.Shunt

Compute the potential coefficients inside a cable from the admittivities of its layers.

# Formulas

## `Formula{:series}`

Potential coefficients inside a cable from the series combination of its annular dielectric
layers [Ametani1980, Ametani2015b](@cite). A layer ``l`` between the radii ``a_l`` and
``b_l`` has the coefficient ``p_l``, and the layers ``\\mathcal{L}_i`` in the interval of
conductor ``i`` add to its own coefficient ``p_i``. In each concentric assembly, with the
conductors numbered ``1, \\ldots, n`` from the inside out, the entry between conductors
``i`` and ``j`` sums the own coefficients of the conductors from the outer of the two to the
outermost one:

```math
p_l = \\frac{s \\ln(b_l / a_l)}{2 \\pi \\kappa_l}, \\qquad
p_i = \\sum_{l \\in \\mathcal{L}_i} p_l, \\qquad
P_{\\mathrm{in}, ij} = \\sum_{k = \\max(i, j)}^{n} p_k,
```

where ``s`` is the complex frequency \\[1/s\\], ``a_l`` and ``b_l`` are in \\[m\\],
``\\kappa_l`` is the admittivity of the layer \\[S/m\\], and the coefficients are in
\\[m/F\\]. A layer with ``a_l = 0`` or ``a_l = b_l`` has ``p_l = 0``. Conductors of
different assemblies have no entry in ``\\mathbf P_\\mathrm{in}``.
"""
module Shunt

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults
using ...Commons: FormulationOptions
using ..Engine: AbstractShunt
import ..Engine: potential!
using ..Engine: CoaxialInput

export Formula, FORMULAS

"Identifiers of the shunt formulas."
const FORMULAS = (:series,)

"""
$(TYPEDEF)

Represent a formula for the potential coefficients inside a cable.
Its options form a record of type `FormulationOptions{Formula{ID}}`, with `ID` the formula
identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, O <: FormulationOptions} <: AbstractShunt
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
        "unknown shunt formula :$ID; the formulas are $FORMULAS"))
end

include("series.jl")

end
