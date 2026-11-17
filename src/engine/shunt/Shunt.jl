"""
    LineCableModels.Engine.Shunt

Compute the potential coefficients inside a cable from the admittivities of its layers.

# Formulas

## `Formula{:series}`

Potential coefficients inside a cable from the series combination of its annular dielectric
layers [Ametani1980, Ametani2015b](@cite). A layer ``l`` between the radii ``a_l`` and
``b_l`` has the coefficient ``p_l``, and the layers ``\\mathcal{L}_i`` in the interval of
conductor ``i`` add to its own coefficient ``p_i``. Each concentric assembly numbers its
conductors ``1, \\ldots, n`` from the inside out. The entry between conductors
``i`` and ``j`` sums the own coefficients from the outer of the two to the outermost one:

```math
p_l = \\frac{s \\ln(b_l / a_l)}{2 \\pi \\kappa_l}, \\qquad
p_i = \\sum_{l \\in \\mathcal{L}_i} p_l, \\qquad
P_{\\mathrm{in}, ij} = \\sum_{k = \\max(i, j)}^{n} p_k,
```

where ``s`` is the complex frequency \\[1/s\\], and ``a_l`` and ``b_l`` are in \\[m\\].
The layer admittivity ``\\kappa_l`` is in \\[S/m\\], and the coefficients are in \\[m/F\\].
The calculation returns zero at ``a_l = 0`` and evaluates the expression otherwise.
At ``a_l = b_l``, the expression has value zero. Conductors of different cables have no
entry in ``\\mathbf P_\\mathrm{in}``.

An assembly whose bore holds members adds the coefficients of the gap between them. The
entries between the conductors of the members ``m`` and ``m'``, each with its members' conductors,
add

```math
p_{mm'} = \\frac{s \\ln(c / d_{mm'})}{2 \\pi \\kappa_f} + \\sum_l p_l + P_\\beta
    + \\frac{s}{2 \\pi \\kappa_f} A_{mm'},
```

where ``c`` is the outer radius of the filler, the region around the members, and
``\\kappa_f`` its admittivity. The radius ``d_{mm'}`` is the distance between the axes of
two members, and ``d_{mm}`` is the outer radius of the interval of ``m``. The sum runs over
the annuli between the filler and the bore, and ``P_\\beta`` is the sum of the own
coefficients from the bore outward. The image matrix ``A`` holds the fields that the members
and the bore return above order 0. The formulation solves it as for the impedances, with
each conductor an equipotential.

### Callable

`formula(r_in, r_ex, κ, s)` returns ``p_l`` \\[m/F\\].

- `r_in`: inner radius of the layer \\[m\\].
- `r_ex`: outer radius of the layer \\[m\\].
- `κ`: admittivity of the layer \\[S/m\\].
- `s`: complex frequency, `j 2π f` \\[1/s\\].
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
