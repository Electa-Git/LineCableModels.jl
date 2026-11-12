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
different cables have no entry in ``\\mathbf P_\\mathrm{in}``.

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

At a harmonic order ``n \\geq 1``, a layer is a line of the field that varies as
``\\cos n\\varphi``. With ``p_n = s / (2 \\pi n \\kappa_l)`` and
``\\rho = (a_l / b_l)^{2n}``, its record is

```math
p_\\mathrm{inner} = p_\\mathrm{outer} = p_n \\frac{1 + \\rho}{1 - \\rho}, \\qquad
p_\\mathrm{transfer} = \\frac{2 p_n (a_l / b_l)^n}{1 - \\rho}.
```

### Callable

`formula(r_in, r_ex, κ, s)` returns ``p_l`` \\[m/F\\], zero when `r_in` is zero or equals
`r_ex`.

- `r_in`: inner radius of the layer \\[m\\].
- `r_ex`: outer radius of the layer \\[m\\].
- `κ`: admittivity of the layer \\[S/m\\].
- `s`: complex frequency, `j 2π f` \\[1/s\\].

`formula(n, r_in, r_ex, κ, s)` returns the record `(; inner, outer, transfer)` of the layer
at the harmonic order `n` \\[m/F\\]. The order `n` \\[dimensionless\\] is a positive
integer, and the first form gives the coefficient at order 0. The radii satisfy
`0 ≤ r_in < r_ex ≤ Inf`, and `r_ex = Inf` gives a medium without an outer bound.
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
