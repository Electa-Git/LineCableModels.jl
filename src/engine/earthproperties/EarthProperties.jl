"""
    LineCableModels.Engine.EarthProperties

Compute the resistivity, relative permittivity, and relative permeability of the earth layers,
and their equivalent earth.

A soil law reads the resistivity ``\\rho_\\mathrm{low}`` \\[Ω·m\\] of a layer, its conductivity
``\\sigma_\\mathrm{low} = 1 / \\rho_\\mathrm{low}`` \\[S/m\\] and the frequency ``f`` \\[Hz\\]. It gives
the resistivity ``\\rho(f)`` \\[Ω·m\\], which is ``1 / \\sigma(f)`` unless the law states otherwise,
and the relative permittivity ``\\varepsilon_r(f)`` \\[dimensionless\\]. The described layer keeps
the relative permeability and the thickness of `layer`. The constant ``\\varepsilon_0`` is [`ε₀`](@ref).

# Formulas

## `Formula{:constant}`

Properties of an earth layer independent of the frequency ``f``:

```math
\\rho(f) = \\rho, \\qquad \\varepsilon_r(f) = \\varepsilon_r, \\qquad \\mu_r(f) = \\mu_r,
```

where ``\\rho`` is the resistivity of the layer \\[Ω·m\\], and ``\\varepsilon_r`` and
``\\mu_r`` are its relative permittivity and permeability \\[dimensionless\\].

### Callable

`formula(layer, f)` returns `layer`, whose properties hold at every frequency.

- `layer`: earth layer.
- `f`: frequency \\[Hz\\].

## `Formula{:bottommost}`

Equivalent earth of the layers below the air: the deepest described layer, the identity for
one layer.

### Callable

`formula(layers, f)` returns the last of `layers`.

- `layers`: the layers below the air, as the law describes them at `f`.
- `f`: frequency \\[Hz\\].

## `Formula{:alipio2014}`

Causal soil law of Alipio and Visacro [Alipio2014](@cite), with
``\\sigma_\\mathrm{low} = 1000 / \\rho_\\mathrm{low}`` \\[mS/m\\] and ``D = a \\sigma_\\mathrm{low}^{b}``:

```math
\\begin{aligned}
\\sigma(f) &= \\sigma_\\mathrm{low} \\left[1 + D \\left(f / 10^6\\right)^\\gamma\\right], \\\\
\\varepsilon_r(f) &= \\varepsilon_\\mathrm{high} + \\tan\\left(\\frac{\\pi \\gamma}{2}\\right)
    \\frac{10^{-3} \\sigma_\\mathrm{low} D f^{\\gamma - 1}}{2 \\pi \\varepsilon_0 10^{6 \\gamma}},
\\end{aligned}
```

and ``\\rho(f) = 1000 / \\sigma(f)``. The options are `exponent` ``\\gamma``,
`epsilon_high` ``\\varepsilon_\\mathrm{high}``, `scale` ``a`` and `sigma_exponent` ``b``.
A computation warns outside 100 Hz to 4 MHz.

### Callable

`formula(layer, f)` returns the layer that the law describes at `f`.

- `layer`: earth layer.
- `f`: frequency \\[Hz\\].

## `Formula{:cigre2019}`

Soil law of CIGRE Technical Brochure 781 [CIGRE2019](@cite):

```math
\\begin{aligned}
\\varepsilon_r(f) &= \\varepsilon_\\mathrm{high} + a \\sigma_\\mathrm{low}^{b} f^{c}, \\\\
\\sigma(f) &= \\sigma_\\mathrm{low} + d \\sigma_\\mathrm{low}^{b} f^{e}.
\\end{aligned}
```

The options are `epsilon_high` ``\\varepsilon_\\mathrm{high}``, `epsilon_scale` ``a``,
`epsilon_sigma_exponent` ``b``, `epsilon_frequency_exponent` ``c``, `sigma_scale`
``d`` and `sigma_frequency_exponent` ``e``. A computation warns outside 100 Hz to 4 MHz.

### Callable

`formula(layer, f)` returns the layer that the law describes at `f`.

- `layer`: earth layer.
- `f`: frequency \\[Hz\\].

## `Formula{:datsios2019}`

Two-limit fit of Datsios and Mikropoulos for sandy soil [Datsios2019](@cite), with
``\\sigma_\\mathrm{low} = 10^4 / \\rho_\\mathrm{low}`` \\[μS/cm\\] at 42 Hz and the option `epsilon_dry`
``\\varepsilon_\\mathrm{dry}``:

```math
\\begin{aligned}
p &= 0.537 \\sigma_\\mathrm{low}^{0.16}, \\\\
\\varepsilon_\\mathrm{high} &= 1.24 \\sigma_\\mathrm{low}^{0.415} \\varepsilon_\\mathrm{dry}, \\\\
\\varepsilon_{3k} &= 4 \\sigma_\\mathrm{low}^{0.463} (2.9 \\varepsilon_\\mathrm{dry} - 3.8), \\\\
\\varepsilon_r(f) &= \\varepsilon_\\mathrm{high}
    + \\left(\\frac{3000}{\\max(f, 3000)}\\right)^p (\\varepsilon_{3k} - \\varepsilon_\\mathrm{high}), \\\\
\\sigma_\\mathrm{high} &= \\sigma_\\mathrm{low} (1 + 0.65 \\sigma_\\mathrm{low}^{-0.57}), \\\\
\\sigma(f) &= 10^{-6} \\sigma_\\mathrm{high} f + 42 \\left[1 - 10^{-6} (f - 42)\\right]
    \\left(\\frac{\\sigma_\\mathrm{low}}{42} - 10^{-6} \\sigma_\\mathrm{high}\\right),
\\end{aligned}
```

and ``\\rho(f) = 10^4 / \\sigma(f)``. A computation warns outside 42 Hz to 1 MHz.

### Callable

`formula(layer, f)` returns the layer that the law describes at `f`.

- `layer`: earth layer.
- `f`: frequency \\[Hz\\].

## `Formula{:longmire1975}`

Thirteen-term relaxation law of Longmire and Smith [Longmire1975](@cite), with the corner
frequency ``f_c = (a \\sigma_\\mathrm{low})^{b}`` and ``f_n = 10^{n - 1} f_c``:

```math
\\begin{aligned}
\\varepsilon_r(f) &= \\varepsilon_\\mathrm{high} + \\sum_{n = 1}^{13} \\frac{a_n}{1 + (f / f_n)^2}, \\\\
\\sigma(f) &= \\sigma_\\mathrm{low} + 2 \\pi f \\varepsilon_0
    \\sum_{n = 1}^{13} \\frac{a_n (f / f_n)}{1 + (f / f_n)^2},
\\end{aligned}
```

where ``a_n`` are the coefficients of the report. The options are `epsilon_high`
``\\varepsilon_\\mathrm{high}``, `corner_scale` ``a`` and `corner_exponent` ``b``. A computation warns
outside 5 Hz to 30 GHz.

### Callable

`formula(layer, f)` returns the layer that the law describes at `f`.

- `layer`: earth layer.
- `f`: frequency \\[Hz\\].

## `Formula{:messier1985}`

Square-root soil law of Messier [Messier1985](@cite):

```math
\\begin{aligned}
\\varepsilon_r(f) &= \\varepsilon_\\mathrm{high}
    + \\sqrt{\\frac{\\sigma_\\mathrm{low} \\varepsilon_\\mathrm{high}}{\\pi f \\varepsilon_0}}, \\\\
\\sigma(f) &= \\sigma_\\mathrm{low} + \\sqrt{4 \\pi f \\sigma_\\mathrm{low} \\varepsilon_0 \\varepsilon_\\mathrm{high}}.
\\end{aligned}
```

The option is `epsilon_high` ``\\varepsilon_\\mathrm{high}``. A computation warns outside 100 Hz to
1 MHz [Messier1985, Pascoalato2022](@cite).

### Callable

`formula(layer, f)` returns the layer that the law describes at `f`.

- `layer`: earth layer.
- `f`: frequency \\[Hz\\].

## `Formula{:portela1999}`

Power-law soil dispersion of Portela [Portela1999](@cite), with ``\\omega = 2 \\pi f``:

```math
\\begin{aligned}
\\sigma(f) &= \\sigma_\\mathrm{low} + 10^{-6} \\beta \\omega^\\alpha, \\\\
\\varepsilon_r(f) &= 10^{-6} \\beta \\tan\\left(\\frac{\\pi \\alpha}{2}\\right)
    \\frac{\\omega^{\\alpha - 1}}{\\varepsilon_0}.
\\end{aligned}
```

The options are `beta` ``\\beta`` and `exponent` ``\\alpha``. A computation warns outside 100 Hz
to 2 MHz [Portela1999, Pascoalato2022](@cite).

### Callable

`formula(layer, f)` returns the layer that the law describes at `f`.

- `layer`: earth layer.
- `f`: frequency \\[Hz\\].

## `Formula{:scott1967}`

Empirical fit of Scott, Carroll and Cunningham for moist rock [Scott1967](@cite), with
``u = \\log_{10}(1000 / \\rho_\\mathrm{low})`` and ``x = \\log_{10} f``:

```math
\\begin{aligned}
\\log_{10} \\varepsilon_r(f) &= 5.491 + 0.946 u - 1.097 x + 0.069 u^2 - 0.114 u x + 0.067 x^2, \\\\
\\log_{10} \\sigma(f) &= 0.028 + 1.098 u - 0.068 x + 0.036 u^2 - 0.046 u x + 0.018 x^2,
\\end{aligned}
```

and ``\\rho(f) = 1000 / \\sigma(f)``, with ``\\sigma`` in \\[mS/m\\]. A computation warns
outside 100 Hz to 1 MHz.

### Callable

`formula(layer, f)` returns the layer that the law describes at `f`.

- `layer`: earth layer.
- `f`: frequency \\[Hz\\].

## `Formula{:visacro1987}`

Empirical soil law of Visacro and Portela [Visacro1987](@cite):

```math
\\begin{aligned}
\\varepsilon_r(f) &= 2.34 \\cdot 10^6 \\sigma_\\mathrm{low}^{0.535} f^{-0.597}, \\\\
\\sigma(f) &= \\sigma_\\mathrm{low} \\left(f / f_0\\right)^{0.072},
\\end{aligned}
```

where the option `normalization_frequency` is ``f_0`` \\[Hz\\]. A computation warns outside 40 Hz
to 2 MHz [Visacro1987, Pascoalato2021](@cite).

### Callable

`formula(layer, f)` returns the layer that the law describes at `f`.

- `layer`: earth layer.
- `f`: frequency \\[Hz\\].

## `Formula{:visacro2012}`

Empirical soil law of Visacro and Alipio [Visacro2012](@cite), with ``f_e = \\max(f, f_0)``:

```math
\\begin{aligned}
\\varepsilon_r(f) &= 1.3 + 7.6 \\cdot 10^3 f_e^{-0.4}, \\\\
\\sigma(f) &= \\sigma_\\mathrm{low} + 1.2 \\cdot 10^{-6} \\sigma_\\mathrm{low}^{0.27} (f_e - f_0)^{0.65},
\\end{aligned}
```

where the option `frequency_boundary` is ``f_0`` \\[Hz\\]. Below it, the law keeps its value at
``f_0``. A computation warns outside 100 Hz to 4 MHz.

### Callable

`formula(layer, f)` returns the layer that the law describes at `f`.

- `layer`: earth layer.
- `f`: frequency \\[Hz\\].

## `Formula{:martins2020}`

Equivalent conductivity of Martins-Britto, Lopes and Rondineau [MartinsBritto2020](@cite). From
``\\sigma_{\\mathrm{eq},N} = \\sigma_N`` at the deepest layer, the recursion climbs to the layer
below the air, where ``\\sigma_\\mathrm{eq} = \\sigma_{\\mathrm{eq},1}``:

```math
\\sigma_{\\mathrm{eq},k} = \\sigma_k \\left[
    \\frac{\\sqrt{\\sigma_k} + \\sqrt{\\sigma_{\\mathrm{eq},k+1}}
        - \\left(\\sqrt{\\sigma_k} - \\sqrt{\\sigma_{\\mathrm{eq},k+1}}\\right) e^{-2 h_k \\sqrt{\\pi f \\mu_k \\sigma_k}}}
    {\\sqrt{\\sigma_k} + \\sqrt{\\sigma_{\\mathrm{eq},k+1}}
        + \\left(\\sqrt{\\sigma_k} - \\sqrt{\\sigma_{\\mathrm{eq},k+1}}\\right) e^{-2 h_k \\sqrt{\\pi f \\mu_k \\sigma_k}}}
    \\right]^2,
```

where layer ``k`` has the conductivity ``\\sigma_k = 1 / \\rho_k`` \\[S/m\\], the thickness ``h_k``
\\[m\\] and the permeability ``\\mu_k = \\mu_0 \\mu_{r,k}`` \\[H/m\\], and ``\\mu_0`` is [`μ₀`](@ref). The
layer it returns has the resistivity ``1 / \\sigma_\\mathrm{eq}`` and takes the relative
permittivity and the relative permeability of the deepest layer. The source validates the
equivalent on overhead conductors from 1 Hz to 2 MHz.

### Callable

`formula(layers, f)` returns one semi-infinite layer.

- `layers`: the layers below the air, as the law describes them at `f`.
- `f`: frequency \\[Hz\\].

## `Formula{:xue2021}`

Equivalent propagation constant of Xue and coauthors for a nonmagnetic earth [Xue2021](@cite). With
``s = \\mathrm{j} \\omega``, ``\\omega = 2 \\pi f`` and the constants ``\\varepsilon_0`` and
``\\mu_0``, layer ``k`` of conductivity ``\\sigma_k`` \\[S/m\\], relative permittivity
``\\varepsilon_{r,k}`` and thickness ``h_k`` \\[m\\] has the constant ``\\gamma_k`` \\[1/m\\] relative
to the air. From ``\\gamma_{\\mathrm{eq},N} = \\gamma_N`` at the deepest layer, the recursion climbs
to the layer below the air, where ``\\gamma_\\mathrm{eq} = \\gamma_{\\mathrm{eq},1}``:

```math
\\begin{aligned}
\\gamma_k^2 &= s \\mu_0 \\left(\\sigma_k + s \\varepsilon_0 \\varepsilon_{r,k}\\right)
    - s^2 \\mu_0 \\varepsilon_0, \\\\
\\gamma_{\\mathrm{eq},k} &= \\gamma_k
    \\frac{\\gamma_k + \\gamma_{\\mathrm{eq},k+1}
        - \\left(\\gamma_k - \\gamma_{\\mathrm{eq},k+1}\\right) e^{-2 h_k \\gamma_k}}
    {\\gamma_k + \\gamma_{\\mathrm{eq},k+1}
        + \\left(\\gamma_k - \\gamma_{\\mathrm{eq},k+1}\\right) e^{-2 h_k \\gamma_k}}.
\\end{aligned}
```

The layer it returns parametrizes ``\\gamma_\\mathrm{eq}^2`` by its resistivity
``1 / \\sigma_\\mathrm{eq}`` and its relative permittivity:

```math
\\sigma_\\mathrm{eq} = \\frac{\\Im\\left(\\gamma_\\mathrm{eq}^2\\right)}{\\omega \\mu_0}, \\qquad
\\varepsilon_{r,\\mathrm{eq}} = 1 - \\frac{\\Re\\left(\\gamma_\\mathrm{eq}^2\\right)}{\\omega^2 \\mu_0 \\varepsilon_0}.
```

This ``\\varepsilon_{r,\\mathrm{eq}}`` is negative where the layers are inductive. The layer takes
the relative permeability of the deepest layer. A computation warns for an earth layer whose
relative permeability differs from 1. The source reports errors below 4 % from 1 Hz to 5 MHz on
overhead conductors.

### Callable

`formula(layers, f)` returns one semi-infinite layer.

- `layers`: the layers below the air, as the law describes them at `f`.
- `f`: frequency \\[Hz\\].
"""
module EarthProperties

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...Commons: defaults, validate
using ...Commons: FormulationOptions
using ..Engine: AbstractEarthProperties, admittivity, ε₀, μ₀
import ..Engine: domain
using ...Earth: EarthLayer, EarthModel
using IntervalSets: (..)

export Formula, FORMULAS

"Identifiers of the earth-properties formulas."
const FORMULAS = (
    :constant, :bottommost, :alipio2014, :cigre2019, :datsios2019, :longmire1975,
    :messier1985, :portela1999, :scott1967, :visacro1987, :visacro2012, :martins2020, :xue2021)

"""
$(TYPEDEF)

Represent a law of the properties of an earth layer, or a model of the equivalent earth.
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
include("alipio2014.jl")
include("cigre2019.jl")
include("datsios2019.jl")
include("longmire1975.jl")
include("messier1985.jl")
include("portela1999.jl")
include("scott1967.jl")
include("visacro1987.jl")
include("visacro2012.jl")
include("martins2020.jl")
include("xue2021.jl")

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
