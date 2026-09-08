# Extension API

Extension code adds methods to the owner that defines their meaning. The public
calculation entry points remain `compute`, `observe`, `observables`,
`report`, `plot`, and `preview`.

`PlotBuilder` owns only optional plotting entry points and the live `UIPlot`
handle. Scientific owners expose observations, quantities, geometry, and
physical property ranges. The Makie extension owns request normalization,
palettes, legend grouping, layout sugar, widgets, and native rendering.

## Input validation

A materialized input owns one direct `validate(::OwnedType)` method. The method
returns its argument unchanged or throws a native Julia exception identifying
the rejected field and value:

```julia
import LineCableModels: validate

struct Annulus{T <: Real}
    r_in::T
    r_ex::T

    function Annulus{T}(r_in::T, r_ex::T) where {T <: Real}
        return validate(new{T}(r_in, r_ex))
    end
end

function validate(value::Annulus)
    value.r_in >= zero(value.r_in) || throw(DomainError(
        value.r_in,
        "Annulus.r_in must be nonnegative"
    ))
    value.r_ex > value.r_in || throw(DomainError(
        value.r_ex,
        "Annulus.r_ex must be greater than r_in"
    ))
    return value
end
```

Do not delegate the method to private checking helpers. Constructors normalize
their admitted grammar; `validate` only checks the completed value and does not
convert, repair, or mutate it.

```@docs
LineCableModels.InputValidation
```

## Grammar, observations, and units

```@autodocs
Modules = [
    LineCableModels.Grammar,
    LineCableModels.Units,
]
Order = [:module, :constant, :type, :function, :macro]
Filter = developer_reference_entry
Public = true
Private = false
```

## Data model and calculations

```@autodocs
Modules = [
    LineCableModels.DataModel,
    LineCableModels.Engine,
    LineCableModels.ParametricBuilder,
    LineCableModels.ImportExport,
]
Order = [:module, :constant, :type, :function, :macro]
Filter = developer_reference_entry
Public = true
Private = false
```

## Uncertainty realisation extensions

An uncertainty method starts from the unresolved points returned by
`points(space)`. `uncertainties(point)` returns the positional
`UncertainValue` tuple in depth-first, left-to-right realisation order. It does
not descend into ordinary tuple- or array-valued deterministic leaves.

Deterministic dependency extensions may add a narrow `materialize` method for
an uncertainty descriptor. A method that already owns concrete physical values
uses the release branch's staged realisation contract:

```julia
descriptors = LineCableModels.uncertainties(point)
physical_values = map(nominal, descriptors)
arguments = LineCableModels.realize_arguments(point, physical_values)
problem = LineCableModels.realize(point, arguments)
```

Cardinality is checked before a stored builder is called. Nested builders run
before their parents, and owned constructor validation remains authoritative.
Extensions must not deduplicate descriptors by identity, add a registry that
selects realisation behavior, or overload these operations with a competing
coordinate model. The PolyChaos extension uses this contract for tensor nodes
and independent validation points while leaving RNG-backed Monte Carlo
realisation unchanged.

## Report definitions

```@autodocs
Modules = [LineCableModels.ReportBuilder]
Order = [:module, :constant, :type, :function, :macro]
Filter = developer_reference_entry
Public = true
Private = false
```

## Plotting shell

```@docs
LineCableModels.PlotBuilder
```

## Index

```@index
Pages = ["extensions.md"]
Order = [:constant, :type, :function, :macro]
```
