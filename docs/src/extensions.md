# Extension API

Extension code adds methods to the owner that defines their meaning. The public
calculation entry points remain `compute`, `observe`, `observables`,
`report`, `plot`, and `preview`.

## Validation rules

A type with field-local constraints may define `Validation.rules(::Type)` and
use the common checker:

```julia
import LineCableModels: validate
import LineCableModels.Validation: Positive, Less, check, rules

struct Annulus{T <: Real}
    r_in::T
    r_ex::T
end

rules(::Type{<:Annulus}) = (Positive(:r_ex), Less(:r_in, :r_ex))
validate(value::Annulus) = check(typeof(value), value)
```

Use a direct `validate(::OwnedType)` method when the invariant is not a reusable
field rule. Validation returns its input unchanged or throws a native Julia
exception.

```@autodocs
Modules = [LineCableModels.Validation]
Order = [:module, :constant, :type, :function, :macro]
Filter = developer_reference_entry
Public = true
Private = true
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

## Plot definitions

```@autodocs
Modules = [LineCableModels.PlotBuilder]
Order = [:module, :constant, :type, :function, :macro]
Filter = developer_reference_entry
Public = true
Private = false
```

## Report definitions

```@autodocs
Modules = [LineCableModels.ReportBuilder]
Order = [:module, :constant, :type, :function, :macro]
Filter = developer_reference_entry
Public = true
Private = false
```

## Index

```@index
Pages = ["extensions.md"]
Order = [:constant, :type, :function, :macro]
```
