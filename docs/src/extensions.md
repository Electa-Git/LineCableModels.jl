# Extension API

Extension code adds methods to the owner that defines their meaning. The public
calculation entry points remain `compute`, `observe`, `observables`, `project`,
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

## Global-sensitivity integration

The core package owns the dependency-neutral [`Sensitivity`](@ref) and
[`SensitivityResult`](@ref) contracts. The GlobalSensitivity extension adds a
narrow `Sensitivity(inner, method::GlobalSensitivity.Sobol, requests; ...)`
convenience and `compute(::ParametricProblem,
::Sensitivity{<:Any,<:GlobalSensitivity.Sobol})`. It retains the dependency's
method and the QuasiMonteCarlo sampler directly; it does not mirror their
registries, define a local Sobol type, or add calculation stages.

An external sensitivity integration follows the same ownership boundary: keep
the method and sampler dependency-owned, specialize `compute` on the concrete
method family, realize physical values through the selected Gridspace point,
and return the core-owned result schema. Native observable requests are the
only output boundary.

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
