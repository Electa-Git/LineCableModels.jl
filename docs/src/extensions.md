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
