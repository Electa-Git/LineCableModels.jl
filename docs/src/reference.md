# API reference

This page documents the public API and the documented implementation surface of
`LineCableModels.jl`.

## Result containers

`CableConstants` stores canonical per-metre R/L/C values. `LineParameters`
stores frequency-dependent Z/Y matrices together with their domain and either a
`:per_length` or `:total` basis. Selecting matrix indices without a frequency
index returns the complete frequency response.

```jldoctest
julia> using LineCableModels

julia> f = [50.0, 100.0];

julia> z = reshape(ComplexF64[1 + 2im, 3 + 4im], 1, 1, 2);

julia> y = reshape(ComplexF64[5 + 6im, 7 + 8im], 1, 1, 2);

julia> parameters = LineParameters(z, y, f);

julia> basis(parameters)
:per_length

julia> R(parameters, 1, 1) == [1.0, 3.0]
true

julia> Z(parameters, 1, 1, 2)
3.0 + 4.0im
```

Monte Carlo calculations return `CableConstantsMC` or `LineParametersMC`.
Use `statistics`, `mean`, `std`, and `quantile` for summaries; `samples` and
`trial` for retained joint trials; `distribution` for retained marginal
histograms; and `surrogate` for the covariance-preserving Measurements.jl
representation. `trial` and `rand` deliberately require retained joint samples.

## Contents

```@contents
Pages = ["reference.md"]
Depth = 3
```

## Core utilities

```@autodocs
Modules = [
    LineCableModels,
    LineCableModels.Commons,
    LineCableModels.UnitHandler,
    LineCableModels.Utils,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Data model and materials

```@autodocs
Modules = [
    LineCableModels.Materials,
    LineCableModels.DataModel,
    LineCableModels.DataModel.BaseParams,
    LineCableModels.EarthProps,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Line-parameter engine

```@autodocs
Modules = [
    LineCableModels.Engine,
    LineCableModels.Engine.EarthAdmittance,
    LineCableModels.Engine.EarthImpedance,
    LineCableModels.Engine.EHEM,
    LineCableModels.Engine.FEM,
    LineCableModels.Engine.InsulationAdmittance,
    LineCableModels.Engine.InsulationImpedance,
    LineCableModels.Engine.InternalImpedance,
    LineCableModels.Engine.Transforms,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Parametric and uncertainty modeling

```@autodocs
Modules = [
    LineCableModels.ParametricBuilder,
    LineCableModels.ParametricBuilder.WirePatterns,
    LineCableModels.UQ,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Plot specifications

```@autodocs
Modules = [
    LineCableModels.PlotBuilder,
    LineCableModels.PlotBuilder.BackendHandler,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Import and export

```@autodocs
Modules = [LineCableModels.ImportExport]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Uncertainty-aware Bessel functions

```@autodocs
Modules = [LineCableModels.UncertainBessels]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Index

```@index
Pages = ["reference.md"]
Order = [:module, :constant, :type, :function, :macro]
```
