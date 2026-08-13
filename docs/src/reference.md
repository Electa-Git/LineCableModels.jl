# API reference

This page documents the public API and the documented implementation surface of
`LineCableModels.jl`.

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
