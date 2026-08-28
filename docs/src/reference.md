# API reference

This page lists the public calculation API by its defining module. Within each
module, Documenter orders constants, types, functions and methods, then macros.
Convenience functions and extension hooks have separate references.

## Calculation grammar

```@autodocs
Modules = [
    LineCableModels,
    LineCableModels.Grammar,
    LineCableModels.Units,
]
Order = [:module, :constant, :type, :function, :macro]
Filter = api_reference_entry
Public = true
Private = false
```

## Materials

```@autodocs
Modules = [LineCableModels.Materials]
Order = [:module, :constant, :type, :function, :macro]
Filter = api_reference_entry
Public = true
Private = false
```

## Data model

```@autodocs
Modules = [
    LineCableModels.DataModel,
    LineCableModels.EarthProps,
]
Order = [:module, :constant, :type, :function, :macro]
Filter = api_reference_entry
Public = true
Private = false
```

## Line and cable calculations

```@autodocs
Modules = [
    LineCableModels.Engine,
    LineCableModels.Engine.EarthAdmittance,
    LineCableModels.Engine.EarthProperties,
    LineCableModels.Engine.EarthImpedance,
    LineCableModels.Engine.EHEM,
    LineCableModels.Engine.InsulationAdmittance,
    LineCableModels.Engine.InsulationImpedance,
    LineCableModels.Engine.InternalImpedance,
    LineCableModels.Engine.Transforms,
]
Order = [:module, :constant, :type, :function, :macro]
Filter = api_reference_entry
Public = true
Private = false
```

## Parameter spaces

```@autodocs
Modules = [
    LineCableModels.ParametricBuilder,
    LineCableModels.ParametricBuilder.Conductor,
    LineCableModels.ParametricBuilder.Insulator,
    LineCableModels.ParametricBuilder.Semiconductor,
    LineCableModels.ParametricBuilder.Filler,
]
Order = [:module, :constant, :type, :function, :macro]
Filter = api_reference_entry
Public = true
Private = false
```

## Uncertainty quantification

```@autodocs
Modules = [LineCableModels.UQ]
Order = [:module, :constant, :type, :function, :macro]
Filter = api_reference_entry
Public = true
Private = false
```

## Import and export

```@autodocs
Modules = [LineCableModels.ImportExport]
Order = [:module, :constant, :type, :function, :macro]
Filter = api_reference_entry
Public = true
Private = false
```

## Plots, reports, and tables

```@docs
LineCableModels.PlotBuilder.UIPlot
LineCableModels.export_svg
LineCableModels.set_backend!
LineCableModels.ReportBuilder.ReportArtifact
LineCableModels.ReportBuilder.TableReportDefinition
LineCableModels.ReportBuilder.XLSXReportDefinition
LineCableModels.report
```

## Index

```@index
Pages = ["reference.md"]
Order = [:module, :constant, :type, :function, :macro]
```
