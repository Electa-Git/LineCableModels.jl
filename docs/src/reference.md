# API reference

This page lists the public calculation API by its defining module. Within each
module, Documenter orders constants, types, functions and methods, then macros.
Convenience functions and extension interfaces have separate references.

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
Modules = [LineCableModels.Materials, LineCableModels.Materials.TemperatureDependent]
Order = [:module, :constant, :type, :function, :macro]
Filter = api_reference_entry
Public = true
Private = false
```

## Data model

```@autodocs
Modules = [
    LineCableModels.DataModel,
    LineCableModels.Earth,
    LineCableModels.Earth.FrequencyDependent,
    LineCableModels.Earth.EquivalentHomogeneous,
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
    LineCableModels.Engine.EarthImpedance,
    LineCableModels.Engine.InsulationAdmittance,
    LineCableModels.Engine.SemiconAdmittance,
    LineCableModels.Engine.InsulationImpedance,
    LineCableModels.Engine.InternalImpedance,
    LineCableModels.Engine.PipeImpedance,
    LineCableModels.Engine.ShuntModel,
    LineCableModels.ModalAnalysis,
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

## PSCAD backend

```@autodocs
Modules = [LineCableModels.PSCAD]
Order = [:module, :constant, :type, :function, :macro]
Filter = api_reference_entry
Public = true
Private = false
```

## Plots, reports, and tables

```@docs
LineCableModels.UIPlot
LineCableModels.plot
LineCableModels.preview
LineCableModels.show_material_scale
LineCableModels.export_svg
LineCableModels.figurecolorbars!
LineCableModels.axisscale!
LineCableModels.resetview!
LineCableModels.addwidget!
LineCableModels.removewidget!
LineCableModels.figurelegend!
LineCableModels.panellegend!
LineCableModels.figuretitle!
LineCableModels.paneltitle!
LineCableModels.plotwindow
LineCableModels.materialcolors
LineCableModels.materialscale!
LineCableModels.ReportBuilder.ReportArtifact
LineCableModels.ReportBuilder.TableReportDefinition
LineCableModels.ReportBuilder.XLSXReportDefinition
LineCableModels.report
```

### Reporting selections and units

Ordinary reporting accepts `report(source; values=nothing, ...)` or
`report(source, selection; ...)`. Sources include completed cable constants,
line parameters, standalone Z/Y results with `frequencies` context, result
collections, parametric studies, UQ results, and retained observations.

```julia
report(constants)
report(constants; values=(R, L, G, C), length_unit=:kilo,
    quantity_units=(R=:base, L=:milli, G=:micro, C=:micro))
report(line_parameters; values=@observe(R[1, 1, 1:12]))
```

Omitted `values`, `nothing`, and `()` select the source owner's defaults or all
retained products. Complete line results default to R/X and G/B; selecting
R/L/G/C requests those four quantities instead. A single selector such as `R`
produces only that quantity's table. Positional and keyword selections cannot
be supplied together. Reporting uses `values`; plotting uses `ydata`.

`units`, `length_unit`, `quantity_units`, and `frequency_unit` express display
units through the observation owner. `freq_unit` is an alternative spelling of
`frequency_unit`; supplying both is an error. Retained observations preserve
recorded units when these options are omitted. Raw-only `clip`, `atol`, and
`frequencies` cannot be supplied for retained inputs.

A separate atomic `reference` is retained without computing comparisons.
`illustration=true` or a callable explicitly requests a plot; its options belong
in `plot_options`, and its `ydata` selection must agree with the report's `values`.
The illustration receives the prepared observations. Default reporting returns
in-memory tables without loading a plotting backend or writing files.

Explicit snapshots are useful for saving detached scientific products. Report
definitions remain useful for specialized operations: `BenchmarkTableDefinition`
organizes comparisons, and `XLSXReportDefinition` requests file output.
`TableReportDefinition` selects retained products within the definition-based
extension workflow; ordinary quantity tables need only `report(source)`.

### Display and quantity-table access

Ordinary quantity reports display their completed tables, with recorded units
and owner-provided gridpoint descriptions. Plain and HTML output provide bounded
previews; an IO context with `:limit => false` shows all reported rows, columns,
and gridpoints. Compact display remains a single-line summary.

`artifact[R]` returns the stored resistance DataFrame for an atomic report, or a
vector in reported-result order for a collection. `artifact[i, R]` returns the
table for result position `i`. Complete transformation and statistical requests
retain their own identities. Lookup preserves the selection made by `values`,
requires an exact match for indexed requests, and performs no new observation or tabulation.
Returned tables are shared with the report; use `copy` for independent edits.
The `.tables` field remains available for inspection of the underlying grouped
tables and specialized summaries.

```@docs
Base.getindex(::LineCableModels.ReportBuilder.ReportArtifact, ::Any)
```

## Index

```@index
Pages = ["reference.md"]
Order = [:module, :constant, :type, :function, :macro]
```
