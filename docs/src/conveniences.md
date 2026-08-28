# Conveniences

These functions support model preparation without changing the line-parameter
calculation grammar.

## Physical cable declarations

The construction namespaces produce ordinary v1 physical objects:

```julia
core = Conductor.Stranded(
    :core,
    DiskDefinition(1.5e-3),
    copper;
    counts=(1, 6, 12, 18, 24),
    lay=(15.0, 13.5, 12.5, 11.0),
)

screen = Conductor.Wires(
    :screen,
    DiskDefinition(0.5e-3),
    copper;
    n=40,
    r=20e-3,
    lay=LayRatio(10),
)

insulation = Insulator.Shell(:insulation, xlpe; t=8e-3)
semiconductor = Semiconductor.Shell(:outer_screen, semicon; t=1e-3)
```

`Conductor.Stranded` returns a `Stack` of physical strand groups. Its `counts`
and `lay` declarations align with the actual noncentral layers; callers do not
expand them into a layer-count/multiplier convention.
`Conductor.Wires`, `Conductor.Strip`, and `Conductor.Tubular` return one
retained conductor `Group`. `Insulator.Shell` and `Semiconductor.Shell` return
contextual physical regions. These values are assembled with `Stack` and
completed by `build(CableDesign, ...)`.

Use `equivalent(design)` only when a homogeneous `CableDesign` is itself the
requested result. Calculation through `LineParametersProblem` performs
formulation adaptation without requiring this extra call.

## Wire-pattern estimates

[`make_stranded`](@ref) and [`make_screened`](@ref) search deterministic wire
patterns and return a [`WireEstimate`](@ref). An infeasible estimate retains
ranked candidates and states which limits were not met.

```julia
estimate = make_stranded(1000.0)
closest = estimate[:match]
fewest_layers = estimate[:layers]
```

The selectors `:match`, `:layers`, `:wires`, and `:diameter` do not change the
stored candidate order.

## VDE designation parsing

The qualified [`LineCableModels.DataModel.vdeparse`](@ref) function decodes the
supported VDE/DIN 0271 and 0276 designation fields:

```julia
fields = LineCableModels.DataModel.vdeparse("N2XS(FL)2Y 1x630/35 76/132 kV RM")
```

The parser returns a dictionary and retains unparsed compact-token text under
`:unparsed_stub`. Parsing a designation does not construct or complete a cable
model.

## Base cable formulas

`DataModel.BaseParams` exposes scalar resistance, inductance, capacitance,
conductance, geometric-mean, and equivalent-material formulas used during
model preparation.

## Reference

```@docs
LineCableModels.DataModel.vdeparse
Base.get
Base.delete!
```

```@autodocs
Modules = [
    LineCableModels.DataModel.BaseParams,
    LineCableModels.ParametricBuilder.WirePatterns,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```
