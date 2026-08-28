# Conveniences

These functions support model preparation without changing the line-parameter
calculation grammar.

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
