# Data entry validation

## Contents

```@contents
Pages = ["validation.md"]
Depth = 3
```

[`validate`](@ref) is the common validation entry point. The function returns its argument
unchanged when the value is valid and throws a native Julia exception when it is not:

```julia
validated = validate(value)
@assert validated === value
```

Validation never converts, fills defaults, rewrites nested values, or mutates its input.
Constructors first promote their complete input with Julia's ordinary `promote` and
`convert` mechanisms, then validate the resulting immutable record. Mutable operations
calculate and validate a complete candidate state before changing the owned collection.

## Validation failures

Construction and mutation validate complete values before returning or changing
them. Validation does not fill defaults, estimate wire counts, rewrite nested
values, or mutate its input.

The failure type communicates the category:

- `ArgumentError` for invalid choices or value kinds.
- `DomainError` for values outside a physical or mathematical domain.
- `DimensionMismatch` for incompatible shapes.
- `MethodError` when no supported operation exists for the supplied types.
- `KeyError` for missing library entries.

## Materialized constructors

Cable-part constructors accept resolved numeric geometry. For example:

```julia
part = Tubular(r_in, r_ex, material)
```

Radius or thickness selection, repetition, and variation belong to the builder
definitions. A builder promotes its complete input before materialising the
first cable object. Materialised objects therefore contain one scalar type and
one resolved geometry.

`add!` on a materialised group, design, earth model, or system accepts only the same
scalar type as the destination. A mixed-scalar insertion throws before mutation. Use an
explicit `convert` or rebuild the complete description when promotion is
required.

Operating temperature is not a cable-part constructor input. Cable designs represent
the common material reference state and reject mixed material reference temperatures.
The line problem owns the operating temperature, and [`compute`](@ref) applies its
correction without mutating the design.

## Packing limits

[`maxfill`](@ref) calculates cable packing limits:

```julia
maxfill(CircStrands, r_in, wire_radius)
maxfill(RectStrands, lay_radius, width)
```

Cable-part validators and WirePatterns estimators call the same methods, so
direct construction and estimation use identical geometric limits.

## Reference

```@docs
LineCableModels.validate
```
