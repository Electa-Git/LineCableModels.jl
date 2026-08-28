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
part = Region(:core, Disk(radius), material)
```

Radius or thickness selection, repetition, and variation use explicit
`Grid` inputs. Scalar-complete calls construct an eager object; varying calls
materialize the same constructor through `Gridspace`. Eager objects therefore
contain one resolved geometry and cannot drift from their declarations.

Mutable libraries and earth models validate a complete candidate before
changing owned state. Cable parts, cable designs, and line systems are
immutable; rebuild the authoritative declaration when it changes.

Operating temperature is not a cable-part constructor input. Cable designs represent
the common material reference state and reject mixed material reference temperatures.
The line problem owns the operating temperature, and [`compute`](@ref) applies its
correction without mutating the design.

## Reference

```@docs
LineCableModels.validate
```
