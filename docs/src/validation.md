# Data entry validation

## Contents

```@contents
Pages = ["validation.md"]
Depth = 3
```

[`validate`](@ref) checks materialized inputs before they enter another
construction or numerical operation. It returns its argument unchanged when
the value is usable and throws a native Julia exception when it is not:

```julia
validated = validate(value)
@assert validated === value
```

Input validation never converts, fills defaults, rewrites nested values, or
mutates its input. Constructors normalize admitted scalar grammar before
checking the resulting object. A selected `Gridpoint` is checked when it
materializes, and computation entry points check their complete problem again.

## Validation failures

Construction and mutation check complete values before returning or changing
them. Rechecking at a numerical boundary detects changes to mutable vectors or
dictionaries retained by an otherwise immutable object.

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
`Grid` inputs. Scalar-complete calls invoke the action directly; varying calls
materialize the same action through `Gridspace`. Completed objects therefore
contain one resolved geometry and cannot drift from their declarations.

Mutable libraries validate a complete candidate before changing owned state.
Earth models, cable designs, and line systems are immutable descriptions;
rebuild the authoritative declaration when it changes.

Operating temperature is not a cable-part constructor input. Cable designs represent
the common material reference state and reject mixed material reference temperatures.
The line problem owns the operating temperature, and [`compute`](@ref) applies its
correction without mutating the design.

## Reference

```@docs
LineCableModels.validate
```
