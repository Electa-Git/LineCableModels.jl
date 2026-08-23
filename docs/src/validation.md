# Data entry validation

## Contents

```@contents
Pages = ["validation.md"]
Depth = 3
```

[`validate`](@ref) is the common validation entry point. It returns its argument
unchanged when the value is valid and throws a native Julia exception when it is not:

```julia
validated = validate(value)
@assert validated === value
```

Validation never converts, fills defaults, rewrites nested values, or mutates its input.
Constructors first promote their complete input with Julia's ordinary `promote` and
`convert` mechanisms, then validate the resulting immutable record. Mutable operations
calculate and validate a complete candidate state before changing the owned collection.

## Declarative rules

Types with simple field constraints specialize
`LineCableModels.Validation.rules(::Type)`. The generic `Validation.check` evaluates
those rules and returns the original value:

```julia
import LineCableModels: validate
import LineCableModels.Validation: Positive, Less, check, rules

struct Annulus{T <: Real}
    r_in::T
    r_ex::T
end

rules(::Type{<:Annulus}) = (Positive(:r_ex), Less(:r_in, :r_ex))
validate(value::Annulus) = check(typeof(value), value)
```

Available primitive rules cover finite, positive, nonnegative, integer, type,
membership, field-ordering, and cable-packing checks. Cross-field rules that require
domain knowledge remain beside their owning type. A direct `validate(::OwnedType)`
method is appropriate when dispatch expresses the check more clearly than a new generic
rule.

The failure type communicates the category:

- `ArgumentError` for invalid choices or value kinds;
- `DomainError` for values outside a physical or mathematical domain;
- `DimensionMismatch` for incompatible shapes;
- `MethodError` when no supported operation exists for the supplied types;
- `KeyError` for missing library entries.

## Materialized constructors

Cable-part constructors accept resolved numeric geometry. For example:

```julia
part = Tubular(r_in, r_ex, material)
```

Radius-versus-thickness intent, repetition, and variation belong to the builder
definitions. The definitions promote all intent before materializing the first cable object. Materialized
objects therefore contain one scalar type and one unambiguous geometry.

`add!` on a materialized group, design, earth model, or system accepts only the same
scalar type as the destination. A mixed-scalar insertion throws before mutation; use an
explicit `convert` or build through a Spec when whole-description promotion is wanted.

Operating temperature is not a cable-part constructor input. Cable designs represent
the common material reference state and reject mixed material reference temperatures.
The line problem owns the operating temperature, and [`compute`](@ref) applies its
correction without mutating the design.

## Packing limits

[`maxfill`](@ref) is the single source for cable packing limits:

```julia
maxfill(CircStrands, r_in, wire_radius)
maxfill(RectStrands, lay_radius, width)
```

The cable-part validators and WirePatterns estimators call these same methods, so direct
construction and best-effort estimation use identical geometric limits.

## API reference

```@autodocs
Modules = [LineCableModels.Validation]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Index

```@index
Pages = ["validation.md"]
Order = [:module, :constant, :type, :function, :macro]
```
